#!/Users/yuzhang/anaconda3/bin/python
# Filename:     gro_common.py
# Description:  This is a python module that includes some common functions
#               to deal with .gro files, including reading coordinates, etc.
#               In order to work with other systems, the units in this gro
#               class uses angstrom, kCal for length and energy.
# Date: 06-14-2018  Created
#       09-26-2018  Add get_atomtypes
#       10-01-2018  Add bondtypes, bonds, angletypes, angles
#       10-02-2018  Add load_ff, and modified subsequent load_* function, and
#                   use generator to read in each parameter files. This helps
#                   to load file from multiple files.
#       1-15-2019   Load bond and angle parameters if they were found in .itp.
#                   Added self.mq2lj to solve the conflicts raise when same LJ
#                   paramters are shared by atoms with different charges.

from collections import Counter

class Gro(object):
    def __init__(self):
        '''
        info:   list of parsed line in .gro file, based on gro file format, 
                including: residue number, residue name, atom name, 
                atom number position(x, y, z)
        atomtypes:  dictionary of atomtype and its corresponding index
        mq:     list of atom mass and charge
        mass:   dictionary of atomtypes and their corresponding mass
        lj:     dictionary of atomtypes and their corresponding lj parameters
                the order is [sigma, epsilon]
        mq2lj:  lookup dictionary for atomtypes of mass and charge to find
                corresponding lj atomtypes (because same lj type could have
                different charges)
        bondtypes:  dictionary of bond types
        angletypes: dictionary of angle types
        dihedraltypes:  dictionary of dihedral types
        bonds:  list of bonds (index1, index2)
        angles: list of angles (index1, index2, index3)
        dihedrals:  list of dihedrals (index1, index2, index3, index4)
        idx:    dictionary of indexes for bondtypes, angletypes and 
                dihedraltypes. This is used to record the index while
                reading force field info.
        '''
        self.box = None
        self.natoms = None
        self.nmols = None
        self.info = []
        self.residuals = Counter()
        self.atomtypes = {}
        self.mq = []
        self.mass = {}
        self.lj = {}
        self.mq2lj = {}
        self.bondtypes = {}
        self.angletypes = {}
        self.dihedraltypes = {}
        self.bonds = []
        self.angles = []
        self.dihedrals = []
        self.idx = {
                'atom': 1, 
                'bond': 1, 
                'angle': 1,
                'dihedral': 1
                }
    def load_gro(self, filename):
        '''read *.gro file, save atom info and box size'''
        with open(filename) as myfile:
            for i, line in enumerate(myfile):
                if i == 1:
                    self.natoms = int(line)
                elif i > 1 and i < self.natoms+2:
                    self.info.append(Gro.read_line(line))
                    self.residuals[self.info[-1][1]] += 1
            self.box = [float(_) for _ in line.split()]
            # box vector for any coordinate system
            while len(self.box) < 9:
                self.box.append(0)
        self.check_itp_gro()
        self.get_atomtypes()
    def load_ff(self, filenames, kcal = True):
        '''Load molar mass and charge infomation such as molar mass, atomic 
        charge, bonds, angles, dihedrals, etc.

        This function uses generator to loop files and read in itp files, line
        by line, and check if a session keyword is in the line. If true, jump
        to the function that reads that specific session.
        The benefit of using generator here is that we can process all sessions 
        by just going through each file only once.

        Additional thoughts: could create a new class say FF that just stores
        the force field related methods, dictionaries. In the Gro class, we 
        only need to define self.ff = FF(**args). Something like that. 
        But gonna just leave it for now.
        '''
        sessions = {
                'atoms': self.load_mq,
                'atomtypes': self.load_nonbonds,
                'bondtypes': self.load_bondtypes,
                'bonds': self.load_bonds, 
                'angletypes': self.load_angletypes,
                'angles': self.load_angles, 
                }
        for filename in filenames:
            file_kcal = kcal
            if filename.endswith('.itp') or filename.endswith('.top'):
                file_kcal = False
            with open(filename) as myfile:
                filegen = (line for line in myfile if not (line.strip() == '' or 
                                            line.strip()[0] in [';', '#']))
                for line in filegen:
                    if '[' in line and ']' in line:
                        currentline = line
                        for item in sessions:
                            if item in currentline:
                                try:
                                    currentline = sessions[item](filegen, file_kcal)
                                except:
                                    currentline = sessions[item](filegen)
                                if not currentline:
                                    break
        self.check_itp_gro()
        self.nmols = self.natoms//len(self.mq)
    def load_mq(self, filegen):
        '''Read lines after [ atomtypes ] session and store molar mass and
        atomic charge into self.mass and self.mq.
        The key for self.atomtypes is its corresponding atom name, so each
        atom type will only have one atom name.
        '''
        for line in filegen:
            if '[' in line and ']' in line:
                return line
            temp = line.split()
            # use atom name instead of atom type for self.mass
            if temp[4] not in self.mass:
                self.mass[temp[4]] = float(temp[-1])
            self.mq.append([float(_) for _ in temp[-2:]])
            if temp[4] not in self.mq2lj:
                self.mq2lj[temp[4]] = temp[1]
            elif temp[1] != self.mq2lj[temp[4]]:
                raise ValueError('Different LJ atomtypes found!')

    def load_bondtypes(self, filegen, kcal = True):
        '''Read lines after [ bondtypes ]  session and store bondtypes as a 
        dictionary, where the key is a tuple (bondtype1, bondtype2) and 
        the corresponding value is the a list of [index, k, b]
        Note that each bond connect by 2 different atomtypes will have 2 keys
        that point to the same list.
        '''
        for line in filegen:
            if '[' in line and ']' in line:
                return line
            self._load_bondtypes(line, kcal)
            temp = line.split()
            k = float(temp[-1])
            b = float(temp[-2])
            if not kcal:
                k /= 100*4.184*2 # kJ to kCal/A^2
                b *= 10 # nm to angstrom
            self.bondtypes[(temp[0], temp[1])] = [self.idx['bond'], k, b]
            self.bondtypes[(temp[1], temp[0])] = self.bondtypes[(temp[0], temp[1])]
            self.idx['bond'] += 1

    def load_bonds(self, filegen, kcal=True):
        '''Note this only returns the relative index within the molecule,
        it may not match the index in gro file.
        If the bonding parameters are found in these session, add them to
        self.bondtypes immediately.
        '''
        print('loading bonds ...')
        for line in filegen:
            if '[' in line and ']' in line:
                return line
            temp = line.split()
            self.bonds.append([int(_) for _ in temp[:2]])
            if len(temp) == 5:
                atom1 = self.info[self.bonds[-1][0]-1][2]
                atom2 = self.info[self.bonds[-1][1]-1][2]
                # default gromacs units are kJ/mol and nm
                # convert them to kCal and angstrom
                if (atom1, atom2) not in self.bondtypes.keys():
                    self.bondtypes[(atom1, atom2)] = [self.idx['bond'], 
                                                    float(temp[-1])/(100*4.184*2), 
                                                    10*float(temp[-2])]
                    self.bondtypes[(atom2, atom1)] = self.bondtypes[(atom1, atom2)]
                    self.idx['bond'] +=1
                    
    def load_angletypes(self, filegen, kcal = True):
        '''Read lines after [ angletypes ]  session and store bondtypes as a 
        dictionary, where the key is a tuple (bondtype1, bondtype2) and 
        the corresponding value is the a list of [index, k, b]
        Note that each bond connect by 2 different atomtypes will have 2 keys
        that point to the same list.
        '''
        for line in filegen:
            if '[' in line and ']' in line:
                return line
            temp = line.split()
            theta = float(temp[-2]) 
            k = float(temp[-1])
            if not kcal:
                k /= 4.184*2 # kJ/rad^2 to kCal/rad^2
            self.angletypes[(temp[0], temp[1], temp[2])] = [self.idx['angle'], k, theta]
            self.angletypes[(temp[2], temp[0], temp[0])] = self.angletypes[(temp[0], temp[1], temp[2])]
            self.idx['angle'] += 1

    def load_angles(self, filegen):
        '''Note this only returns the relative index within the molecule,
        it may not match the index in gro file.
        If the parameters for angle are found in this session, add them 
        to self.angletypes.
        '''
        print('Loading angles ...')
        for line in filegen:
            if '[' in line and ']' in line:
                return line
            temp = line.split()
            self.angles.append([int(_) for _ in temp[:3]])
            if len(temp) ==  6:
                angle1 = self.info[self.angles[-1][0]-1][2]
                angle2 = self.info[self.angles[-1][1]-1][2]
                angle3 = self.info[self.angles[-1][2]-1][2]
                if (angle1, angle2, angle3) not in self.angletypes.keys():
                    self.angletypes[(angle1, angle2, angle3)] = [self.idx['angle'],
                                                                float(temp[-1])/(4.184*2),
                                                                float(temp[-2])]
                    self.angletypes[(angle3, angle2, angle2)] = self.angletypes[(angle1, angle2, angle3)]
                    self.idx['angle'] += 1

    def load_dihedraltypes(self):
        raise NotImplementedError()

    def load_dihedrals(self):
        raise NotImplementedError()

    def load_nonbonds(self, filegen, kcal = True):
        '''Read gromacs topology file and find the Lennard-Jones parameter
        and write it to self.lj'''
        ljtypes = set(self.mq2lj.values())
        for line in filegen:
            if '[' in line and ']' in line:
                return line
            temp = line.split()
            if temp[0] in ljtypes:
                sigma = float(temp[-2])
                epsilon = float(temp[-1])
                if not kcal:
                    sigma *= 10
                    epsilon /= 4.184
                self.lj[temp[0]] = [epsilon, sigma]

    def get_atomtypes(self):
        for item in self.info:
            if item[2] not in self.atomtypes:
                self.atomtypes[item[2]] = self.idx['atom']
                self.idx['atom'] += 1

    def check_itp_gro(self):
        '''Check if the number of atoms in gro file is a multiple of
        the number of atoms in itp file'''
        if self.natoms and self.mq and self.natoms%len(self.mq):
            raise ValueError('Number of atoms in gro file does not match itp file!')
    
    def write(self, filename, filetype = 'gro'):
        raise NotImplementedError()


    def _write_line(self, infolist, f):
        '''Write a line of atom info into gro file format
        infolist: a list of atom info
        f: file to be written 
        '''
        f.write("%5s%5s%5s%5s%8.3f%8.3f%8.3f\n" %tuple(infolist))

    def split_residuals(self):
        '''Split the gro files into several gro files according to their residual name'''
        file_dict = {}
        for item in self.info:
            if item[1] not in file_dict:
                f = open(item[1]+'.gro', 'w')
                print('Writing to file %s' %(item[1]+'.gro'))
                file_dict[item[1]] = f
                f.write('Residual %s\n' %item[1])
                f.write('%5d\n' %(self.residuals[item[1]]))
            self._write_line(item, file_dict[item[1]])
        for f in file_dict.values():
            f.write('%10.5f%10.5f%10.5f\n' %(self.box[0], self.box[1], self.box[2]))
            f.close()


    @staticmethod
    def read_line(string):
        '''gro file format '%5d%5s%5s%5d%8.3f%8.3f%8.3f'
        '''
        results = []
        for i in range(4):
            results.append(string[i*5:(i+1)*5].strip())
        results.extend([float(_) for _ in string[20:44].split()])
        return results
