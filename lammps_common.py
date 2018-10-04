#!/Users/yuzhang/anaconda3/bin/python
# Filename:     lammps_common.py
# Description:  This is a python script that has a class for lammps data file
#               that stores info for atom, bond, angle, etc.
# Date: 10-02-2018  Created
#       10-03-2018  Recursively read in data files using _read_currentline, 
#                   add check_consistency function to double check

class LAMMPS(object):
    def __init__(self):
        numnames = ['atoms', 'bonds', 'angles', 'dihedrals',
                'atom types', 'bond types', 'angle types', 'dihedral types']
        self.nums = {item: 0 for item in numnames}
        fftypenames = ['Masses', 'Bond Coeffs', 'Angle Coeffs', 'Dihedral Coeffs']
        self.fftypes = {item: [] for item in fftypenames}
        self.atoms= []
        conn_names = ['Bonds', 'Angles', 'Dihedrals']
        self.connections = {item: [] for item in conn_names}
        self.box = []
        self.lj = []
        self.nmols = 0

    def add(self, other):
        trans = {'Masses': 'atom types',
                'Bond Coeffs': 'bond types',
                'Angle Coeffs': 'angle types',
                'Dihedral Coeffs': 'dihedral types',
                }
        base = self.nums.copy()
        for name in self.nums:
            self.nums[name] += other.nums[name]
        for name in other.fftypes:
            for item in other.fftypes[name]:
                temp = item.copy()
                temp[0] += base[trans[name]]
                self.fftypes[name].append(temp)
        for atom in other.atoms:
            temp = atom.copy()
            temp[0] += base['atoms']
            temp[1] += self.nmols
            temp[2] += base['atom types']
            self.atoms.append(temp)
        self.nmols += other.nmols
        for name in self.connections:
            for item in other.connections[name]:
                temp = item.copy()
                temp[0] += base[name.lower()]
                temp[1] += base[name.lower()[:-1]+' types']
                for i in range(2, len(temp)):
                    temp[i] += base['atoms']
                self.connections[name].append(temp)
        for lj in other.lj:
            temp = lj.copy()
            temp[0] += base['atom types']
            temp[1] += base['atom types']
            self.lj.append(temp)
        if self.box:
            for i in range(3):
                self.box[i][0] = min(self.box[i][0], other.box[i][0])
                self.box[i][1] = max(self.box[i][1], other.box[i][1])
        else:
            self.box = other.box.copy()

    def read_data(self, filenames):
        for filename in filenames:
            with open(filename) as myfile:
                filegen = (line for line in myfile if not (line.strip() == '' or
                                line.strip()[0] is '#'))
                for line in filegen:
                    self._read_currentline(line, filegen)
        self.check_consistency()
                
    def _read_currentline(self, currentline, filegen):
        temp = currentline.split()
        if temp[-1] in self.nums:
            self.nums[temp[-1]] = int(temp[0])
        elif temp[-1] == 'types':
            self.nums[' '.join(temp[-2:])] = int(temp[0])
        elif temp[-1] in ['xhi', 'yhi', 'zhi']:
            self.box.append([float(_) for _ in temp[:2]])
        elif temp[-1] in ['Masses', 'Coeffs']:
            typename = ' '.join(temp)
            self._read_types(typename, filegen)
        elif temp[-1] == 'Atoms':
            self._read_atoms(filegen)
        elif temp[-1] in self.connections:
            self._read_connections(temp[-1], filegen)
        elif temp[0] == 'pair_coeff':
            self.lj.append([int(_) for _ in temp[1:3]])
            self.lj[-1].extend([float(_) for _ in temp[3:]])

    def _read_types(self, typename, filegen):
        for line in filegen:
            temp = line.split()
            try:
                self.fftypes[typename].append([int(temp[0])])
                self.fftypes[typename][-1].extend([float(_) for _ in temp[1:]])
            except:
                self._read_currentline(line, filegen)
                return

    def _read_atoms(self, filegen):
        for line in filegen:
            temp = line.split()
            try:
                self.atoms.append([int(_) for _ in temp[:3]])
                self.atoms[-1].extend([float(_) for _ in temp[3:]])
                self.nmols = max(self.nmols, int(temp[1]))
            except:
                self._read_currentline(line, filegen)
                return

    def _read_connections(self, typename, filegen):
        for line in filegen:
            temp = line.split()
            try:
                self.connections[typename].append([int(_) for _ in temp])
            except:
                self._read_currentline(line, filegen)
                return

    def check_consistency(self):
        if self.nums['atoms'] != len(self.atoms):
            raise ValueError('%d atoms do not match the number of coordinates: %d' 
                    %(self.nums['atom'], len(self.atoms)))
        if self.nums['bonds'] != len(self.connections['Bonds']):
            raise ValueError('%d bonds do not match the number of bonds in connections: %d'
                    %(self.nums['bonds'], len(self.connections['Bonds'])))
        if self.nums['angles'] != len(self.connections['Angles']):
            raise ValueError('%d angles do not match the number of angles in connections: %d'
                    %(self.nums['angles'], len(self.connections['Angles'])))
        if self.nums['dihedrals'] != len(self.connections['Dihedrals']):
            raise ValueError('%d dihedrals do not match the number of dihedrals in connections: %d'
                    %(self.nums['dihedrals'], len(self.connections['Dihedrals'])))
        if self.nums['atom types'] != len(self.fftypes['Masses']):
            raise ValueError('%d atom types do not match the number of Masses: %d'
                    %(self.nums['atom types'], len(self.fftypes['Masses'])))
        if self.nums['bond types'] != len(self.fftypes['Bond Coeffs']):
            raise ValueError('%d bond types do not match the number of Bond Coeffs: %d'
                    %(self.nums['bond types'], len(self.fftypes['Bond Coeffs'])))
        if self.nums['angle types'] != len(self.fftypes['Angle Coeffs']):
            raise ValueError('%d angle types do not match the number of Angle Coeffs: %d'
                    %(self.nums['angle types'], len(self.fftypes['Angle Coeffs'])))
        if self.nums['dihedral types'] != len(self.fftypes['Dihedral Coeffs']):
            raise ValueError('%d dihedral types do not match the number of Dihedral Coeffs: %d'
                    %(self.nums['dihedral types'], len(self.fftypes['Dihedral Coeffs'])))
    
    def write_data(self, datafile, ljfile):
        with open(datafile, 'w') as myfile:
            myfile.write('Lammps data file\n\n')
            for name in self.nums:
                if self.nums[name]:
                    myfile.write('%10d %s\n' %(self.nums[name], name))
            myfile.write('\n')
            if self.box:
                for i, direction in enumerate('xyz'):
                    myfile.write('%12.6f%12.6f %s %s\n' 
                            %tuple(self.box[i]+[direction+'lo', direction+'hi']))
            else:
                for i, direction in enumerate('xyz'):
                    myfile.write('%12.6f%12.6f %s %s\n' %(0, 0))

            for name in self.fftypes:
                if self.fftypes[name]:
                    myfile.write('\n%s\n\n' %name)
                    fmt = '{:5d}'+'{:12.6f}'*(len(self.fftypes[name][0])-1)+'\n'
                    for item in self.fftypes[name]:
                        myfile.write(fmt.format(*item))

            myfile.write('\nAtoms\n\n')
            for item in self.atoms:
                myfile.write('%5d%5d%5d%12.5f%8.3f%8.3f%8.3f\n' %tuple(item))

            for name in self.connections:
                if self.connections[name]:
                    myfile.write('\n%s\n\n' %name)
                    fmt = '{:10d}'*len(self.connections[name][0])+'\n'
                    for item in self.connections[name]:
                        myfile.write(fmt.format(*item))
        with open(ljfile, 'w') as myfile:
            for item in self.lj:
                myfile.write('pair_coeff\t%d\t%d%12.6f%12.6f\n' %tuple(item))

def test():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', nargs = '*', help = 'input file')
    args = parser.parse_args()
    if args.input:
        lmp = LAMMPS()
        lmp.read_data(args.input)
        print(lmp.nums)
        print(lmp.fftypes)
        print('Atoms:', len(lmp.atoms))
        for item in lmp.connections:
            print(item, len(lmp.connections[item]))
        print('LJ:',lmp.lj)
        print('Box sizes:', lmp.box)
if __name__ == '__main__':
    test()
