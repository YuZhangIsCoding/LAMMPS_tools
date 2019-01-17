#!/Users/yuzhang/anaconda3/bin/python
# Filename: gro2lammps.py
# Description:  This is a python script that converts the gromacs file to
#               LAMMPS file
# Date:     09-25-2018, created

import argparse
from gro_common import Gro
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument('-g', '--gro', help = 'grofile')
parser.add_argument('-ff', '--forcefields', nargs = '*', type = str, 
                    help = 'force field and supplementary files')
parser.add_argument('-o', '--output', default = 'data.lammps', help = 'output file')
parser.add_argument('-ljout', '--ljout', default = 'data.lj', 
                    help = 'output file for lj parameters')
parser.add_argument('-kcal', '--kcal', action = 'store_true', help = 'Units for energy')
try:
    __IPYTHON__
    args = parser.parse_args([])
except NameError:
    args = parser.parse_args()

if args.kcal:
    print('Note: reading additional inputs as kCal, angstrom, etc.')
else:
    print('Note: reading additional inputs as kJ, nm, etc.')

####
# for class gro
gro = Gro()
gro.load_gro(args.gro)
gro.load_ff(args.forcefields, args.kcal)
####

with open(args.output, 'w') as outfile:
    outfile.write('Generate lammps data file from gromacs file %s and force filed files (%s)\n\n' %(args.gro, ', '.join(args.forcefields)))

    outfile.write('%10d atoms\n' %(gro.natoms))
    if gro.bonds:
        outfile.write('%10d bonds\n' %(len(gro.bonds)*gro.nmols))
    if gro.angles:
        outfile.write('%10d angles\n' %(len(gro.angles)*gro.nmols))

    outfile.write('\n')
    outfile.write('%10d atom types\n' %(len(gro.atomtypes)))
    if gro.bondtypes:
        outfile.write('%10d bond types\n' %(len(gro.bondtypes)//2))
    if gro.angletypes:
        outfile.write('%10d angle types\n' %(len(gro.angletypes)//2))
    outfile.write('\n')

    for i, direction in enumerate('xyz'):
        outfile.write('%12.6f%12.6f %s %s\n' %(0, 10*gro.box[i], direction+'lo', direction+'hi'))

    outfile.write('\nMasses\n\n')
    for atom in gro.atomtypes:
        outfile.write('%5d%12.6f\n' %(gro.atomtypes[atom], gro.mass[atom]))

    if gro.bondtypes:
        outfile.write('\nBond Coeffs\n\n')
        appeared = set()
        for bond in gro.bondtypes:
            if gro.bondtypes[bond][0] not in appeared:
                appeared.add(gro.bondtypes[bond][0])
                outfile.write('%5d%8.3f%8.3f\n' %(gro.bondtypes[bond][0], gro.bondtypes[bond][1], gro.bondtypes[bond][2]))

    if gro.angletypes:
        outfile.write('\nAngle Coeffs\n\n')
        appeared = set()
        for angle in gro.angletypes:
            if gro.angletypes[angle][0] not in appeared:
                appeared.add(gro.angletypes[angle][0])
                outfile.write('%5d%8.3f%8.3f\n' %(gro.angletypes[angle][0], gro.angletypes[angle][1], gro.angletypes[angle][2]))

    outfile.write('\nAtoms\n\n')
    for i, item in enumerate(gro.info):
        chargeinx = i%len(gro.mq)
        outfile.write('%8d%8d%5d%12.5f%8.3f%8.3f%8.3f\n'
                %(i+1, i//len(gro.mq)+1, gro.atomtypes[item[2]],
                    gro.mq[i%len(gro.mq)][0], item[-3]*10, item[-2]*10, item[-1]*10))

    if gro.bonds:
        outfile.write('\nBonds\n\n')
        num_per_mol = gro.natoms//gro.nmols
        cnt = 1
        for bond in gro.bonds:
            for i in range(gro.nmols):
                seq = gro.bondtypes[tuple([gro.info[_-1][2] for _ in bond])][0]
                outfile.write('{:10d}{:10d}{:10d}{:10d}\n'.format(cnt, seq, *[_+i*num_per_mol for _ in bond]))
                cnt += 1

    if gro.angles:
        outfile.write('\nAngles\n\n')
        cnt = 1
        for angle in gro.angles:
            for i in range(gro.nmols):
                seq = gro.angletypes[tuple([gro.info[_-1][2] for _ in angle])][0]
                outfile.write('{:10d}{:10d}{:10d}{:10d}{:10d}\n'.format(cnt, seq, *[_+i*num_per_mol for _ in angle]))
                cnt += 1

                
## add pair_coeffs
with open(args.ljout, 'w') as outfile:
    outfile.write('#Pair Coeffs\n\n')
    typelist = list(gro.atomtypes.keys())
    for i in range(len(typelist)):
        for j in range(i, len(typelist)):
            atom1 = gro.mq2lj[typelist[i]]
            atom2 = gro.mq2lj[typelist[j]]
            epsilon = np.sqrt(gro.lj[atom1][0]*gro.lj[atom2][0])
            sigma = (gro.lj[atom1][1]+gro.lj[atom2][1])/2
            outfile.write('pair_coeff\t%d\t%d%12.6f%12.6f\n' 
                    %(gro.atomtypes[typelist[i]], gro.atomtypes[typelist[j]], epsilon, sigma))
