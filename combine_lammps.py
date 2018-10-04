#!/Users/yuzhang/anaconda3/bin/python
# Filename: gro2lammps.py
# Description:  This is a python script that combines lammps data files
# Date:     10-03-2018, created

import argparse
from lammps_common import LAMMPS

parser = argparse.ArgumentParser()
parser.add_argument('-d', '--datafiles', nargs = '*', help = 'lammps data file')
parser.add_argument('-lj', nargs = '*', help = 'lennard-jones parameters')
parser.add_argument('-o', '--output', default = 'data.lammps', help = 'output file name')
parser.add_argument('-ljout', '--ljout', default = 'data.lj', help = 'output file for lj parameters')

try:
    __IPYTHON__
    args = parser.parse_args([])
except NameError:
    args = parser.parse_args()

lmp = LAMMPS()
for i in range(len(args.datafiles)):
    temp = LAMMPS()
    temp.read_data([args.datafiles[i], args.lj[i]])
    lmp.add(temp)
lmp.write_data(args.output, args.ljout)
