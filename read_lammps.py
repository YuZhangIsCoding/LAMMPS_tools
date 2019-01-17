#!/Users/yuzhang/anaconda3/bin/python
# Filename:     read_lammps.py
# Description:  This is a python script that can be used to exact physical
#               infomations from lammps simulation, such as charge, atoms
#               coordinates, etc.
# Dates:        09-29-2018  Created
#               01-16-2018  Record charges for multiple atoms

import argparse
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input', dest = 'input', help = 'input file')
parser.add_argument('-o', '--output', dest = 'output', default = 'read_lammps_out.txt', help = 'output file')
parser.add_argument('-c', '--choice', choices = ('charge', 'xyz'),
        help = 'select the properties: charge or xyz coordinates')

try:
    __IPYTHON__
    args = parser.parse_args([])
except NameError:
    args = parser.parse_args()

def read_charge(args):
    time_pre = None
    time_flag = False
    charge_flag = False
    charge = []
    with open(args.input) as myfile:
        for line in myfile:
            if 'TIMESTEP' in line:
                time_flag = True # flag to record current time
                charge_flag = False # flag to record current charge
                continue
            if time_flag:
                time = int(line)
                if time == time_pre:
                    # make sure same timestep won't be record twice
                    time = None
                else:
                    time_pre = time
                    charge.append([time])
                time_flag = False
                continue
            if time is not None and 'ITEM: ATOMS id type q' in line:
                # after this line, the charge will be appended to the list
                charge_flag = True
                continue
            if charge_flag:
                charge[-1].append(float(line.split()[-1]))
    charge = np.array(charge)
    np.savetxt(args.output, charge, fmt = ['%10d']+['%12.3f']*(charge.shape[1]-1))

if args.choice == 'charge':
    read_charge(args)
else:
    raise NotImplementedError()
