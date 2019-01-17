#!/Users/yuzhang/anaconda3/bin/python

import argparse
from gro_common import Gro

parser = argparse.ArgumentParser()
parser.add_argument('-g', '--gro', help = 'grofile')

try:
    __IPYTHON__
    args = parser.parse_args([])
except NameError:
    args = parser.parse_args()

gro = Gro()
gro.load_gro(args.gro)
gro.split_residuals()
