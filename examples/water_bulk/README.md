# Build bulk simulation of water

This repo gives an example of building Lammps input files from Gromacs-like files for bulk simulations of water.

Prerequisites:

* LAMMPS
* Packmol
* GROMACS (optional)

Files needed:

* in.water: Lammps input file
* water.pdb: pdb file of water that stores atomic positions
* pack_water.inp: input file for packmol to pack water molecules in a simulation box
* water.itp: Gromacs itp file for water, which contains infomation such as atomic mass, charge, bonds, angles, etc.
* supplement.txt: additional force field parameters that stored in gromacs-like format, such as atomtypes (LJ parameters), bondtypes, angletypes, etc.

Steps:

1. Use Packmol to pack water molecules into a simulation box. This will generate bulk.pdb that contains the coordinates of 800 water molecules in a 3-nm-long simulation box.

    `packmol < pack_water.inp`

2. (Optional) Use Gromacs's editconf command to convert bulk.pdb file to gromacs format.

    `gmx_mpi editconf -f bulk.pdb -o bulk.gro -box 3`

    You can use other tools to achieve this, such as VMD.

3. Use gro2lammps.py to generate input files for Lammps.

    `python gro2lammps.py -g bulk.gro -ff water.itp supplement.txt -o data.water -ljout lj.water -kcal`

    To checkout the usuage of gro2lammps.py, please use 

    `python gro2lammps.py -h`

    This step will generate 2 output files: 
    
    * data.water: contains the coordinates and bonding infomation
    * lj.water: contains the Lennard-Jones parameters, which could be written in the lammps input file (in.water).

4. Run Lammps simulation:

    Serial:

    `lmp_serial < in.water`

    Mpi:

    `mpirun -np 16 lmp_mpi < in.water`
