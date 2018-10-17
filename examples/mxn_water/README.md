# Simulation of freezed MXene in water

This repo gives an example of building Lammps input files from Gromacs-like files for two MXene electrodes in pure water.

Prerequisites:

* LAMMPS
* Packmol
* GROMACS (optional)

Files needed:

* in.combined: Lammps input file
* water.pdb: pdb file of water that stores atomic positions.
* pack_water.inp: input file for packmol to pack water molecules around the electrodes.
* water.itp: Gromacs itp file for water, which contains infomation such as atomic mass, charge, bonds, angles, etc.
* supplement.txt: additional force field parameters that stored in gromacs-like format, such as atomtypes (LJ parameters), bon    dtypes, angletypes, etc.
* MXene.gro: Gromacs gro file for MXene
* MXene.itp: Gromacs itp file for MXene
* MXene_nb.itp: supplemental force field parameters for MXene, which contains LJ parameters for each atomtype

Steps:   

1. Use packmol to generate water phase around MXene, and convert the formats to Lammps input formats. Check the water_bulk example for more details.

    ```packmol < pack_water.inp
    gmx_mpi editconf -f bulk.pdb -o bulk.gro
    python gro2lammps.py -g bulk.gro -ff water.itp supplement.txt -o data.water -ljout lj.water -kcal
    ```

    The commands above will generate two files that are in Lammps format:

    * data.water: contains the coordinates and bonding infomation.
    * lj.water: contains the Lennard-Jones parameters.

2. Use gro2lammps.py to convert MXene files to Lammps format.

    `python gro2lammps.py -g MXene.gro -ff MXene.itp MXene_nb.txt -o data.MXene -ljout lj.MXene`

    This will generate 2 MXene files in Lammps format:

    * data.MXene: coordinates of MXene
    * lj.MXene: LJ parameters of MXene

3. Use combine_lammps.gro to combine the data files and LJ files.

    `python combine_lammps.py -d data.MXene data.water -lj lj.MXene lj.water -o data.combined -ljout lj.combined`

    To check out the usage of combine_lammps.gro, please use

    `python combine_lammps.gro -h`

    Similarly, this will generate 2 files for the hybrid system:

    * data.combined: coordinates and bonding info.
    * lj.MXene: all LJ parameters, which could be copied into the input file (in.combined).

4. Run Lammps simulation

    Serial:

    `lmp_serial < in.combined`

    Mpi:

    `mpirun -np 16 lmp_mpi < in.combined`
