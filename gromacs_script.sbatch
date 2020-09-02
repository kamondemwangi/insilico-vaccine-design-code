#!/usr/bin/env bash
#SBATCH -p batch
#SBATCH -J gromacs
#SBATCH -n 4

#Load modules
module load gromacs/2018.2

#Input = pdb file  # Converting protein.pdb to gromacs input. Force field : OPLS, water type: spce
# hydrogen bonds in the pdb file are ignored.

echo -e "5\n1" | gmx pdb2gmx -f 679.pdb -o protein_processed.gro -water spce -ignh

# Adding simulation box (cubic, min dist from wall 1.0 nm)
gmx editconf -f protein_processed.gro -o protein_newbox.gro -c -d 1.0 -bt cubic

# Adding water mols (water type: spc216)
gmx solvate -cp protein_newbox.gro -cs spc216.gro -o protein_solv.gro -p topol.top

# Adding ions (salt concentration = 0.1)
gmx grompp -f ions.mdp -c protein_solv.gro -p topol.top -o ions.tpr -maxwarn 1


echo -e "13\n1" | gmx genion -s ions.tpr -o protein_solv_ions.gro -p topol.top -neutral

# Energy minimization of protein + water (min.mdp -> boundary = xyz)
gmx grompp -f minim.mdp -c protein_solv_ions.gro -p topol.top -o em.tpr

gmx mdrun -v -deffnm em -nt 4

echo "10" | gmx energy -f em.edr -o potential.xvg

# NVT Temperature control
gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr

gmx mdrun -v -deffnm nvt -nt 4

# NPT Pressure control
gmx grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr

gmx mdrun -v -deffnm npt -nt 4

# Real MD 
gmx grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -o md_0_1.tpr

gmx mdrun -v -deffnm md_0_1 -nt 4

#Analysis generation of RMSD files and Gyration files
echo "0" |  gmx trjconv -s md_0_1.tpr -f md_0_1.xtc -o md_0_1_noPBC.xtc -pbc mol -ur compact

echo "4" "4" | gmx rms -s md_0_1.tpr -f md_0_1_noPBC.xtc -o rmsd.xvg -tu ns

echo "4" | gmx gyrate -s md_0_1.tpr -f md_0_1_noPBC.xtc -o gyrate.xvg
#RMSF
echo "4" "4" | gmx rmsf -f md_0_1.xtc -s md_0_1.tpr -o rmsf-per-residue.xvg -ox average.pdb -res
