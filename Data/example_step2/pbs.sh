
#!/bin/bash
#PBS -N TILT56D
#PBS -l walltime=999::00
#PBS -l nodes=1:ppn=1
#PBS -q batch
#PBS -j oe

cd "$PBS_O_WORKDIR"
source ~/.bash_profile
export PATH=$PATH:/usr/local/bin
source ~/work/software/Amber/amber20_src/amber.sh


../../Step2/Antibody_TI.pl -lig antibody.pdb -rec receptor.pdb -watlig antibody_WAT.pdb -amber 20 -boxsize 6 -pos 173 -length_lig1 117 -length_lig 224 -aawt THR -aamt ASP -step 1 -ssfile ./temp_sslink -length_receptor 68	
	
