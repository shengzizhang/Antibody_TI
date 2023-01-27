# Antibody_TI

This pipeline incorporates FoldX, Rosetta, and thermodynamics integration (TI) to identify beneficial mutations in antibody optimization. The pipeline first uses FoldX and Rosetta to perform saturation mutagenesis and to remove deleterious mutations. TI is then used to identify beneficial mutations from a list of potential beneficial mutations predicted by FoldX and/or Rosetta.

Required software:
FoldX, Rosetta, Torque, Amber v20 or above, python2, perl with modules threads and Cwd, R with ggplot2 module. The absolute paths to these software should be added to the scripts. The pipeline is tested on Linux only.

Step 1:
A pdb file of an antibody/antigen complex or other protein complex is required. An example pdb file is in the Data/example_step1 folder. Water, ion, and other post translation modifications should be removed. For antibody, the heavy and light chains should be labeled H and L in the pdb file.

Type saturation_mutagenesis.pl without parameters to see options.

saturation_mutagenesis.pl -pdb ../Data/ example_step1/2BDN.pdb -mut “H25,H26” -o 2BDN -p foldx -c HL,A

saturation_mutagenesis.pl -pdb ../Data/ example_step1/2BDN.pdb -mut “H25,H26” -t 8 -o 2BDN -p rosetta -c HL,A
To run rosetta, a script, cart2.script in the /Data/ example_step1/ folder is required.

Step 2:
To prepare structures for TI simulation, the original structure should be reprocessed using pdb4amber. Waters should be named WAT. The processed structure should be separated to antibody only and receptor only files. If water from the antibody chain will be included (optional), they need to be removed from the antibody structure and included in a separate file. 

pdb4amber -i ../Data/ example_step1/2BDN.pdb --reduce --add-missing-atoms -o temp.pdb

Antibody_TI.pl -lig antibody.pdb -rec receptor.pdb -watlig antibody_WAT.pdb -pos 173 -length_lig1 117 -length_lig 224 -aawt THR -aamt ASP -step 1 -ssfile ./temp_sslink -length_receptor 68

After the TI simulation, change the directory to the free_energy folder.
Free_energy_calculation.pl


