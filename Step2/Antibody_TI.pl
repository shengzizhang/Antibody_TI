#!/usr/bin/perl
use strict;
use Cwd qw(abs_path);
#use threads;
# This script is used to calculate free energy change

if(@ARGV%2>0||@ARGV==0){die "
Usage: TI_Amber_softcore.pl 
	-lig ligand alone structure positions must start from 1. 
	-rec receptor alone structure, must from complex without coordinate change
	-pos position in lig for mutation
	-length_lig number of residues in ligand
	-aawt widetype residues in three letter
	-aamt mutated residues
	-step 1 for equilibration 2 for production
	-length_lig1 length of the first chain if multiple ligand chains and the mutation is in the first chain
	-cuda gpu or cpu for density calculation. 1 for gpu, default. 2 for CPU
	-folder protein or complex for calculation
	-nsteps nanoseconds for production run. default: 3ns
	-ssfile disulfid bond file from pdb4amber
	-boxsize size of water box, default 6A
	-windows default 12, alternative 16
	-watlig waters for antibody
	-length_receptor aa lenght of receptor
	-glycanfile glycan link to Asn. the order will be postion of ASN+tab+pos glycan. 
	\n";
	}
	
my %para=@ARGV;
if(!$para{'-lig'} || !$para{'-rec'} || !$para{'-pos'} || !$para{'-length_lig'}|| !$para{'-length_lig1'}|| !$para{'-length_receptor'}){die "parameters are missing\n";}
if(!$para{'-dev'}){$para{'-dev'}='0,1';}
my $lig=$para{'-lig'};
my $complex='complex_H.pdb';
my $mut='mutation.pdb';
my $receptor=$para{'-rec'};
my $pos1=$para{'-pos'};
my $lengthlig=$para{'-length_lig'};
my $pos2=$pos1+$lengthlig;
my $pathc=abs_path();
my $mut=$para{"-aamt"}.$pos1;
my $common='@N,C,O,CA,H,HA';

my $SSprotein='';
my $SScomplex='';
my $restrain='';
my $exchange='';
if(!$para{'-step'}){$para{'-step'}=1;}
if(!$para{'-MDt'}){$para{'-MDt'}=1;}
if(!$para{'-cuda'}){$para{'-cuda'}=1;}
if(!$para{'-SSC'}){$para{'-SSC'}=2;}
if(!$para{'-boxsize'}){$para{'-boxsize'}=6;}
if(!$para{'-nsteps'}){$para{'-nsteps'}=3000000;}
if(!$para{'-nsteps_equil6'}){$para{'-nsteps_equil6'}=250000;}
my $perexchange=2000;
$exchange=$para{'-nsteps'}/$perexchange;

my $SScomplex='';
 if(!$para{'-MDenv'}){$para{'-MDenv'}='NVT';}
 
my $command="@ARGV";
my 	@lam=(0.00922,0.04794, 0.11505, 0.20634, 0.31608, 0.43738,  0.56262, 0.68392, 0.79366, 0.88495, 0.95206,0.99078);
my $amber="source ~/work/software/Amber/amber20_src/amber.sh\n";
my $antibody_TI="~/work/software/TI/Step2/Antibody_TI.pl";
if(!$amber ||!$antibody_TI){die "need to specify the absolute paths to amber source file or antibody_TI\n";}
if($para{'-aawt'} =~/PRO/ || $para{'-aamt'} =~/PRO/ ){$common='@N,C,O,CA,HA';}
elsif($para{'-aawt'} =~/GLY/ || $para{'-aamt'} =~/GLY/ ){$common='@N,C,O,CA,H';}

if(($para{'-aawt'} =~/PRO/ || $para{'-aamt'} =~/PRO/) && ($para{'-aawt'} =~/GLY/ || $para{'-aamt'} =~/GLY/)){$common='@N,C,O,CA';}

if($para{'-ssfile'} && $para{'-step'}==1){
	&generate_SS();
	}
my $MDe='
  ntb     = 1,
  ';
my $MDt="
  ntt     = 1, 
  tautp   = 2.0,
  ";

my $SSC='
    gti_lam_sch = 1,
    scalpha = 0.2, 
    scbeta = 50.0,
   ';	

if($para{'-step'}==1){
		&mutation($lig,$pos1,$para{'-aawt'},$para{'-aamt'});
		&tleapcombine(1,1);
		my @nion=&ions();
		&tleapcombine(@nion);
		
		&Timerge();
        system("ambpdb -p merged_protein.prmtop -c merged_protein.inpcrd -conect >merged_protein.pdb");
        system("ambpdb -p merged_complex.prmtop -c merged_complex.inpcrd -conect >merged_complex.pdb");
		$pos2=$lengthlig+1;
		if($pos1 <= $para{'-length_lig1'}){$pos2=$para{'-length_lig1'}+1;}
		&relaxation();
		
		$command=~s/\-step[ \t]+1/\-step 2 /;
		$command=~s/\-cuda[ \t]+\d/\-cuda 1 /;
		$command="$antibody_TI $command -step 2";	
		chdir($pathc);	
		&checkrelaxation();

	}
else{
		$pos2=$lengthlig+1;
		if($pos1 <= $para{'-length_lig1'}){$pos2=$para{'-length_lig1'}+1;}

			&free_energy($para{'-folder'});

}
##########
sub checkrelaxation{
	my $p=0;
	my $co=0;
    my $total=0;
	while($total<2){
		if(-e "$pathc/relaxation/protein/equil$lam[0].out" && $p==0){
			my $pp=`tail -1 $pathc/relaxation/protein/equil$lam[0].out`;
			if($pp=~/Total wall time:/){
				$p=1;
					for(my $i=0;$i<@lam;$i++){
						my $j=$i+1;
						if($j<10){$j="00$j";}
						else{$j="0$j";}
						system("mv $pathc/relaxation/protein/equil$lam[$i].rst7 $pathc/relaxation/protein/equil$j.rst");
					}
				chdir($pathc); 
				&pbs("$command -folder protein",'FElig');}
		}
	   if(-e "$pathc/relaxation/complex/equil$lam[0].out" && $co==0){
			my $pp=`tail -1 $pathc/relaxation/complex/equil$lam[0].out`;
			if($pp=~/Total wall time:/){
				$co=1;
					for(my $i=0;$i<@lam;$i++){
						my $j=$i+1;
						if($j<10){$j="00$j";}
						else{$j="0$j";}
						system("mv $pathc/relaxation/complex/equil$lam[$i].rst7 $pathc/relaxation/complex/equil$j.rst");
					}				
				chdir($pathc);&pbs("$command -folder complex",'FEcom');
				}
		}
		$total=$p+$co;
		
		sleep(30);
   }
}


sub relaxation{
	
	mkdir("relaxation");
	chdir("relaxation");
	system("rm -r *");
    mkdir("protein");
	chdir("protein");
	&min1();
    &min2();
	&heat();
	&density();
    #&density1();
	&equil_rel(1,5);
	&equil_rel(2,4);	
	&equil_rel(3,3);	
	&equil_rel(4,2);	
	&equil_rel(5,1);
	&equil_rel6(6);	
	my $eq=&equilibration(\@lam,'protein');	

	&pbs_rel('../../merged_protein',$eq);
	system("qsub pbs.sh");
	
    chdir("..");
    mkdir("complex");
	chdir("complex");
	&min1();
    &min2();
	&heat();
	&density();
    #&density1();
	&equil_rel(1,5);
	&equil_rel(2,4);	
	&equil_rel(3,3);	
	&equil_rel(4,2);	
	&equil_rel(5,1);
	&equil_rel6(6);		
	$eq=&equilibration(\@lam,'complex');
	&pbs_rel('../../merged_complex',$eq);
	system("qsub pbs.sh");
}
sub free_energy{
	my ($f)=@_;
	if($f!~/\w/){die "folder not specificied\n";}
    mkdir('free_energy');
    chdir('free_energy');
	my $m=6;
    mkdir($f);
	chdir($f);
	system("rm *");

    #prepair input file
		open IN,">input.dat";
		print IN "HAMILTONIAN\nHamiltonian Replica Exchange with different mdins\n";
		foreach(@lam){print IN "$_\n";}
		close IN;
	#prepair group file
		open GR,">groupfile.ref";
		print GR "\# Replica REPNUM\n-O -i mdin.rep.REPNUM -p $pathc/merged_$f\.prmtop -c $pathc/relaxation/$f/equilREPNUM.rst -o mdout.rep.REPNUM -r rst7.rep.REPNUM -e prod.en.REPNUM";
		close GR;
	# prepair mdin
		open MD,">mdin.ref";
		print MD &HREMD();
		close MD;
	# generate input files for MD
	system("genremdinputs.py -inputs input.dat -groupfile groupfile.ref -i mdin.ref -O");
	my $thread=@lam;
	system(" mpirun -np $thread pmemd.cuda.MPI -ng $thread -rem 3 -groupfile groupfile -remlog HREMD.log");
	#&qsubstr($task,$mut.'_'.$f);
	chdir("../..");
   
}

############
sub HREMD{
	my $str="
NVT MD 
 \&cntrl
  ntx     = 5,
  irest   = 1,
  ntpr    = 50000,
  ntwx    = 0,
  ntwr    = 500000,
  ntwe	  = 2000,
  numexchg= ".$exchange.",
  nstlim  = ".$perexchange.",
    
  ig	  = RANDOMNUM,	
  ntf     = 1,
  cut     = 8.0,
 
  nscm    = 0,
  t       = 0.0,
  dt      = 0.001,
  iwrap   = 1,
  
  temp0   = 298.0, !temperature control
  tempi  =298,
  gremd_acyc =1,
  ".$MDe."
  ntt     = 3, 
  tautp   = 2.0,
  gamma_ln = 2.0, 
  
  ntc     = 1,
  tol     = 0.000001,
 
   icfe = 1, 
   ifsc=1, 
   clambda = HAMILTONIAN,
   logdvdl = 0,".$SSC."
   timask1 = \'\:$pos1\', 
    
   timask2 = \'\:$pos2\', 
   
   scmask1 = \':$pos1 & !$common\', 
   scmask2 = \':$pos2 & !$common\',
/

";
return $str;
}

###equil
sub equilibration{
	my ($lam,$f)=@_;
    my $m=int(@$lam/2);
    my @lam=@$lam;
    my $command='';
	$command=&equil_singlepath($lam[$m],$pos1,$pos2,"$pathc/merged_$f\.prmtop","equi6.rst7");
	for(my $i=$m+1;$i<@lam;$i++){		
		$command.=&equil_singlepath($lam[$i],$pos1,$pos2,"$pathc/merged_$f\.prmtop","equil".$lam[$i-1].'.rst7');	
	}
	for(my $i=$m-1;$i>=0;$i--){		
		$command.=&equil_singlepath($lam[$i],$pos1,$pos2,"$pathc/merged_$f\.prmtop","equil".$lam[$i+1].'.rst7');
	}	
	

	return $command;
}
##########
sub qsubstr{
	my ($file,$title)=@_;

open HH,">pbs.sh";	
print HH "
#!/bin/bash
#PBS -N EFprod$title
#PBS -l walltime=999::00
#PBS -l nodes=1:ppn=24:gpus=4
#PBS -q batch
#PBS -j oe

cd \"\$PBS_O_WORKDIR\"
source ~/.bash_profile
export PATH=\$PATH:/usr/local/bin
$amber

$file

";
system("qsub pbs.sh");
}

sub pbs{
	my ($file,$stage)=@_;
my $th=@lam;
open HH,">pbs$stage.sh";	
print HH "
#!/bin/bash
#PBS -N FEM$mut\_$stage
#PBS -l walltime=999::00
#PBS -l nodes=1:ppn=$th\:gpus=4
#PBS -q batch
#PBS -j oe
#PBS -p 1023

cd \"\$PBS_O_WORKDIR\"
source ~/.bash_profile
export PATH=\$PATH:/usr/local/bin

$amber

$file
";
system("qsub pbs$stage.sh");
}

sub tleapcombine{
	 my @ions=@_;
open HH,">tleapcombine.in";
print HH '
source leaprc.GLYCAM_06j-1
source leaprc.protein.ff14SB
source leaprc.DNA.OL15
source leaprc.lipid17
source leaprc.water.tip3p
source leaprc.gaff2

wtlig = loadpdb ',$lig,' #Load PDB file for protein-ligand complex
mtlig = loadpdb mutation.pdb
';
if($para{'-watlig'}) {
	print HH "watlig= loadpdb $para{'-watlig'}
receptor= loadpdb $receptor

protein = combine \{wtlig mtlig watlig\}
complex = combine \{wtlig mtlig receptor watlig\}
";	
	}
else{
	print HH '
receptor= loadpdb ',$receptor,'

protein = combine {wtlig mtlig }
complex = combine {wtlig mtlig receptor}';
}

print HH '
',$SSprotein,$SScomplex,'
set default nocenter on
';
if($ions[0]==1){
print HH '
# create protein in solution
solvateBox protein TIP3PBOX ',$para{'-boxsize'},' 0.75
#addions2 protein Cl- 0 #Add Cl- ions to neutralize the system
#addions2 protein Na+ 0
#addIons2 protein Cl- ',$ions[0],' Na+ ',$ions[0],'

savepdb protein protein.pdb
saveamberparm protein protein.prmtop protein.inpcrd

# create complex in solution
solvateBox complex TIP3PBOX ',$para{'-boxsize'},' 0.75
#addions2 complex Cl- 0 #Add Cl- ions to neutralize the system
#addions2 complex Na+ 0
#addIons2 complex Cl- ',$ions[1],' Na+ ',$ions[1],'
savepdb complex complex.pdb
saveamberparm complex complex.prmtop complex.inpcrd
quit #Quit tleap program
';		
	}
	else{
		print HH '
# create protein in solution
solvateBox protein TIP3PBOX ',$para{'-boxsize'},' 0.75
addions protein Cl- 0 #Add Cl- ions to neutralize the system
addions protein Na+ 0
addIons2 protein Cl- ',$ions[0],' Na+ ',$ions[0],'
savepdb protein protein.pdb
saveamberparm protein protein.prmtop protein.inpcrd

# create complex in solution
solvateBox complex TIP3PBOX ',$para{'-boxsize'},' 0.75
addions complex Cl- 0 #Add Cl- ions to neutralize the system
addions complex Na+ 0
addIons2 complex Cl- ',$ions[1],' Na+ ',$ions[1],'
savepdb complex complex.pdb
saveamberparm complex complex.prmtop complex.inpcrd
quit #Quit tleap program
';
}
close HH;
system("tleap -f tleapcombine.in >tleapout.txt ");

}
sub ions{
	open HH,"tleapout.txt" or die 'tleap file not exist\n';
	my @box=();
	my $ionl=0;
	my $ionc=150;#150nM
	while(<HH>){
		if($_=~/  Volume: ([\d\.]+)/){
			push @box,$1;
		}
	}
	foreach(0,1){
		if(!$box[$_]){die "$_ box size is zero\n";}
		$box[$_]=sprintf("%.0f",$ionc*6.022*$box[$_]*1e-7);
	}
	print "ions added @box\n";
	return @box;
}
sub pbs_rel{
	my ($file,$eq)=@_;
	my $ncore=1;
	my $ngpu=':gpus=1';
	my $dens="pmemd.cuda -i density.in -p $file\.prmtop -c heat.rst7 -ref $file\.inpcrd -O -o density.out  -r density.rst7 -x density.mdcrd";
if($para{'-cuda'}>1){
	$ncore=32;
	$dens="mpirun -np 32 pmemd.MPI -i density.in -p $file\.prmtop -c heat.rst7 -ref $file\.inpcrd -O -o density.out  -r density.rst7 -x density.mdcrd";
	$ngpu='';
	}
open HH,">pbs.sh";	
print HH "
#!/bin/bash
#PBS -N $mut\_relax
#PBS -l walltime=999::00
#PBS -l nodes=1:ppn=$ncore$ngpu
#PBS -q batch
#PBS -j oe


cd \"\$PBS_O_WORKDIR\"
source ~/.bash_profile
export PATH=\$PATH:/usr/local/bin
$amber

pmemd.cuda -i min1.in -p $file\.prmtop -c $file\.inpcrd -ref $file\.inpcrd -O -o min1.out -r min1.rst7 
pmemd.cuda -i min2.in -p $file\.prmtop -c min1.rst7 -ref $file\.inpcrd -O -o min2.out  -r min2.rst7 
pmemd.cuda -i heat.in -p $file\.prmtop -c min2.rst7 -ref $file\.inpcrd -O -o heat.out  -r heat.rst7 
$dens
pmemd.cuda -i equi1.in -p $file\.prmtop -c density.rst7 -ref density.rst7 -O -o equi1.out  -r equi1.rst7 
pmemd.cuda -i equi2.in -p $file\.prmtop -c equi1.rst7 -ref density.rst7 -O -o equi2.out  -r equi2.rst7 
pmemd.cuda -i equi3.in -p $file\.prmtop -c equi2.rst7 -ref density.rst7 -O -o equi3.out  -r equi3.rst7
pmemd.cuda -i equi4.in -p $file\.prmtop -c equi3.rst7 -ref density.rst7 -O -o equi4.out  -r equi4.rst7 
pmemd.cuda -i equi5.in -p $file\.prmtop -c equi4.rst7 -ref density.rst7 -O -o equi5.out  -r equi5.rst7 
pmemd.cuda -i equi6.in -p $file\.prmtop -c equi5.rst7 -O -o equi6.out  -r equi6.rst7 -x equi6.mdcrd
$eq

";
}


sub Timerge{
	my $p1="1\-$lengthlig";
	my $p2=($lengthlig+1).'-'.$lengthlig*2;
	system("parmed protein.prmtop <<_EOF
	loadRestrt protein.inpcrd
setOverwrite True
tiMerge :$p1 :$p2 :$pos1 :$pos2
 
outparm merged_protein.prmtop merged_protein.inpcrd
quit	
	");
	
system("parmed complex.prmtop <<_EOF
	loadRestrt complex.inpcrd
setOverwrite True
tiMerge :$p1 :$p2 :$pos1 :$pos2

outparm merged_complex.prmtop merged_complex.inpcrd
quit	
	");
		
}
sub mutation{
	my ($file,$posi,$aaw,$aam)=@_;
	open HH,"$file" or die;
	open YY, ">mutation.pdb";
	my $checkp=0;
    while(<HH>){
	  if($_=~/ATOM|HETATM/){
	   my $aanum=substr($_,length('ATOM   3326  OXT GLU L'),4);
	   $aanum=~s/ +//g;
	   my $aachain=substr($_,length('ATOM      7  O   SER '),1);
	   if($aaw=~/HIE|HID|HIP|HIS/){$aaw='HIS|HIE|HID|HIP';}
	   if($aaw=~/GLH|GLU/){$aaw='GLU|GLH';}
	   if($aaw=~/ASH|ASP/){$aaw='ASP|ASH';}
	   if($aaw=~/CYX|CYS/){$aaw='CYS|CYX';}
       if($aanum==$posi && $_=~/$aaw/ ){
		   ~s/$aaw/$aam/;
		   $checkp=1;
		  my $atom=substr($_,length('ATOM    901  '),3);
		  $atom=~s/ +//g;
		  
	      if($atom=~/^(N|C|O|CA|H|HA)$/){
			  if($para{'-aamt'}=~/PRO/ && ($atom eq 'H' )){}
			  elsif($para{'-aamt'}=~/GLY/ && ($atom eq 'HA' )){}
			  else{print YY $_;}
			  
			  }
	   }
       else{print YY $_;}

	}
	else{print YY $_;}

}
if($checkp==0){die "$file,$posi,$aaw,$aam not found\n";}
}


sub min1{
open HH,">min1.in";	
print HH "minimisation1
 \&cntrl
  imin   = 1,

  ntpr   = 100,
  ntwe = 100,
  ntb    = 1,
  cut    = 8.0,
  nsnb   = 10,

  ntr    = 1,
  restraint_wt = 25.00,
  restraintmask=\'\!:WAT &!:Na+ & !:Cl-\',
  maxcyc = 10000,
  ntmin  = 2,
  ncyc   = 50,
  dx0    = 0.1,
  icfe = 1, ifsc = 1, clambda = 0.5, ",$SSC,"
   logdvdl = 0,
   timask1 = \':$pos1\', timask2 = \':$pos2\',
   scmask1 = \':$pos1 & !$common\', 
   scmask2 = \':$pos2 & !$common\',
 /
 
";	
	
}
sub min2{
open HH,">min2.in";	
print HH "minimisation2
 \&cntrl
  imin   = 1,

  ntpr   = 100,
  ntwe = 100,
  ntb    = 1,
  cut    = 8.0,
  nsnb   = 10,

  ntr    = 1,
  restraint_wt = 5.00,
  restraintmask=\'\!:WAT &!:Na+ & !:Cl-\',
  maxcyc = 10000,
  ntmin  = 2,
  ncyc   = 50,
  dx0    = 0.1,
  icfe = 1, ifsc = 1, clambda = 0.5, ",$SSC,"
   logdvdl = 0,
   timask1 = \':$pos1\', timask2 = \':$pos2\',
   scmask1 = \':$pos1 & !$common\', 
   scmask2 = \':$pos2 & !$common\',
 /
 
";	
	
}

sub heat{
open HH,">heat.in";	
print HH "NVT MD w/ position restraints (5 kcal/molA) and PME (sander)
 \&cntrl
  nmropt = 1,

  ntpr   = 1000,
  ntwr   = 10000,
  
  ntf    = 1,
  ntb    = 1,
  cut    = 8.0,
  nsnb   = 10,

  ntr    = 1,
   restraint_wt = 5.00,
   restraintmask=\'\!:WAT &!:Na+ & !:Cl-\',
  nstlim = 100000,
  nscm   = 1000,
  dt     = 0.002,
   icfe = 1, ifsc = 1, clambda = 0.5, ",$SSC,"
   logdvdl = 0,
   timask1 = \':$pos1\', timask2 = \':$pos2\',
   scmask1 = \':$pos1 & !$common\', 
   scmask2 = \':$pos2 & !$common\',
",$MDt,"   

  tempi  = 100.0,
  temp0  = 298,

  ntc    = 2,
  tol    = 0.000001,
 /

 &ewald
 / 

 &wt
   type='TEMP0',
   istep1 = 0, istep2 = 50000,                                      
   value1 = 0.0, value2 = 298.0
 /

 &wt type = 'END'
 /


";	
	
}

sub density{
open HH,">density.in";	
print HH "NPT MD w/ position restraints (5 kcal/molA) and PME (sander)
 \&cntrl
  ntx    = 5,
  irest  = 1,
  ntpr   = 10000,
  ntwx   = 10000,
  ntwr  = 10000,
  ntf    = 1,
  ntb    = 2,
  cut    = 8.0,
  nsnb   = 10,

  ntr    = 1,
  restraint_wt = 5.00,
  restraintmask=\'\!:WAT &!:Na+ & !:Cl-\',
  nstlim = 200000,
  nscm   = 0,
  dt     = 0.002,

 ",$MDt," 
  tempi  = 298.0,   
  temp0  = 298,
  ntp    = 1,
  pres0  = 1.0,
  taup   = 5,
  rgbmax = 15,
  ntc    = 2,
  tol    = 0.000001,
        icfe = 1, ifsc = 1, clambda = 0.5, ",$SSC,"
   logdvdl = 0,
   timask1 = \':$pos1\', timask2 = \':$pos2\',
   scmask1 = \':$pos1 & !$common\', 
   scmask2 = \':$pos2 & !$common\',
/
 
";	
	
}

sub equil_rel{
	my ($n,$restr)=@_;
open HH,">equi$n.in";	
print HH "NVT MD w/ position restraints (5 kcal/molA) and PME (sander)
 &cntrl
  ntx    = 5,
  irest  = 1, 
  ntpr   = 10000,
  ntwx   = 20000,
  ntwr   = 0,

  ntf    = 1,
  ntb    = 2,
  cut    = 8.0,
  nsnb   = 10,

  ntr    = 1, 
  restraint_wt = $restr,
  restraintmask=\'\!:WAT &!:Na+ & !:Cl-\',
  nstlim = 250000,
  nscm   = 0,

  dt     = 0.002,
  ",$MDt,"
",$MDe,"
  temp0   = 298.0, !temperature control
  tempi  = 298.0,   

  ntc    = 2,
  tol    = 0.000001,
   icfe = 1, ifsc = 1, clambda = 0.5, ",$SSC,"
   logdvdl = 0,
   timask1 = \':$pos1\', timask2 = \':$pos2\',
   scmask1 = \':$pos1 & !$common\', 
   scmask2 = \':$pos2 & !$common\',
   
   /
 
";	
	
}

sub equil_rel6{
	my ($n)=@_;
open HH,">equi$n.in";	
print HH "NVT MD w/No position restraints and PME (sander)
 \&cntrl
  ntx    = 5,
  irest  = 1,
  ntpr   = 20000,
  ntwx   = 25000,
  ntwr   = 0,
  iwrap  = 1,
  ntf    = 1,
  ntb    = 2,
  cut    = 8.0,
  nsnb   = 10,

  nstlim = $para{'-nsteps_equil6'},
  nscm   = 0,
  dt     = 0.002,

  temp0   = 298.0, !temperature control
  tempi  = 298.0,   
  ",$MDt,"
",$MDe,"
  ntc    = 2,
  tol    = 0.000001,
   icfe = 1, ifsc = 1, clambda = 0.5, ",$SSC,"
   logdvdl = 0,
   timask1 = \':$pos1\', timask2 = \':$pos2\',
   scmask1 = \':$pos1 & !$common\', 
   scmask2 = \':$pos2 & !$common\',
 /
";	
	
}

sub equil_singlepath{
	my ($clam,$pos1,$pos2,$parm,$rst)=@_;
	open HH,">equil$clam.in";
print HH "
NVT MD w/No position restraints and PME (sander)
 \&cntrl
  ntx     = 5,
  irest   = 0,
  ntpr    = 100000,
  ntwx    = 0,
  ntwr    = 100000,
  ntwe	  = 100000,
  
  ntf     = 1,
  cut     = 8.0,
  $restrain
  nstlim  = 1000000, 
  nscm    = 0,
  t       = 0.0,
  dt      = 0.001,
  iwrap   = 1,
  
  temp0   = 298.0, !temperature control
  tempi  =298,

  ",$MDt,"
",$MDe,"
  ntc     = 1,
  tol     = 0.000001,
 
   icfe = 1, 
   ifsc=1, 
   clambda = $clam,
",$SSC,"
   logdvdl = 0,
   timask1 = \'\:$pos1\', 
    
   timask2 = \'\:$pos2\', 
   
   scmask1 = \':$pos1 & !$common\', 
   scmask2 = \':$pos2 & !$common\',
/

";	
	return "pmemd.cuda -i equil$clam.in -p $parm -c $rst -ref $rst -O -o equil$clam.out  -r equil$clam.rst7 \n";
}


sub generate_SS{
	open HH,"$para{'-ssfile'}" or die "SS file not found\n";
	while(<HH>){
		chomp;
		my @l=split/[\t ]+/,$_;
		if($l[1]<=$para{'-length_lig'}){
			 $SSprotein.="bond protein.$l[0].SG protein.$l[1].SG\n";
			 $SSprotein.="bond protein.".($l[0]+$para{'-length_lig'}).".SG protein.".($l[1]+$para{'-length_lig'}).".SG\n";
			 $SScomplex.="bond complex.$l[0].SG complex.$l[1].SG\n";
			 $SScomplex.="bond complex.".($l[0]+$para{'-length_lig'}).".SG complex.".($l[1]+$para{'-length_lig'}).".SG\n";
		}
		else{
			 
			 $SScomplex.="bond complex.".($l[0]+$para{'-length_lig'}).".SG complex.".($l[1]+$para{'-length_lig'}).".SG\n";
		}
	}
		if($para{'-glycanfile'}){
	  open HH,"$para{'-glycanfile'}" or die "glycan file not found\n";
	  while(<HH>){
		chomp;
		my @l=split/[\t ]+/,$_;
		if($l[1]<=$para{'-length_lig'}){
			 $SSprotein.="bond protein.$l[0].ND2 protein.$l[1].C1\n";
			 $SSprotein.="bond protein.".($l[0]+$para{'-length_lig'}).".ND2 protein.".($l[1]+$para{'-length_lig'}).".C1\n";
			 $SScomplex.="bond complex.$l[0].SG complex.$l[1].SG\n";
			 $SScomplex.="bond complex.".($l[0]+$para{'-length_lig'}).".ND2 complex.".($l[1]+$para{'-length_lig'}).".C1\n";
		}
		else{
			 
			 $SScomplex.="bond complex.".($l[0]+$para{'-length_lig'}).".ND2 complex.".($l[1]+$para{'-length_lig'}).".C1\n";
		}
		}
	}
}


