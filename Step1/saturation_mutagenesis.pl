#!/usr/bin/perl
use strict;
use threads;
use Cwd;
my $dir = getcwd;
if(!@ARGV||@ARGV%2>0){
die'
Usage:
-pdb  pdb file absolute path
-mut list of mutations: "H25,H27,L34". Postions should be ordered by 
-t number of threads, rosetta only
-o output file prefix
-p program for saturation mutagenesis: rosetta or foldx 
-s sequence of the antibody heavy and light chain concatenated, must match PDB file order
-c chain for antibody and antigen; e.g., HL,A
';
}

my %para=@ARGV;
my $name='Ab';
if($para{'-o'}){$name=$para{'-o'};}
my $n=1;
my $pdb=$para{'-pdb'};
my $alo='rosetta';
if($para{'-p'}){$alo=$para{'-p'};}
if(!$para{'-c'}){$para{'-c'}='HL,A';}
my %seq_pdb=();&pdb_seq();

if(!$para{'-mut'}){die "no mutation specified\n";}

my @postest=split/\,/,$para{'-mut'};
my $rosetta_folder='~/work/software/Rosetta/';
my $foldX='~/work/software/foldxLinux64/foldx';
if(!$foldX||!$rosetta_folder){die "please set the absolute path to foldX and rosseta\n";}
if($alo=~/rosetta/){
	system("$rosetta_folder/main/source/bin/relax.static.linuxgccrelease  -use_input_sc -ignore_unrecognized_res -nstruct 1 -relax:cartesian -score:weights ref2015_cart -relax:min_type lbfgs_armijo_nonmonotone -relax:script cart2.script -overwrite -multithreading:total_threads 12  -fa_max_dis 9.0 -s $pdb");
$pdb=~s/.pdb$/_0001.pdb/;
}
elsif($alo=~/foldx/){
  	&pdb_repair();
  	$pdb=~s/.pdb$/_Repair.pdb/;
}

my %aa=('A',1,'D',1,'E',1,'F',1,'H',1,'I',1,'K',1,'L',1,'M',1,'N',1,'P',1,'Q',1,'R',1,'S',1,'T',1,'V',1,'W',1,'Y',1,'G',1);
my @mutants=();
open IND,">individual_list_$name\.txt";
foreach(@postest){
	  my $pos=$_; $pos=~s/\'//g;
       foreach(sort keys %aa){
		   if($seq_pdb{$pos}!~/[C$_]/){
			   print IND "$seq_pdb{$pos}$pos$_;\n";
			   push @mutants,"$seq_pdb{$pos}$pos$_";
		   }
	   }
}
close IND;
if($alo=~/foldx/){&foldx_energy();die();}

open HH,">rosetta_$name\.txt";
print HH "mutation";
	
foreach(1..$n){
			print HH "\tddG_rosetta_R$_";
}		
print HH "\tmean\n";

my $count_thread=0;
my $current=0;  
while($count_thread<@mutants){# 
	my $m=$_;
	if($current>$para{'-t'}){		
		if(threads->list()){
		   foreach(threads->list(threads::joinable)){
			my $result=$_->join();
			print  HH $result; 
			$current--;    
			} 
		}
		sleep(1);}
	else{
     threads->create({'context' => 'list'},\&rosetta,$mutants[$count_thread],$n,$pdb);
             	  print "processing thread $count_thread $mutants[$count_thread]\n";
             	  $count_thread++;
             	  $current++; 
             	  	sleep(1);       
	}
	
}   

while(threads->list()){
    foreach(threads->list(threads::joinable)){
        my $result=$_->join();
        print  HH $result;
        
    }
    sleep(1);
}


sub rosetta{	
	my ($mut,$n,$pdb)=@_;
		print "working on $mut\n";
		chdir("$dir");
		system("mkdir temp$mut");
		chdir("temp$mut");		

		my $shift=0;

		##rosetta
		if(substr($mut,1,1) eq 'H'){&mutation($mut,0);}
		elsif(substr($mut,1,1) eq 'L'){&mutation($mut,$shift);}
		chdir($dir);
		system("cd $dir/temp$mut;$rosetta_folder/main/source/bin/cartesian_ddg.static.linuxgccrelease -database $rosetta_folder/main/database  -ddg:mut_file mutation_file.txt  -ddg:iterations $n -force_iterations false  -ddg::score_cutoff 1.0 -ddg::cartesian -ddg::dump_pdbs true -bbnbrs 1 -fa_max_dis 9.0 -interface_ddg 2  -ddg::legacy false  -score:weights talaris2014_cart  -overwrite -mute all -s ../$pdb ");
		chdir($dir);
		my @mutfile=<$dir/temp$mut/MUT*bj*.pdb>;
		my @WT=<$dir/temp$mut/WT*bj*.pdb>;
		my $p=@WT;
		my $total=0;
		my $result=$mut;

	    open YY,"$dir/temp$mut/mutation_file.ddg" or die "$mut rosetta ddG file not found\n";
	    my %wt=();
	    my %mut=();
	    $total=0;
	    foreach(<YY>){
			if($_=~/COMPLEX\:   Round(\d+)\: WT_\: +([^ ]+)/){$wt{"R$1"}=$2;}
			elsif($_=~/COMPLEX\:   Round(\d+)\: MUT[^ \t]+\: +([^ ]+)/){$mut{"R$1"}=$2;}
		}
		foreach(1..$n){
			my $t=$mut{"R$_"}-$wt{"R$_"};
			$total+=$t;
			$result.="\t$t";
		}
		$result.="\t".($total/$n)."\n";
		system("rm -Rf ./temp$mut");
        return $result;
		#system("rm WT*pdb MUT*pdb *.ddg *.log");
}
sub pdb_repair{	
    open CF,">foldx_commands.txt";
    print CF "<TITLE>FOLDX_runscript;\n<JOBSTART>#;\n<PDBS>$para{'-pdb'};\n<BATCH>#;\n<COMMANDS>FOLDX_commandfile;\n<RepairPDB>#,repaired;\n<END>#;\n<OPTIONS>FOLDX_optionfile;\n<Temperature>298;\n<R>#;\n<pH>7;\n<IonStrength>0.050;\n<rmsd_pdb>true;\n<water>-CRYSTAL;\n<metal>-CRYSTAL; \n<VdWDesign>0;\n<repair_Interface>ONLY; \n<OutPDB>true; \n<numberOfRuns>1;\n<pdb_hydrogens>false; \n<complex_with_DNA> false; \n<END>#;\n<JOBEND>#; \n<ENDFILE>#;";
	close CF;
	system("$foldX --command=RepairPDB --pdb=$para{'-pdb'}");
 
	unlink "foldx_commands.txt";
}
sub foldx_energy{
	my ($mut_file)=@_;
	open ZZ,">$para{'-o'}\_freeenergy_change.txt";
	print ZZ "Mutations\tInteraction_energy_change\tGroup1_stability_change\tGroup2_stability_change\n";
    my $file=$pdb;
    $file=~s/\.pdb//;
	&foldx();
	
	unlink <*$file*fxout>,<*$file*_*.pdb>, 'rotabase.txt';
	rmdir 'molecules';
}
sub pdb_seq{
	open HH,"$para{'-pdb'}" or die "can not open $para{'-pdb'}\n";
	my %aaname_revt=reverse ('A','ALA','C','CYS','D','ASP','E','GLU','F','PHE','G','GLY','H','HIS','I','ILE','K','LYS','L','LEU','M','MET','N','ASN','P','PRO','Q','GLN','R','ARG','S','SER','T','THR','V','VAL','W','TRP','Y','TYR');
	while(<HH>){
		if($_=~/ATOM..(.{5})..C .{2}(...) (.)(....)/){
		my $pose=$4; 
		my $chain=$3;
		my $aa=$2;	
		$pose=~s/ +//g;
		 $seq_pdb{"$chain$pose"}=$aaname_revt{$aa};
		}
	}
	
}

sub mutation{
	my ($mut1,$shift)=@_;
	my @muts=split/\,/,$mut1;
	my $line='';
	
	foreach(@muts){
		~/(\d+)/;
		my $po=$1+$shift;
		$line.=substr($_,0,1).' '.$po.' '.substr($_,-1,1)."\n";
	}
	open M,">mutation_file.txt";
	my $num=@muts;
	print M "total $num\n$num\n$line";
	
}
sub foldx{
    #my ($pdbfile,$chain)=@_;	
    print "pdbfile $pdb\n";
    open ML,"individual_list_$name\.txt" or die "mutation list not found\n";
    my @mutants=<ML>;
    close ML;

    system("$foldX --command=BuildModel --pdb=$pdb --mutant-file=individual_list_$name\.txt");   
    my $i=1;
    my $file=$pdb;
    $file=~s/\.pdb//;
    foreach(@mutants){
    	chomp;
    	$_=~s/\;//;		
    	if($_!~/[A-Z0-9a-z]/){next;}
    	my $m=$_;
    my $mutant=$file."_$i";
    my $pdb_wt='WT_'.$mutant;
    print "$mutant $pdb_wt\n";
	  system("$foldX --command=AnalyseComplex --pdb=$mutant\.pdb --analyseComplexChains=$para{'-c'}");
    
    system("$foldX --command=AnalyseComplex --pdb=$pdb_wt\.pdb --analyseComplexChains=$para{'-c'}");
    my @mut_energy=&energy("Summary_$mutant\_AC.fxout");
    my @wt_energy=&energy("Summary_$pdb_wt\_AC.fxout");
	  print ZZ "$m\t",$mut_energy[0]-$wt_energy[0],"\t",$mut_energy[1]-$wt_energy[1],"\t",$mut_energy[2]-$wt_energy[2],"\n";
	  $i++;
  }
}


sub energy{
	my ($file,$id)=shift;
	  open HH,"$file" or die "foldx running wrong\n";	  
	  my $wt_energy='';
	  my $wt_stability='';
	  my $mu_energy='';
	  my $mu_stability='';
	  my $wt_stability_grp2='';
	  my $mu_stability_grp2='';
	  while(<HH>){
	  	chomp;
	  	if($_=~/$id/){
	  		 my @line=split/\t/,$_;
	  		   $wt_energy=$line[5];
	  		   $wt_stability=$line[6];
	  		   $wt_stability_grp2=$line[7];
	  		}
	  	}
	  	return $wt_energy,$wt_stability,$wt_stability_grp2;
	
}
