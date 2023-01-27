#!/usr/bin/perl
use strict;

use GD::Graph::linespoints;
use GD::Graph::Data;
if(@ARGV%2>0||@ARGV==0){die "
Usage: free_energy_calculation.pl 
no parameter is required. If any lambda window need to be excluded, use '0.11505|0.20634'
	\n";
	}
my @folders=();
opendir my $dh, "./"
  or die "$0: opendir: $!";

while (defined(my $name = readdir $dh)) {
  next unless -d "./$name";
  if($name ne '.' && $name ne '..' ){push @folders,$name;}
}

my $offset=251;
my %energy=();
my %dG=();	
my %stageddG=();
my $frame=1500;
my $nst=3;
my $ns=$frame/$nst;
my @lam=(0.00922,0.04794, 0.11505, 0.20634, 0.31608, 0.43738,  0.56262, 0.68392, 0.79366, 0.88495, 0.95206,0.99078);
my $path_getdvdl='~/work/software/TI/Step2/getdvdl.py';
my $path_dvdlR='~/work/software/TI/Step2/dadl.R';
my $path_python2='/usr/bin/python2';
if(!$path_getdvdl||!$path_dvdlR||!$path_python2){die "please set the absolute path to getdvdl and dvdl.R\n";}
foreach('complex','protein'){
	chdir($_);
	my @f=();
	my $folder=$_;
  
	foreach(<prod.en*>){
		my $file=$_;
		$file=~/^prod.en\.0+(.+)/;
		my $lam=$1-1;
		if($lam[$lam]!~/$ARGV[0]/){
		   system("cp $file prod$lam[$lam].en");
	    }	   
	    elsif($lam[$lam]=~/$ARGV[0]/){
		   system("mv prod$lam[$lam].en prod$lam[$lam].en1");
	    }
	}
	$frame++;
	for(my $i=0;($i+1)*$ns<$frame;$i++){
		my $start=$i*$ns+1;
		my $end=($i+1)*$ns+1;
		`$path_python2 $path_getdvdl $start \'prod*.en\' $end >dadl_$i.txt`;
		$stageddG{$folder}{'g'.$i}=`tail -1 dadl_$i.txt`;
		$stageddG{$folder}{'g'.$i}=~s/.+\, (.+)\)/$1/;
		chomp $stageddG{$folder}{'g'.$i};
		`$path_python2 $path_getdvdl $offset \'prod*.en\' $end >dadl_l$i.txt`;
		$stageddG{$folder}{'gl'.$i}=`tail -1 dadl_l$i.txt`;
		$stageddG{$folder}{'gl'.$i}=~s/.+\, (.+)\)/$1/;
		chomp $stageddG{$folder}{'gl'.$i};
	}

	`$path_python2 $path_getdvdl $offset \'prod*.en\' 999999 >dadl.txt`;

	open HH,"dadl.txt" or die;
	@f=<HH>;
	my @lams=();
	my @stdl=();
	my @stdh=();
	my @mean=();
	foreach(@f){
		chomp;
        if($_=~/^(\d.+)/){
			my @l=split/ /,$1;
			push @lams,$l[0];
			push @mean,$l[1];
			push @stdl,$l[1]-$l[3];
			push @stdh,$l[1]+$l[3];
			
		}
	   if($_=~/dG =\'\, (.+)\)/){
		   $energy{$folder}=$1;
		   
		   }
   }
   
    open combined, ">combined_energy.txt";
    open combineddvdl, ">combined_dvdl.txt";   
    my $steps=20;
	foreach(<prod*.en>){
		open ZZ,"$_";
		my $outf=$_;
		my @data1=();
		my @data3=();
		my @data2=();
		my @data=();
		my $step1=1;
		$outf=~s/en/energy.txt/;		
		my $i=10;
		while($i>0){
			<ZZ>;
			$i--;
		}
		while(<ZZ>){
			chomp;
			if($_=~/^L0/){
			    my @ll=split/[ \t]+/,$_;
			    
			    push @data1,$ll[1];
			    push @data3, $ll[3];
			    print combined "$folder\t$outf\t$steps\t$ll[3]\n";
			    
			    $steps+=20;
		    }
			elsif($_=~/^L9/){
				my @ll=split/[ \t]+/,$_;			    	
			    push @data2,$ll[$#ll];			    
			    print combineddvdl "$folder\t$outf\t$step1\t$ll[$#ll]\n";
			    $step1++;
			}
		}
		@data=([@data1],[@data2]);
	}

	system("dadl.R");
	
	chdir("..");
}

print "infered ddG is\t$energy{'complex'}\t$energy{'protein'}\t",$energy{'complex'}-$energy{'protein'}, "\n";











