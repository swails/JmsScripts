#!/usr/bin/perl -w

# find GaussRead.pm
use FindBin qw($RealBin);
use Getopt::Long;
use Data::Dumper;
use lib $RealBin;
#use lib $exe_dir;
use GaussRead;
use strict;
use warnings;

my (@file,$mol);
my($jobs,$crds);

my(%optctl);
my($file, $do_dump,$skip_check);
%optctl = (
           file => \$file,
           dump => \$do_dump,
           skip => \$skip_check, # This should remain undocumented to prevent people from 
                                 # doing stupid things
          );

&GetOptions(\%optctl, "file=s", "dump!", "skip!");

if (! $file){
    print "You must specify a filename!\n";
    print "Usage:\n";
    print "Energy.pl -file FILENAME [-dump]\n";
    print 'Where FILENAME is a Gaussian output file.',"\n";
    print '-dump for a dump to screen of job and crd hashes.',"\n";
    print '-skip will skip the consistency check. DONT USE if the file is to be added to the database';
    print "\n";
    exit 1;
}

#remove extraneous .gjfs .out or .llod
($mol = $file) =~ s/\.out|\.gjf|\.llod//g; 
# remove path name
$mol =~ s!.*/!!;
#chop off anything else after the first . that is alone
$mol =~ s/([+\-\w])\.[+\-\w].*/$1/;

print "Using \'$mol\' as york name\n";

($jobs,$crds)=processGauss($file,$mol);   

if($do_dump){
print Dumper($jobs);
print Dumper($crds);
}

#do consistency check after we dump stuff out
if(!$skip_check){
  die if not(&ConsistencyCheck($jobs,$crds,$file));
}

my ($data,$i,%calcd,@hl,@solv,$nsolv,%dG);
my ($au2kcal) = '627.5093916934901017157244496';
my ($hlvl_basis)='6-311++G(3df,2p)';
my ($hlvl_theory)='B3LYP';
# my ($hlvl_theory)='M062X';
# my ($hlvl_theory)='MP2';
# examine gas phase jobs
#print "jobs=",$jobs,"\n";
print "at jobs=",$#{$jobs},"\n";
#print "num at jobs=",$#{@{$jobs}},"\n";
for $i (0 .. $#{$jobs}){
  $data=$jobs->[$i];
#  print "$i:   ",$data->{energy}*$au2kcal,"\n";
  next if ($data->{SolvType} ne 'GAS');
  if ($data->{type} eq 'MIN' or $data->{type} eq 'TS'){
    die "More than one gas phase optimization\n" if(defined($calcd{E_ll}));
    $calcd{E_ll}=$data->{energy};
  }
  if ($data->{route}=~/freq/i){
    die "more than one GasPhase frequency Job\n" if (defined($calcd{zvpe}));
    $calcd{zvpe}=$data->{zpve};
    $calcd{ThermU}=$data->{Therm_Ucorr};
    $calcd{ThermH}=$data->{Therm_Hcorr};
    $calcd{ThermG}=$data->{Therm_Gcorr};
  }
  if (uc($data->{basis_set}) eq uc($hlvl_basis) and uc($data->{theory}) eq uc($hlvl_theory)){
    if (defined($data->{energy})){
      if (not defined($calcd{E_hl})){
	$calcd{E_hl}=$data->{energy};
	$calcd{E_hldat}=$i;
      }
      push (@hl, $i);
    }
  }
}
for $i (1 .. $#hl){
  die "High level energies not consistent ($calcd{E_hl} != $jobs->[$hl[$i]]->{energy})\n" 
	if($calcd{E_hl}!=$jobs->[$hl[$i]]->{energy});
}
if(defined($calcd{zvpe}) and defined($calcd{E_hl})){
  $calcd{E_0}=($calcd{zvpe}+$calcd{E_hl})*$au2kcal;
  $calcd{U}=($calcd{ThermU}+$calcd{E_hl})*$au2kcal;
  $calcd{H}=($calcd{ThermH}+$calcd{E_hl})*$au2kcal;
  $calcd{G}=($calcd{ThermG}+$calcd{E_hl})*$au2kcal;
  $calcd{TS}=$calcd{H}-$calcd{G};
}else{
  $calcd{E_0}='';
  $calcd{U}='';
  $calcd{H}='';
  $calcd{TS}='';
  $calcd{G}='';
}
if (! defined($calcd{E_hl})){$calcd{E_hl}=''}else{$calcd{E_hl}*=$au2kcal};
if (! defined($calcd{E_ll})){$calcd{E_ll}=''}else{$calcd{E_ll}*=$au2kcal};

# examine solvation jobs
$nsolv=0;
for $i (0 .. $#{$jobs}){
  $data=$jobs->[$i];
  next if (uc($data->{SolvType}) eq 'GAS');
  $solv[$nsolv]{job}=$i;
  if(defined($data->{dGsolv})){ #Short Cut for SM5 calcs
    $solv[$nsolv]{dGsolv}=$data->{dGsolv};
  }else{
    if(defined($data->{'<psi(0)|H|psi(0)>'})){
      $solv[$nsolv]{'<psi(0)|H|psi(0)>'}=$data->{'<psi(0)|H|psi(0)>'}*$au2kcal;
    }else{
      if(uc($jobs->[$calcd{E_hldat}]->{basis_set}) ne uc($data->{basis_set}) 
	 or uc($jobs->[$calcd{E_hldat}]->{theory}) ne uc($data->{theory})){
	print "Warning: solvation Basis/thoery inconsistant with HL Basis/theory\n";
      }
      $solv[$nsolv]{'<psi(0)|H|psi(0)>'}=$calcd{E_hl};
    } 
    # we calculate and store this information so people can print it later (if they want)
    $solv[$nsolv]{'<psi(f)|H|psi(f)>'}=$data->{'<psi(f)|H|psi(f)>'}*$au2kcal;
    $solv[$nsolv]{SolutePolarization}=($solv[$nsolv]{'<psi(f)|H|psi(f)>'}-$solv[$nsolv]{'<psi(0)|H|psi(0)>'});
    $solv[$nsolv]{'PolarSolute-Solvent'}=$data->{'PolarSolute-Solvent'};
    $solv[$nsolv]{TotalElectro}=$solv[$nsolv]{'PolarSolute-Solvent'}+$solv[$nsolv]{SolutePolarization};
    $solv[$nsolv]{RepulsionE}=$data->{RepulsionE};
    $solv[$nsolv]{DispersionE}=$data->{DispersionE};
    $solv[$nsolv]{CavitationE}=$data->{CavitationE};
    $solv[$nsolv]{TotalNonElectro}=$solv[$nsolv]{RepulsionE}+$solv[$nsolv]{DispersionE}+$solv[$nsolv]{CavitationE};
    $solv[$nsolv]{dGsolv}=$solv[$nsolv]{TotalNonElectro}+$solv[$nsolv]{TotalElectro};  # this is what gets printed
  }
  $nsolv++;
}

# create a hash of unique deltaG calculations
my $index;
for $i (0 .. $#solv){ 
  $data=$jobs->[$solv[$i]{job}];
  $index=$data->{theory}.'/'.$data->{basis_set}.':'.$data->{SolvModel}.':'.$data->{SolvRad};
  $dG{$index} = $i unless defined($dG{$index});
} 

my $gform='@<<<<<<<<<<< @<<<<<<<<<<< @<<<<<<<<<<< @<<<<<<<<<<< @<<<<<<<<<<< @<<<<<<<<<<< @<<<<<<<<<<<';
print '*' x 80,"\n";
print " Gas Phase Analysis: #energies in kcal/mol#\n";
print " Energies follow convention: E(HL) + E(thermal correction)\n"; 
print "  SCF(LL)       SCF(HL)        E_0          U         Enthalpy        TS      Free Energy\n";
$^A="";
formline($gform,$calcd{E_ll},$calcd{E_hl},$calcd{E_0},$calcd{U},$calcd{H},$calcd{TS},$calcd{G});
print "$^A\n";
print '*' x 80,"\n";
print "Solvation Analysis: #energies in kcal/mol#\n";
print "  THEORY                  Model            Radii          dG(solv)\n";
my $dgform='@<<<<<<<<<<<<<<<<<<<<<<<  @<<<<<<<<<<<<<<  @<<<<<<<<<<  @<<<<<<<<<<<<';
foreach $i (sort keys(%dG)){
  $nsolv = $solv[$dG{$i}]{job};
  $^A="";
  $data=$jobs->[$nsolv];
  formline($dgform,$data->{theory}.'/'.$data->{basis_set},$data->{SolvModel},$data->{SolvRad},$solv[$dG{$i}]{dGsolv});
  print "$^A\n";
}
print '*' x 80,"\n";
