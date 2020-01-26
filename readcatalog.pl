#!/usr/bin/perl

use strict;
use warnings;
use Math::Trig;
use Math::SigFigs;
use Data::Dumper;
$Data::Dumper::Sortkeys = 1;
use Array::Utils qw(:all); # TODO: figure out how to not need require this dependency.

# readcatalog.pl - script that reads catalog.txt and outputs the binary star catalog.
#
# This Perl script was custom-made for the catalog.txt. There are several idiosyncrasies that
# are specific to the way the handles systems. Notably, the order of the line in each catalog
# matters. This determines the order of the systems as outputted in the .stc file, and also
# helps the script determine what counts as a component of a previous system and what counts
# as a new system entirely.
#
# This script frequently references component designations and suffixes, (e.g. the "A" in "Sirius
# A"). I use "component designations" and "suffixes" interchangeably to refer to them.
#
# There is one optional command-line argument, which is --verbose or -v; when using this argument
# the script will print a status update, and count the number of systems and stars.

# ===================== CONVERSION FACTORS AND OTHER IMPORTANT CONSTANTS ===================== #

my $LY_TO_PARSEC = 3.26167; # internal value in Celestia
my $OBLIQUITY = 23.4392911;
my $MJ_TO_MSOL = 0.0009543;
my $SOLAR_RADIUS = 695700;

# ============================ STELLAR PARAMETERS AND NOMENCLATURE =========================== #

# This hash stores the order of each component designation as it should be in the .stc file.
# Some have the same value, because they would not be seen in the same system (e.g. A-BC, AB-C).
my %csort = (
'' => '00', 'AB-C' => '00', 'A-BC' => '00', 'AB-CD' => '00',
'AB-CE' => '00', 'AB-DE' => '00', 'AB' => '01', 'AC' => '01',
'A' => '02', 'Aa' => '03', 'Aa1' => '04', 'Aa2' => '05',
'Ab' => '06', 'Ab1' => '07', 'Ab2' => '08', 'BC' => '09',
'B' => '10', 'Ba' => '11', 'Ba1' => '12', 'Ba2' => '13',
'Bb' => '14', 'Bb1' => '15', 'Bb2' => '16', 'CD' => '17',
'CE' => '17', 'DE' => '17', 'C' => '18', 'Ca' => '19',
'Cb' => '20', 'D' => '21', 'Da' => '22', 'Db' => '23',
'E' => '24', 'Ea' => '25', 'Eb' => '26',
);

# Table matching spectral types to effective temperatures (for main-sequence stars)
# http://www.pas.rochester.edu/~emamajek/EEM_dwarf_UBVIJHK_colors_Teff.txt
my %SpTeff = (
'O5V' => 41500,
'O6V' => 39000,
'O7V' => 36500,
'O8V' => 34500,
'O9V' => 32500,
'B0V' => 31500,
'B1V' => 26000,
'B2V' => 20600,
'B3V' => 17000,
'B4V' => 16700,
'B5V' => 15700,
'B6V' => 14500,
'B7V' => 14000,
'B8V' => 12500,
'B9V' => 10700,
'A0V' => 9700,
'A1V' => 9200,
'A2V' => 8840,
'A3V' => 8550,
'A4V' => 8270,
'A5V' => 8080,
'A6V' => 8000,
'A7V' => 7800,
'A8V' => 7500,
'A9V' => 7440,
'F0V' => 7220,
'F1V' => 7030,
'F2V' => 6810,
'F3V' => 6720,
'F4V' => 6640,
'F5V' => 6510,
'F6V' => 6340,
'F7V' => 6240,
'F8V' => 6170,
'F9V' => 6060,
'G0V' => 5920,
'G1V' => 5880,
'G2V' => 5770,
'G3V' => 5720,
'G4V' => 5680,
'G5V' => 5660,
'G6V' => 5590,
'G7V' => 5530,
'G8V' => 5490,
'G9V' => 5340,
'K0V' => 5280,
'K1V' => 5170,
'K2V' => 5040,
'K3V' => 4830,
'K4V' => 4600,
'K5V' => 4410,
'K6V' => 4230,
'K7V' => 4070,
'K8V' => 4000,
'K9V' => 3940,
'M0V' => 3870,
'M1V' => 3700,
'M2V' => 3550,
'M3V' => 3410,
'M4V' => 3200,
'M5V' => 3030,
'M6V' => 2850,
'M7V' => 2650,
'M8V' => 2500,
'M9V' => 2400,
);

# Table matching spectral types to masses. For main-sequence stars O7V and later, the data
# come from Eric Mamajek's stellar parameter table. For other stars, the data come from
# Straizys and Kuriliene (1981), Ap&SS 80, p.353
# "Fundamental stellar parameters derived from the evolutionary tracks"
# https://ui.adsabs.harvard.edu/abs/1981Ap%26SS..80..353S/abstract
my %SpMass = (
'O5V' => 10**1.81, 'O5IV' => 10**1.85, 'O5III' => 10**1.89, 'O5II' => 10**1.90, 'O5Ib' => 10**1.92, 'O5Iab' => 10**1.99,
'O6V' => 10**1.70, 'O6IV' => 10**1.76, 'O6III' => 10**1.80, 'O6II' => 10**1.80, 'O6Ib' => 10**1.87, 'O6Iab' => 10**1.91, 'O6Ia' => 10**2.00,
'O7V' => 28, 'O7IV' => 10**1.65, 'O7III' => 10**1.68, 'O7II' => 10**1.71, 'O7Ib' => 10**1.76, 'O7Iab' => 10**1.83, 'O7Ia' => 10**1.92,
'O8V' => 22.9, 'O8IV' => 10**1.54, 'O8III' => 10**1.60, 'O8II' => 10**1.65, 'O8Ib' => 10**1.72, 'O8Iab' => 10**1.76, 'O8Ia' => 10**1.90,
'O9V' => 19.7, 'O9IV' => 10**1.45, 'O9III' => 10**1.49, 'O9II' => 10**1.58, 'O9Ib' => 10**1.66, 'O9Iab' => 10**1.72, 'O9Ia' => 10**1.83,
'B0V' => 17.5, 'B0IV' => 10**1.34, 'B0III' => 10**1.40, 'B0II' => 10**1.40, 'B0Ib' => 10**1.48, 'B0Iab' => 10**1.56, 'B0Ia' => 10**1.70,
'B1V' => 11, 'B1IV' => 10**1.18, 'B1III' => 10**1.23, 'B1II' => 10**1.28, 'B1Ib' => 10**1.38, 'B1Iab' => 10**1.46, 'B1Ia' => 10**1.64,
'B2V' => 7.3, 'B2IV' => 10**1.04, 'B2III' => 10**1.08, 'B2II' => 10**1.18, 'B2Ib' => 10**1.30, 'B2Iab' => 10**1.38, 'B2Ia' => 10**1.54,
'B3V' => 5.4, 'B3IV' => 10**0.88, 'B3III' => 10**0.94, 'B3II' => 10**1.11, 'B3Ib' => 10**1.23, 'B3Iab' => 10**1.32, 'B3Ia' => 10**1.45,
'B4V' => 5.0,
'B5V' => 4.6, 'B5IV' => 10**0.72, 'B5III' => 10**0.74, 'B5II' => 10**1.00, 'B5Ib' => 10**1.18, 'B5Iab' => 10**1.26, 'B5Ia' => 10**1.40,
'B6V' => 4.0, 'B6IV' => 10**0.64, 'B6III' => 10**0.68, 'B6II' => 10**0.94, 'B6Ib' => 10**1.15, 'B6Iab' => 10**1.26, 'B6Ia' => 10**1.38,
'B7V' => 3.9, 'B7IV' => 10**0.57, 'B7III' => 10**0.60, 'B7II' => 10**0.91, 'B7Ib' => 10**1.11, 'B7Iab' => 10**1.23, 'B7Ia' => 10**1.36,
'B8V' => 3.4, 'B8IV' => 10**0.49, 'B8III' => 10**0.52, 'B8II' => 10**0.88, 'B8Ib' => 10**1.08, 'B8Iab' => 10**1.20, 'B8Ia' => 10**1.34,
'B9V' => 2.8, 'B9IV' => 10**0.45, 'B9III' => 10**0.49, 'B9II' => 10**0.85, 'B9Ib' => 10**1.04, 'B9Iab' => 10**1.20, 'B9Ia' => 10**1.32,
'A0V' => 2.3, 'A0IV' => 10**0.39, 'A0III' => 10**0.43, 'A0II' => 10**0.81, 'A0Ib' => 10**1.04, 'A0Iab' => 10**1.18, 'A0Ia' => 10**1.30,
'A1V' => 2.15, 'A1IV' => 10**0.36, 'A1III' => 10**0.41, 'A1II' => 10**0.78, 'A1Ib' => 10**1.00, 'A1Iab' => 10**1.18, 'A1Ia' => 10**1.30,
'A2V' => 2.05, 'A2IV' => 10**0.34, 'A2III' => 10**0.39, 'A2II' => 10**0.75, 'A2Ib' => 10**0.98, 'A2Iab' => 10**1.15, 'A2Ia' => 10**1.30,
'A3V' => 2.00, 'A3IV' => 10**0.32, 'A3III' => 10**0.36, 'A3II' => 10**0.75, 'A3Ib' => 10**0.97, 'A3Iab' => 10**1.11, 'A3Ia' => 10**1.30,
'A4V' => 1.90,
'A5V' => 1.85, 'A5IV' => 10**0.29, 'A5III' => 10**0.33, 'A5II' => 10**0.74, 'A5Ib' => 10**0.95, 'A5Iab' => 10**1.11, 'A5Ia' => 10**1.30,
'A6V' => 1.83,
'A7V' => 1.76, 'A7IV' => 10**0.26, 'A7III' => 10**0.30, 'A7II' => 10**0.73, 'A7Ib' => 10**0.94, 'A7Iab' => 10**1.15, 'A7Ia' => 10**1.32,
'A8V' => 1.76,
'A9V' => 1.67,
'F0V' => 1.59, 'F0IV' => 10**0.20, 'F0III' => 10**0.23, 'F0II' => 10**0.72, 'F0Ib' => 10**0.93, 'F0Iab' => 10**1.20, 'F0Ia' => 10**1.38,
'F1V' => 1.50,
'F2V' => 1.44, 'F2IV' => 10**0.16, 'F2III' => 10**0.20, 'F2II' => 10**0.72, 'F2Ib' => 10**0.93, 'F2Iab' => 10**1.20, 'F2Ia' => 10**1.40,
'F3V' => 1.43,
'F4V' => 1.39,
'F5V' => 1.33, 'F5IV' => 10**0.13, 'F5III' => 10**0.18, 'F5II' => 10**0.72, 'F5Ib' => 10**0.93, 'F5Iab' => 10**1.26, 'F5Ia' => 10**1.40,
'F6V' => 1.25,
'F7V' => 1.21,
'F8V' => 1.18, 'F8IV' => 10**0.11, 'F8II' => 10**0.72, 'F8Ib' => 10**0.93, 'F8Iab' => 10**1.28, 'F8Ia' => 10**1.41,
'F9V' => 1.14,
'G0V' => 1.08, 'G0IV' => 10**0.10, 'G0II' => 10**0.72, 'G0Ib' => 10**0.93, 'G0Iab' => 10**1.30, 'G0Ia' => 10**1.43,
'G1V' => 1.07,
'G2V' => 1.02, 'G2IV' => 10**0.10, 'G2III' => 10**0.33, 'G2II' => 10**0.72, 'G2Ib' => 10**0.93, 'G2Iab' => 10**1.30, 'G2Ia' => 10**1.45,
'G3V' => 1.00,
'G4V' => 0.99,
'G5V' => 0.98, 'G5IV' => 10**0.08, 'G5III' => 10**0.39, 'G5II' => 10**0.73, 'G5Ib' => 10**0.94, 'G5Iab' => 10**1.32, 'G5Ia' => 10**1.46,
'G6V' => 0.97,
'G7V' => 0.96,
'G8V' => 0.94, 'G8IV' => 10**0.08, 'G8III' => 10**0.42, 'G8II' => 10**0.76, 'G8Ib' => 10**0.94, 'G8Iab' => 10**1.32, 'G8Ia' => 10**1.46,
'G9V' => 0.90,
'K0V' => 0.87, 'K0IV' => 10**0.11, 'K0III' => 10**0.46, 'K0II' => 10**0.78, 'K0Ib' => 10**0.96, 'K0Iab' => 10**1.30, 'K0Ia' => 10**1.45,
'K1V' => 0.85, 'K1IV' => 10**0.13, 'K1III' => 10**0.46, 'K1II' => 10**0.78, 'K1Ib' => 10**0.96, 'K1Iab' => 10**1.30, 'K1Ia' => 10**1.45,
'K2V' => 0.78, 'K2III' => 10**0.45, 'K2II' => 10**0.79, 'K2Ib' => 10**0.98, 'K2Iab' => 10**1.28, 'K2Ia' => 10**1.43,
'K3V' => 0.75, 'K3III' => 10**0.38, 'K3II' => 10**0.80, 'K3Ib' => 10**1.00, 'K3Iab' => 10**1.30, 'K3Ia' => 10**1.43,
'K4V' => 0.72, 'K4III' => 10**0.36,
'K5V' => 0.68, 'K5III' => 10**0.37, 'K5II' => 10**0.83, 'K5Ib' => 10**1.08, 'K5Iab' => 10**1.30, 'K5Ia' => 10**1.45,
'K6V' => 0.65,
'K7V' => 0.63,
'K8V' => 0.59,
'K9V' => 0.56,
'M0V' => 0.55, 'M0III' => 10**0.48, 'M0II' => 10**0.83, 'M0Ib' => 10**1.15, 'M0Iab' => 10**1.32, 'M0Ia' => 10**1.46,
'M1V' => 0.49, 'M1III' => 10**0.54, 'M1II' => 10**0.83, 'M1Ib' => 10**1.18, 'M1Iab' => 10**1.34, 'M1Ia' => 10**1.48,
'M2V' => 0.44, 'M2III' => 10**0.54, 'M2II' => 10**0.81, 'M2Ib' => 10**1.18, 'M2Iab' => 10**1.36, 'M2Ia' => 10**1.50,
'M3V' => 0.36, 'M3III' => 10**0.52, 'M3II' => 10**0.84, 'M3Ib' => 10**1.20, 'M3Iab' => 10**1.38, 'M3Ia' => 10**1.56,
'M4V' => 0.22, 'M4III' => 10**0.51,
'M5V' => 0.16, 'M5III' => 10**0.41,
'M6V' => 0.10, 'M6III' => 10**0.40,
'M7V' => 0.090,
'M8V' => 0.082,
'M9V' => 0.079,
);

# Table matching spectral types to absolute magnitudes (for main-sequence stars)
# http://www.pas.rochester.edu/~emamajek/EEM_dwarf_UBVIJHK_colors_Teff.txt
my %SpAbsMag = (
'O3V' => -5.7,
'O4V' => -5.5,
'O5V' => -5.4,
'O6V' => -5.1,
'O7V' => -4.8,
'O8V' => -4.5,
'O9V' => -4.2,
'B0V' => -4.0,
'B1V' => -3.1,
'B2V' => -1.7,
'B3V' => -1.1,
'B4V' => -1.0,
'B5V' => -0.9,
'B6V' => -0.5,
'B7V' => -0.4,
'B8V' => -0.2,
'B9V' => 0.7,
'A0V' => 1.11,
'A1V' => 1.34,
'A2V' => 1.48,
'A3V' => 1.55,
'A4V' => 1.76,
'A5V' => 1.84,
'A6V' => 1.89,
'A7V' => 2.07,
'A8V' => 2.29,
'A9V' => 2.30,
'F0V' => 2.51,
'F1V' => 2.79,
'F2V' => 2.99,
'F3V' => 3.08,
'F4V' => 3.23,
'F5V' => 3.40,
'F6V' => 3.70,
'F7V' => 3.87,
'F8V' => 4.01,
'F9V' => 4.15,
'G0V' => 4.45,
'G1V' => 4.50,
'G2V' => 4.79,
'G3V' => 4.86,
'G4V' => 4.94,
'G5V' => 4.98,
'G6V' => 5.13,
'G7V' => 5.18,
'G8V' => 5.32,
'G9V' => 5.55,
'K0V' => 5.76,
'K1V' => 5.89,
'K2V' => 6.19,
'K3V' => 6.57,
'K4V' => 6.98,
'K5V' => 7.36,
'K6V' => 7.80,
'K7V' => 8.15,
'K8V' => 8.47,
'K9V' => 8.69,
'M0V' => 8.91,
'M1V' => 9.69,
'M2V' => 10.30,
'M3V' => 11.14,
'M4V' => 12.80,
'M5V' => 14.30,
'M6V' => 16.62,
'M7V' => 17.81,
'M8V' => 18.84,
'M9V' => 19.36,
);

my $verbose = 0; # Verbose mode will print out more information.
if (defined $ARGV[0]) {
    $verbose = 1 if ($ARGV[0] eq "-v" || $ARGV[0] eq "--verbose");
}

# ================================= READ EXISTING .STC FILES ================================= #

my @existing_bins = ();

open(VISUALBINS, '<', 'visualbins.stc') or die $!;

print "Reading default data files...\n" if $verbose;

while (my $line = <VISUALBINS>) {
    if ($line =~ /Barycenter (\d+)/) {
        push @existing_bins, $1;
    }
}

# ================================ READ INTO ASSOCIATIVE ARRAY =============================== #

my $filename = 'catalog.txt';
my $output = 'binarycatalog.stc';

my %binaries; # will store all the data for each line in the catalog
my %objects; # will store all the data for each object in Celestia

open(FILE, '<', $filename) or die $!;

print "Reading table...\n" if $verbose;

# read into associative array
while (my $line = <FILE>) {
    # skip the header
    next if $. < 2;
    chomp $line;
    
    # Stores data for each pair
    my @eachdata = split('\t', $line, -1);
    
    # convert RA, Dec from sexagesimal to degrees
    my @raarray = split(' ', $eachdata[12]);
    my $radeg = $raarray[0] * 15 + $raarray[1] * 0.25 + $raarray[2] / 240;
    my @decarray = split(' ', $eachdata[13]);
    my $decdeg = (substr($decarray[0], 1) + $decarray[1] / 60 + $decarray[2] / 3600);
    $decdeg = -$decdeg if (substr($decarray[0], 0, 1) eq "-");
    
    my $rabdeg = "";
    my $decbdeg = "";
    if ($eachdata[14]) {
        my @rabarray = split(' ', $eachdata[14]);
        $rabdeg = $rabarray[0] * 15 + $rabarray[1] * 0.25 + $rabarray[2] / 240;
    }
    if ($eachdata[15]) {
        my @decbarray = split(' ', $eachdata[15]);
        $decbdeg = (substr($decbarray[0], 1) + $decbarray[1] / 60 + $decbarray[2] / 3600);
        $decbdeg = -$decbdeg if (substr($decbarray[0], 0, 1) eq "-");
    }
    
    # Pad the index with leading zeros, so it can sort properly
    my $index = sprintf("%03d", $eachdata[0]);
    
    # Populate the hash with binaries
    $binaries{$index}{"HIP"} = $eachdata[1];
    $binaries{$index}{"HD"} = $eachdata[2];
    $binaries{$index}{"ADS"} = $eachdata[3];
    $binaries{$index}{"CCDM"} = $eachdata[4];
    $binaries{$index}{"ProperNames"} = $eachdata[5];
    $binaries{$index}{"OtherNames"} = $eachdata[6];
    $binaries{$index}{"Names1"} = $eachdata[7];
    $binaries{$index}{"Names2"} = $eachdata[8];
    $binaries{$index}{"BAY"} = $eachdata[9];
    $binaries{$index}{"FLA"} = $eachdata[10];
    $binaries{$index}{"mult"} = $eachdata[11];
    $binaries{$index}{"RA"} = $radeg;
    $binaries{$index}{"Dec"} = $decdeg;
    $binaries{$index}{"RA B"} = $rabdeg;
    $binaries{$index}{"Dec B"} = $decbdeg;
    $binaries{$index}{"plx"} = $eachdata[16];
    $binaries{$index}{"D"} = $eachdata[18];
    $binaries{$index}{"V A"} = $eachdata[20];
    $binaries{$index}{"V B"} = $eachdata[21];
    $binaries{$index}{"V type"} = $eachdata[22];
    $binaries{$index}{"SpT A"} = $eachdata[23];
    $binaries{$index}{"SpT B"} = $eachdata[24];
    $binaries{$index}{"M A"} = $eachdata[25];
    $binaries{$index}{"M B"} = $eachdata[27];
    $binaries{$index}{"R A"} = $eachdata[29];
    $binaries{$index}{"R B"} = $eachdata[31];
    $binaries{$index}{"L A"} = $eachdata[33];
    $binaries{$index}{"L B"} = $eachdata[35];
    $binaries{$index}{"L log?"} = $eachdata[37];
    $binaries{$index}{"Teff A"} = $eachdata[38];
    $binaries{$index}{"Teff B"} = $eachdata[40];
    $binaries{$index}{"T log?"} = $eachdata[42];
    $binaries{$index}{"Prot A"} = $eachdata[43];
    $binaries{$index}{"P"} = $eachdata[45];
    $binaries{$index}{"Pu"} = $eachdata[47];
    $binaries{$index}{"a"} = $eachdata[48];
    $binaries{$index}{"a0"} = $eachdata[50];
    $binaries{$index}{"au"} = $eachdata[51];
    $binaries{$index}{"e"} = $eachdata[52];
    $binaries{$index}{"i"} = $eachdata[54];
    $binaries{$index}{"node"} = $eachdata[56];
    $binaries{$index}{"arg"} = $eachdata[58];
    $binaries{$index}{"argc"} = $eachdata[60];
    $binaries{$index}{"T_0"} = $eachdata[61];
    $binaries{$index}{"Tmin"} = $eachdata[63];
    $binaries{$index}{"Tu"} = $eachdata[64];
    $binaries{$index}{"flip"} = $eachdata[65];
    $binaries{$index}{"K1"} = $eachdata[66];
}

# ===================================== PROCESS THE DATA ===================================== #

# Initialize these variables here, because they will be used across different rows
my $previousHIP = ""; # stores previously handled HIP designation
my @systemnames = (); # stores names of previously handled system
my %partialnames; # stores names that apply to part of a system, but not all (e.g. GAM2 And)
my $systemmult = "";
my $systemindex = "";
my $dist = 0;
my $newsystem = 0;
my $newHIP = 0;
my $systemscount = 0;

print "Processing data...\n" if $verbose;

# get stars from list of pairs
foreach my $index (sort keys %binaries)
{
    # Get easier-to-use variables from each hash element.
    my $HIP = $binaries{$index}{"HIP"};
    my $HD = $binaries{$index}{"HD"};
    my $ADS = $binaries{$index}{"ADS"};
    my $CCDM = $binaries{$index}{"CCDM"};
    my $ProperNames = $binaries{$index}{"ProperNames"};
    my $OtherNames = $binaries{$index}{"OtherNames"};
    my $Names1 = $binaries{$index}{"Names1"};
    my $Names2 = $binaries{$index}{"Names2"};
    my $BAY = $binaries{$index}{"BAY"};
    my $FLA = $binaries{$index}{"FLA"};
    my $mult = $binaries{$index}{"mult"};
    my $RA = $binaries{$index}{"RA"};
    my $Dec = $binaries{$index}{"Dec"};
    my $RA_B = $binaries{$index}{"RA B"};
    my $Dec_B = $binaries{$index}{"Dec B"};
    my $plx = $binaries{$index}{"plx"};
    my $D = $binaries{$index}{"D"};
    my $V_A = $binaries{$index}{"V A"};
    my $V_B = $binaries{$index}{"V B"};
    my $V_type = $binaries{$index}{"V type"};
    my $SpT_A = $binaries{$index}{"SpT A"};
    my $SpT_B = $binaries{$index}{"SpT B"};
    my $M_A = $binaries{$index}{"M A"};
    my $M_B = $binaries{$index}{"M B"};
    my $R_A = $binaries{$index}{"R A"};
    my $R_B = $binaries{$index}{"R B"};
    my $L_A = $binaries{$index}{"L A"};
    my $L_B = $binaries{$index}{"L B"};
    my $L_log = $binaries{$index}{"L log?"};
    my $Teff_A = $binaries{$index}{"Teff A"};
    my $Teff_B = $binaries{$index}{"Teff B"};
    my $T_log = $binaries{$index}{"T log?"};
    my $Prot_A = $binaries{$index}{"Prot A"};
    my $P = $binaries{$index}{"P"};
    my $Pu = $binaries{$index}{"Pu"};
    my $a = $binaries{$index}{"a"};
    my $a0 = $binaries{$index}{"a0"};
    my $au = $binaries{$index}{"au"};
    my $e = $binaries{$index}{"e"};
    my $i = $binaries{$index}{"i"};
    my $node = $binaries{$index}{"node"};
    my $arg = $binaries{$index}{"arg"};
    my $argc = $binaries{$index}{"argc"};
    my $T_0 = $binaries{$index}{"T_0"};
    my $Tmin = $binaries{$index}{"Tmin"};
    my $Tu = $binaries{$index}{"Tu"};
    my $flip = $binaries{$index}{"flip"};
    my $K1 = $binaries{$index}{"K1"};
    
    # Create list of names per system. This array is only used to check if we've moved on to
    # a new system or not.
    my @names = ();
    push @names, split(':', $ProperNames) if ($ProperNames);
    push @names, $BAY if ($BAY);
    push @names, $FLA if ($FLA);
    push @names, split(':', $OtherNames) if ($OtherNames);
    push @names, "ADS " . $ADS if ($ADS);
    push @names, "CCDM J" . $CCDM if ($CCDM);

    # To check if we've moved on to a new system, strip the component designations off of the
    # ADS/CCDM designations, and then compare the names to the current system's names.
    # If any of the names in @names are in @systemnames, then the two are part of the same
    # system.
    my @genericsysnames = StripComponents(@systemnames);
    my @newdesignation = StripComponents(@names);
    if (scalar intersect(@genericsysnames, @newdesignation) == 0) {
        $newsystem = 1;
        $systemindex = $index;
    } else {
        $newsystem = 0;
    }
    
    # Check if the multiplicity string is the same as the system's multiplicity string. If so,
    # then better to treat this as a new system.
    if (!$newsystem && $systemmult eq $mult) {
        $newsystem = 1;
        $systemindex = $index;
    }
    
    # Set the system name.
    @systemnames = @names if ($newsystem);
    
    # The Names1 Names2 columns contain a variety of different designations that have different
    # priorities, so sort them into the proper bin.
    my $HIP_AB = "";
    my ($HIP_A, $HIP_B) = ("", "");

    my @PrimaryNames = ();
    my @PrimaryDesignations = ();
    my @SecondaryNames = ();
    my @SecondaryDesignations = ();
    my @Names1 = split(':', $Names1);
    for (my $i = 0; $i < scalar @Names1; $i++) {
        if ($Names1[$i] =~ /^(HD|SAO)/) {
            push @PrimaryDesignations, $Names1[$i];
        } elsif ($Names1[$i] =~ /TYC (\d+)-(\d+)-(\d+)/) {
            $HIP_A = $3 * 1000000000 + $2 * 10000 + $1;
        } else {
            push @PrimaryNames, $Names1[$i];
        }
    }
    my @Names2 = split(':', $Names2);
    for (my $i = 0; $i < scalar @Names2; $i++) {
        if ($Names2[$i] =~ /^(HD|SAO)/) {
            push @SecondaryDesignations, $Names2[$i];
        } elsif ($Names2[$i] =~ /HIP (\d+)/) {
            $HIP_B = $1;
        } elsif ($Names2[$i] =~ /TYC (\d+)-(\d+)-(\d+)/) {
            $HIP_B = $3 * 1000000000 + $2 * 10000 + $1;
        } else {
            push @SecondaryNames, $Names2[$i];
        }
    }
    
    # Now create list of names per barycenter, primary, and secondary
    my @names_AB = ();
    my @names_A = ();
    my @names_B = ();
    my ($primary, $secondary) = GetComponents($mult, $mult);
    if ($ProperNames) {
        my @ProperNames = split(':', $ProperNames);
        for (my $i = 0; $i < scalar @ProperNames; $i++) {
            # Check each name to see if it applies to the sub-component only.
            # If it does, then the primary and secondary designations should
            # be A and B.
            # Note that currently, this has no way to add sub-components
            # (e.g. Aa, Ab).
            if (grep(/$ProperNames[$i]/, @systemnames)) {
                if (($mult eq "AB" && $newsystem) || $mult =~ /(.+)-(.+)/) {
                    push @names_AB, $ProperNames[$i];
                } else {
                    push @names_AB, $ProperNames[$i] . " " . $mult;
                }
                push @names_A, $ProperNames[$i] . " " . $primary;
                push @names_B, $ProperNames[$i] . " " . $secondary;
            } elsif ($partialnames{$ProperNames[$i]}) {
                # Compare current mult to immediately higher mult
                my $highermult = $partialnames{$ProperNames[$i]};
                my ($curprim, $cursec) = GetComponents($highermult, $highermult);
                if ($mult eq $curprim) {
                    push @names_AB, $ProperNames[$i] . " " . "A";
                    push @names_A, $ProperNames[$i] . " " . "Aa";
                    push @names_B, $ProperNames[$i] . " " . "Ab";
                } elsif ($mult eq $cursec) {
                    push @names_AB, $ProperNames[$i] . " " . "B";
                    push @names_A, $ProperNames[$i] . " " . "Ba";
                    push @names_B, $ProperNames[$i] . " " . "Bb";
                }
            } else {
                # If this name doesn't apply to the highest level, record it in this hash
                # so it can be remembered later
                $partialnames{$ProperNames[$i]} = $mult;
                push @names_AB, $ProperNames[$i];
                push @names_A, $ProperNames[$i] . " " . "A";
                push @names_B, $ProperNames[$i] . " " . "B";
            }
        }
    }
    for (my $i = 0; $i < scalar @PrimaryNames; $i++) {
        push @names_A, $PrimaryNames[$i];
    }
    for (my $i = 0; $i < scalar @SecondaryNames; $i++) {
        push @names_B, $SecondaryNames[$i];
    }
    if ($BAY) {
        if (grep(/$BAY/, @systemnames)) {
            if (($mult eq "AB" && $newsystem) || $mult =~ /(.+)-(.+)/) {
                push @names_AB, $BAY;
            } else {
                push @names_AB, $BAY . " " . $mult;
            }
            if ($BAY =~ / (A|B|C)$/) {
                push @names_A, $BAY . "a";
                push @names_B, $BAY . "b";
            } else {
                push @names_A, $BAY . " " . $primary;
                push @names_B, $BAY . " " . $secondary;
            }
        } elsif ($partialnames{$BAY}) {
            # Compare current mult to immediately higher mult
            my $highermult = $partialnames{$BAY};
            my ($curprim, $cursec) = GetComponents($highermult, $highermult);
            if ($mult eq $curprim) {
                push @names_AB, $BAY . " " . "A";
                push @names_A, $BAY . " " . "Aa";
                push @names_B, $BAY . " " . "Ab";
            } elsif ($mult eq $cursec) {
                push @names_AB, $BAY . " " . "B";
                push @names_A, $BAY . " " . "Ba";
                push @names_B, $BAY . " " . "Bb";
            }
        } else {
            $partialnames{$BAY} = $mult;
            push @names_AB, $BAY;
            push @names_A, $BAY . " " . "A";
            push @names_B, $BAY . " " . "B";
        }
    }
    if ($FLA) {
        if (grep(/$FLA/, @systemnames)) {
            if (($mult eq "AB" && $newsystem) || $mult =~ /(.+)-(.+)/) {
                push @names_AB, $FLA;
            } else {
                push @names_AB, $FLA . " " . $mult;
            }
            if ($FLA =~ / (A|B|C)$/) {
                push @names_A, $FLA . "a";
                push @names_B, $FLA . "b";
            } else {
                push @names_A, $FLA . " " . $primary;
                push @names_B, $FLA . " " . $secondary;
            }
        } elsif ($partialnames{$FLA}) {
            # Compare current mult to immediately higher mult
            my $highermult = $partialnames{$FLA};
            my ($curprim, $cursec) = GetComponents($highermult, $highermult);
            if ($mult eq $curprim) {
                push @names_AB, $FLA . " " . "A";
                push @names_A, $FLA . " " . "Aa";
                push @names_B, $FLA . " " . "Ab";
            } elsif ($mult eq $cursec) {
                push @names_AB, $FLA . " " . "B";
                push @names_A, $FLA . " " . "Ba";
                push @names_B, $FLA . " " . "Bb";
            }
        } else {
            $partialnames{$FLA} = $mult;
            push @names_AB, $FLA;
            push @names_A, $FLA . " " . "A";
            push @names_B, $FLA . " " . "B";
        }
    }
    if ($OtherNames) {
        my @LaterNames = split(':', $OtherNames);
        for (my $i = 0; $i < scalar @LaterNames; $i++) {
            if (grep(/$LaterNames[$i]/, @systemnames)) {
                if (($mult eq "AB" && $newsystem) || $mult =~ /(.+)-(.+)/) {
                    push @names_AB, $LaterNames[$i];
                } else {
                    push @names_AB, $LaterNames[$i] . " " . $mult;
                }
                push @names_A, $LaterNames[$i] . " " . $primary;
                push @names_B, $LaterNames[$i] . " " . $secondary;
            } elsif ($partialnames{$LaterNames[$i]}) {
                # Compare current mult to immediately higher mult
                my $highermult = $partialnames{$LaterNames[$i]};
                my ($curprim, $cursec) = GetComponents($highermult, $highermult);
                if ($mult eq $curprim) {
                    push @names_AB, $LaterNames[$i] . " " . "A";
                    push @names_A, $LaterNames[$i] . " " . "Aa";
                    push @names_B, $LaterNames[$i] . " " . "Ab";
                } elsif ($mult eq $cursec) {
                    push @names_AB, $LaterNames[$i] . " " . "B";
                    push @names_A, $LaterNames[$i] . " " . "Ba";
                    push @names_B, $LaterNames[$i] . " " . "Bb";
                }
            } else {
                $partialnames{$LaterNames[$i]} = $mult;
                push @names_AB, $LaterNames[$i];
                push @names_A, $LaterNames[$i] . " " . "A";
                push @names_B, $LaterNames[$i] . " " . "B";
            }
        }
    }
    if ($ADS) {
        if ($ADS =~ /^(\d+) (\w+)?$/) {
            ($primary, $secondary) = GetComponents($2, $mult);
        }
        push @names_AB, "ADS " . $ADS;
        push @names_A, "ADS " . $1 . " " . $primary;
        push @names_B, "ADS " . $1 . " " . $secondary;
    }
    if ($CCDM) {
        if ($CCDM =~ /^(\d+[-\+]\d+)([A-Za-z]+)+?$/) {
            ($primary, $secondary) = GetComponents($2, $mult);
        }
        push @names_AB, "CCDM J" . $CCDM;
        push @names_A, "CCDM J" . $1 . $primary;
        push @names_B, "CCDM J" . $1 . $secondary;
    }
    for (my $i = 0; $i < scalar @PrimaryDesignations; $i++) {
        push @names_A, $PrimaryDesignations[$i];
    }
    for (my $i = 0; $i < scalar @SecondaryDesignations; $i++) {
        push @names_B, $SecondaryDesignations[$i];
    }
    
    # print "$names_AB[0]\n";
    # print "@names_A\n";
    # print "@names_B\n";
    
    # In the catalog, if the spectral type is listed as "*", then it is a barycenter
    my ($isabarycenter, $isbbarycenter) = (0, 0);
    $isabarycenter = 1 if ($SpT_A eq "*");
    $isbbarycenter = 1 if ($SpT_B eq "*");
    
    # Get distances - distance column is preferred, then use parallax column.
    # We are using the "naive" 1/parallax method, since the errors are (mostly) low.
    # If the distance or parallax is not given in the column, use the last successful value.
    # TODO: provide safeguard that it's not using previous system's distance.
    if ($D) {
        $dist = $D * $LY_TO_PARSEC;
    } elsif ($plx) {
        $dist = $LY_TO_PARSEC * (1000 / $plx);
    }
    
    # First, convert from logarithms to actual values, if necessary.
    if ($L_log eq "log") {
        $L_A = 10**$L_A if $L_A;
        $L_B = 10**$L_B if $L_B;
    }
    if ($T_log eq "log") {
        $Teff_A = 10**$Teff_A if $Teff_A;
        $Teff_B = 10**$Teff_B if $Teff_B;
    }
    
    # Trim spectral types with dashes or slashes, so that luminosity classes get read by Celestia
    # For simplicity, the earliest spectral type in the string is used.
    $SpT_A = TrimSpType($SpT_A);
    $SpT_B = TrimSpType($SpT_B);
    # Estimate spectral types if missing.
    # Use the temperature, then the mass (assuming main-sequence), then the magnitude.
    my ($SpT_ANote, $SpT_BNote) = ("", "");
    if (!$SpT_A) {
        if ($Teff_A) {
            $SpT_A = Teff_to_SpT($Teff_A);
            $SpT_ANote = " # Estimate from temperature";
        } elsif ($M_A) {
            $SpT_A = Mass_to_SpT($M_A);
            $SpT_ANote = " # Estimate from mass";
        } elsif ($V_A && $V_type eq "MV") {
            $SpT_A = AbsMag_to_SpT($V_A);
            $SpT_ANote = " # Estimate from absolute magnitude";
        } elsif ($V_A && $V_type eq "mV") {
            my $absmag = App_to_AbsMag($V_A, $dist);
            $SpT_A = AbsMag_to_SpT($absmag);
            $SpT_ANote = " # Estimate from apparent magnitude";
        }
    }
    if (!$SpT_B) {
        if ($Teff_B) {
            $SpT_B = Teff_to_SpT($Teff_B);
            $SpT_BNote = " # Estimate from temperature";
        } elsif ($M_B) {
            $SpT_B = Mass_to_SpT($M_B);
            $SpT_BNote = " # Estimate from mass";
        } elsif ($V_B && $V_type eq "MV") {
            $SpT_B = AbsMag_to_SpT($V_B);
            $SpT_BNote = " # Estimate from absolute magnitude";
        } elsif ($V_B && $V_type eq "mV") {
            my $absmag = App_to_AbsMag($V_B, $dist);
            $SpT_B = AbsMag_to_SpT($absmag);
            $SpT_BNote = " # Estimate from apparent magnitude";
        }
    }
    
    # Extract magnitudes and estimate values, if missing.
    my ($AppMagA, $AppMagB) = ("", "");
    my ($AbsMagA, $AbsMagB) = ("", "");
    my ($Mag_ANote, $Mag_BNote) = ("", "");
    if ($V_A && $V_type eq "mV") {
        $AppMagA = $V_A;
    } elsif ($V_A && $V_type eq "MV") {
        $AbsMagA = $V_A;
    } elsif (($L_A || $R_A) && $Teff_A) {
        # Calculate using Stefan-Boltzmann law, and bolometric correction
        # Mass-radius relationship doesn't work for white dwarfs, so we exclude them
        if (!$L_A) {
            $L_A = SBL_to_Lum($R_A, $Teff_A);
            $Mag_ANote = " # From radius and temperature";
        } else {
            $Mag_ANote = " # From luminosity and temperature";
        }
        my $MBolA = 4.74 - 2.5 * log10($L_A);
        my $BC = Teff_to_BC($Teff_A);
        $AbsMagA = $MBolA - $BC;
    } elsif ($SpT_A =~ /([OBAFGKM])([0-9])+[\/.-]?[0-9]?([-V*]+)/ && !$SpT_ANote) {
        # Check if spectral type corresponds to a main-sequence star with subclass
        # If so, use it to guess the magnitude
        my $SpTKey = substr($1, 0, 1) . substr($2, 0, 1) . substr($3, 0, 1);
        $AbsMagA = $SpAbsMag{$SpTKey};
        if ($SpT_ANote eq " # Estimate from mass") {
            $Mag_ANote = " # Estimate from mass";
        } else {
            $Mag_ANote = " # Estimate from spectral type";
        }
    } elsif ($SpT_A =~ /^D/) {
        # TODO: estimate magnitudes for white dwarfs
        $AbsMagA = 13;
        $Mag_ANote = " # Guess, for a white dwarf";
    } elsif ($M_A) {
        # Assume that spectral types are main-sequence or close to it
        # (white dwarfs taken care of), and estimate mass
        my $SpTKey = Mass_to_SpT($M_A);
        $AbsMagA = $SpAbsMag{$SpTKey};
        $Mag_ANote = " # Guess, from mass";
    }
    if ($V_B && $V_type eq "mV") {
        $AppMagB = $V_B;
    } elsif ($V_B && $V_type eq "MV") {
        $AbsMagB = $V_B;
    } elsif (($L_B || $R_B) && $Teff_B) {
        # Calculate using Stefan-Boltzmann law, and bolometric correction
        # Mass-radius relationship doesn't work for white dwarfs, so we exclude them
        if (!$L_B) {
            $L_B = SBL_to_Lum($R_B, $Teff_B);
            $Mag_BNote = " # From radius and temperature";
        } else {
            $Mag_BNote = " # From luminosity and temperature";
        }
        my $MBolB = 4.74 - 2.5 * log10($L_B);
        my $BC = Teff_to_BC($Teff_B);
        $AbsMagB = $MBolB - $BC;
    } elsif ($SpT_B =~ /([OBAFGKM])([0-9])+[\/.-]?[0-9]?([-V*]+)/ && !$SpT_BNote) {
        # Check if spectral type corresponds to a main-sequence star with subclass
        # If so, use it to guess the magnitude
        my $SpTKey = substr($1, 0, 1) . substr($2, 0, 1) . substr($3, 0, 1);
        $AbsMagB = $SpAbsMag{$SpTKey};
        if ($SpT_BNote eq " # Estimate from mass") {
            $Mag_BNote = " # Estimate from mass";
        } else {
            $Mag_BNote = " # Estimate from spectral type";
        }
    } elsif ($SpT_B =~ /^D/) {
        # TODO: estimate magnitudes for white dwarfs
        $AbsMagB = 13;
        $Mag_BNote = " # Guess, for a white dwarf";
    } elsif ($M_B) {
        # Assume that spectral types are main-sequence or close to it
        # (white dwarfs taken care of), and estimate mass
        my $SpTKey = Mass_to_SpT($M_B);
        $AbsMagB = $SpAbsMag{$SpTKey};
        $Mag_BNote = " # Guess, from mass";
    }
    
    if ($HIP ne $previousHIP) {
        $newHIP = 1;
        $HIP_AB = $HIP;
        
        # Case for if the there are two HIP designations
        if ($HIP =~ /(\d+)\/(\d+)/) {
            $HIP_AB = "";
            $HIP_A = $1;
            $HIP_B = substr($1, 0, -length($2)) . $2;
        }
        # If the primary and secondary have their own TYC designations, then apply the HIP
        # designation to the barycenter. Otherwise, if the secondary has a HD designation
        # but the first does not, then the HIP designation should be on the primary to
        # prevent the secondary from effectively receiving two HIP designations.
        if (!($HIP_A && $HIP_B) && !grep(/HD/, @names_A) && grep(/HD/, @names_B)) {
            $HIP_AB = "";
            $HIP_A = $HIP;
        }
        # If only the secondary has a TYC designation, then apply the HIP designation to
        # the primary.
        if ($HIP && $HIP_B && !$HIP_A) {
            $HIP_AB = "";
            $HIP_A = $HIP;
        }
    } else {
        $newHIP = 0;
    }
    
    # Convert radii to km, for use in Celestia.
    my ($RadiusA, $RadiusB) = ("", "");
    $RadiusA = FormatSigFigs($R_A * $SOLAR_RADIUS, 4) if ($R_A);
    $RadiusB = FormatSigFigs($R_B * $SOLAR_RADIUS, 4) if ($R_B);
    
    # Estimate masses, using spectral types as a guess.
    # Mass ratios are approximate because Straizys & Kuriliene's table is very old, and for giant
    # stars, the masses may vary for the same spectral type.
    my ($M_ANote, $M_BNote) = ("", "");
    if (!$M_A && $SpT_A) {
        $M_A = SpT_to_Mass($SpT_A);
        $M_ANote = " # Guess";
    }
    if (!$M_B && $SpT_B) {
        $M_B = SpT_to_Mass($SpT_B);
        $M_BNote = " # Guess";
    }
    
    # Write notes containing the mass ratios. Barycenters are marked with asterisks for now; the
    # component masses are added later in the script.
    my $AxisNote = "";
    my ($PrintedM_A, $PrintedM_B) = ($M_A, $M_B);
    $PrintedM_A = "*" if ($isabarycenter);
    $PrintedM_B = "*" if ($isbbarycenter);
    if ($M_ANote || $M_BNote) {
        $AxisNote = " # Mass ratios approximate";
    } else {
        $AxisNote = " # Mass ratio $PrintedM_A:$PrintedM_B";
    }
    my $AxisNote_AB = ""; # Mass ratio comments, when on barycenters, are handled differently.

    # Initialize variables that are used for "positioning". Most are "positioned" using orbits,
    # but some are just placed as stationary stars.
    # Parameters in Celestia:
    my $period = $P;
    my $Epoch = "";
    my $MeanAnomaly = "";
    my $Atot = "";
    my ($A1, $A2) = ("", "");
    my ($Arg1, $Arg2) = ("", "");
    my ($Inclination, $AscendingNode, $QuotedPeri) = ("", "", "");
    
    # Comments that appear next to the parameters:
    my $missingparam = " # Plane-of-sky parameter unknown";
    my $InclinationNote = "";
    my $NodeNote = "";
    my $ArgNote = "";
    my $EpochNote = "";
    my $FullySpecified = "";
    
    # RA and Dec values for stationary stars.
    my $RA_AB = $RA;
    my $Dec_AB = $Dec;
    my $RA_A = $RA;
    my $Dec_A = $Dec;

    # If an orbit is defined, calculate the parameters for Celestia.
    if ($P) {
        # Standardize periods to use years.
        $period = $P / 365.25 if ($Pu eq "d");
        
        # If T_0 is given in Julian date, use the Epoch parameter in Celestia. Otherwise, use the
        # MeanAnomaly parameter.
        if ($T_0) {
            if ($Tu eq "JD") {
                $Epoch = $T_0;
            } elsif ($Tu eq "B") {
                $MeanAnomaly = mod((2000 - $T_0) / $period, 1) * 360;
            }
        }
        
        # Convert angular semimajor axes to physical values, taking into account the orbit type.
        if (!$a) {
            # If no semimajor axis, calculate it using Kepler's third law
            $Atot = (($period**2) * ($M_A+$M_B))**(1/3);
            $A1 = $Atot * $M_B / ($M_A + $M_B);
            $A2 = $Atot * $M_A / ($M_A + $M_B);
        } elsif (!$a0) {
            if ($au eq "a") {
                $Atot = $a;
            } elsif ($au eq "r") {
                $Atot = $a / 215.032;
            } elsif ($au eq "s") {
                $Atot = Arcsec_to_AU($a, $dist);
            } elsif ($au eq "m") {
                $Atot = Arcsec_to_AU($a / 1000, $dist);
            }
            $A1 = $Atot * $M_B / ($M_A + $M_B);
            $A2 = $Atot * $M_A / ($M_A + $M_B);
        } elsif ($a0) {
            # If the a0 column is marked "a0", then the quoted value refers to the orbit of A
            # around the barycenter (aka its photocentric orbit). This is typically the case
            # for astrometric binaries.
            if ($au eq "a") {
                $A1 = $a;
            } elsif ($au eq "r") {
                $A1 = $a / 215.032;
            } elsif ($au eq "s") {
                $A1 = Arcsec_to_AU($a, $dist);
            } elsif ($au eq "m") {
                $A1 = Arcsec_to_AU($a/1000, $dist);
            }
            $A2 = $A1 / ($M_B / $M_A);
            $Atot = $A1 + $A2;
        }
        
        # Estimate inclination for spectroscopic binaries without visual/astrometric detection.
        # This is done using the velocity semi-amplitude of the primary (K1), period, and eccentricity
        # to calculate the a*sin(i). The true semi-major axis is already estimated, using component
        # masses and Kepler's third law.
        # Equation 5 from Bischoff et al. (2017), AN Vol.338, Issue 6, pg.671
        # "Radial velocity measurements and orbit determination of eight single-lined spectroscopic
        # binary systems"
        if (!$i && $K1) {
            my $asini = $K1 * ($period * 365.25 * 86400) * sqrt(1 - $e**2) / (2 * pi) / 149597870.7;
            if ($asini <= $A1) {
                $i = sprintf("%u", rad2deg(asin($asini / $A1)));
                $InclinationNote = "\n\t\t# Plane-of-sky inclination is about $i degrees, estimated from K1";
            } else {
                # If K1 value is missing, assume 45 degrees for a spectroscopic binary.
                # An inclination of 90 degrees would cause eclipses, while an inclination of 0
                # degrees would be undetectable.
                $i = 45;
                $InclinationNote = $missingparam;
            }
        }
        
        # Because different astronomers use different conventions (i.e. there is a degeneracy between
        # i and -i), some of these (marked in the "flip" column) have to be modified to work for
        # Celestia.
        $i = -$i if ($flip eq "yes");
        
        # Transform orbits.
        ($Inclination, $AscendingNode, $QuotedPeri) = OrbitTransform($RA, $Dec, $period, $i, $node, $arg);
        
        # Apply argument of periastron to correct components, using the argc column.
        # Usually the quoted value refers to the secondary (i.e. argc is flagged as 2), but for
        # astrometric orbits and single-lined spectroscopic orbits, the value applies to the primary
        # (the one that is actually detected.) These are flagged as 1 in the argc column.
        if (!$argc) {
            if ($T_0) {
                # If there is an epoch, an argument of periastron value is still needed to place
                # the objects in the correct place at the correct time.
                $Arg2 = $QuotedPeri;
                $Arg1 = mod(($QuotedPeri + 180), 360);
            } else {
                $Arg2 = 180; # for cases where no argument of periastron is given
            }
        } elsif ($argc == 1) {
            $Arg1 = $QuotedPeri;
            $Arg2 = mod(($QuotedPeri + 180), 360);
        } elsif ($argc == 2) {
            $Arg2 = $QuotedPeri;
            $Arg1 = mod(($QuotedPeri + 180), 360);
        }
        
        # Make a note if the plane-of-sky parameter is not given. (The $InclinationNote variable is
        # used earlier, to keep note if the plane-of-sky inclination is estimated from K1.)
        $InclinationNote = $missingparam if (!$i);
        $NodeNote = $missingparam if (!$node);
        $ArgNote = $missingparam if ($arg eq "");
        $ArgNote = "" if ($Arg2 == 180);
        
        # If $Tmin is marked as "Tmin", then the quoted epoch refers to the epoch of primary minimum.
        if ($Tmin) {
            # Shortcut: if circular orbit, then the MeanAnomaly value is 90.
            if (!$e) {
                $MeanAnomaly = 90;
                $ArgNote = "";
            } else {
                my $longperi = $arg;
                $longperi = 0 if (!$arg);
                my $trueanomaly = deg2rad(mod(90 - $longperi, 360));
                my $eanomaly = acos((cos($trueanomaly)+$e) / (1+$e*cos($trueanomaly)));
                my $meananomaly = rad2deg($eanomaly - $e * sin($eanomaly));
                $meananomaly = 360 - $meananomaly if ($trueanomaly > deg2rad(180));
                $MeanAnomaly = $meananomaly;
            }
            $EpochNote = " # Epoch of primary minimum";
        }
        
        # Format orbital elements using the correct (or at least reasonable) number of sig figs.
        $period = FormatSigFigs($period, CountSigFigs($P));
        my $asigfigs = 3;
        $asigfigs = CountSigFigs($a) if ($a);
        $A1 = FormatSigFigs($A1, $asigfigs);
        $A2 = FormatSigFigs($A2, $asigfigs);
        my $idigits = 2;
        $idigits = CountDecimalDigits($i) if ($i);
        $Inclination = sprintf("%.*f", $idigits, $Inclination);
        $AscendingNode = sprintf("%.*f", $idigits, $AscendingNode);
        my $argdigits = 0;
        $argdigits = CountDecimalDigits($arg) if ($arg);
        $argdigits = $idigits if (!$arg && $T_0);
        
        $Arg1 = sprintf("%.*f", $argdigits, $Arg1) if ($Arg1);
        $Arg2 = sprintf("%.*f", $argdigits, $Arg2) if ($Arg2);
        $MeanAnomaly = sprintf("%.*f", $idigits, $MeanAnomaly) if ($MeanAnomaly);
        
        # If no epoch value, then don't use it
        if (!$T_0) {
            $MeanAnomaly = "";
            $Epoch = "";
        }
        
        # If there is no epoch, inclination, or longitude of ascending node, then don't bother with
        # positioning at all
        if (!$T_0 && !$i && !$node) {
            $Inclination = "";
            $AscendingNode = "";
        }
        
        # Mark orbits that are fully specified (all orbital elements known, and radial velocity
        # resolves mirror ambiguity)
        $FullySpecified = "       # Fully specified orientation" if ($node && $flip);

    # Otherwise, calculate RA and Dec for each component.
    } elsif ($RA_B) {
        $RA_A = $RA;
        $Dec_A = $Dec;
        $RA_AB = (($RA * $M_A) + ($RA_B * $M_B)) / ($M_A + $M_B);
        $Dec_AB = (($Dec * $M_A) + ($Dec_B * $M_B)) / ($M_A + $M_B);
        
        # Put mass ratio comments on the barycenter.
        $AxisNote_AB = $AxisNote;
        $AxisNote = "";
    } else {
        $RA_B = $RA;
        $Dec_B = $Dec;
    }
    
    # Get rotation period and convert from days to hours.
    my $RotationPeriod = "";
    $RotationPeriod = $Prot_A * 24 if ($Prot_A);
    
    # Check to see if a system is already in a default .stc file.
    my $duplicate = 0;
    if (grep(/$HIP/, @existing_bins) && $mult eq "AB") {
        $duplicate = 1;
    }
    
    # The key for %objects is a five-digits key. It doesn't mean anything, but is used
    # internally to sort all the objects for the .stc file.
    #
    # Redefine $primary, $secondary to use the mult column, not the ADS/CCDM names.
    ($primary, $secondary) = GetComponents($mult, $mult);
    
    # Create new entries for %objects. First, the barycenter:
    if (!exists $objects{$systemindex . $csort{$mult}}) {
        $objects{$systemindex . $csort{$mult}} = {
            'Names' => [@names_AB],
            'Barycenter' => 1,
            'RA' => $RA_AB,
            'Dec' => $Dec_AB,
            'HIP' => $HIP_AB,
            'Parent' => '',
            'Distance' => $dist,
            'mult' => $mult,
            'AxisNote' => $AxisNote_AB,
            'Duplicate' => $duplicate,
        };
    } else {
        # Add new HIP designation, if it exists
        $objects{$systemindex . $csort{$mult}}{'HIP'} = $HIP_AB if (!$HIP_A && $HIP_AB);
        $objects{$systemindex . $csort{$mult}}{'HIP'} = $HIP_AB if ($HIP_A && $HIP_B && $HIP_AB);
        
        my $curnames = $objects{$systemindex . $csort{$mult}}{'Names'};
        # Decide whether or not the names should be replaced. Go with the array with more
        # entries. If they have the same number of entries, replace it.
        if (scalar @$curnames <= scalar @names_AB) {
            $objects{$systemindex . $csort{$mult}}{'Names'} = [@names_AB];
        }
        
        # Add component masses to higher-level pair
        my $currentmult = $objects{$systemindex . $csort{$mult}}{'mult'};
        my ($curprim, $cursec) = GetComponents($currentmult, $currentmult);
        
        my $MassRatioInsert = "*";
        if ($AxisNote =~ /Mass ratio ([\*\(\d].*[\*\(\d])/) {
            $MassRatioInsert = "(" . $1 . ")";
        }
        $objects{$systemindex . $csort{$currentmult}}{'AxisNote'} =~ s/\*/$MassRatioInsert/;
        $objects{$systemindex . $csort{$curprim}}{'AxisNote'} =~ s/\*/$MassRatioInsert/;
        $objects{$systemindex . $csort{$cursec}}{'AxisNote'} =~ s/\*/$MassRatioInsert/;
        
        # Use mult from this level to clarify higher-level pair, for ADS and CCDM
        my $ADS_child = "none"; # it cannot be "", because this is what we are looking for
        my $CCDM_child = "none";
        if ($ADS =~ /\d+ (.+)/) {
            $ADS_child = $1;
        }
        if ($CCDM =~ /\d+(\+|-)\d+([A-Z]*)/) {
            $CCDM_child = $2;
        }
        my $currentmultnames = $objects{$systemindex . $csort{$currentmult}}{'Names'};
        my $ADS_parent = "none"; # it cannot be "", because this is what we are looking for
        my $CCDM_parent = "none";
        for (my $i = 0; $i < scalar @$currentmultnames; $i++) {
            if ($$currentmultnames[$i] =~ /ADS \d+ (.+)/) {
                $ADS_parent = $1;
            } elsif ($$currentmultnames[$i] =~ /CCDM J\d+(\+|-)\d+([A-Z]*)/) {
                $CCDM_parent = $2;
            }
        }
        my $ADS_sibling = GetSibling($ADS_parent, $ADS_child);
        my $CCDM_sibling = GetSibling($CCDM_parent, $CCDM_child);
        my $curprimnames = $objects{$systemindex . $csort{$curprim}}{'Names'};
        my $cursecnames = $objects{$systemindex . $csort{$curprim}}{'Names'};
        for (my $i = 0; $i < scalar @$curprimnames; $i++) {
            if ($$curprimnames[$i] =~ /ADS \d+ $/) {
                $$curprimnames[$i] = $$curprimnames[$i] . $ADS_sibling;
            } elsif ($$curprimnames[$i] =~ /CCDM J\d+(\+|-)\d+$/) {
                $$curprimnames[$i] = $$curprimnames[$i] . $CCDM_sibling;
            }
        }
        for (my $i = 0; $i < scalar @$cursecnames; $i++) {
            if ($$cursecnames[$i] =~ /ADS \d+ $/) {
                $$cursecnames[$i] = $$cursecnames[$i] . $ADS_sibling;
            } elsif ($$cursecnames[$i] =~ /CCDM J\d+(\+|-)\d+$/) {
                $$cursecnames[$i] = $$cursecnames[$i] . $CCDM_sibling;
            }
        }
        
        # If there is a case where both AB and Ab exist, change AB to A-B. This is because
        # Celestia cannot differentiate between the two, and so it may put the stars in orbit
        # around the wrong object.
        if ($mult eq "A" && $currentmult eq "AB") {
            my $badbarycenter = $objects{$systemindex . $csort{$curprim}}{'Parent'};
            if (lc($names_AB[0] . "b") eq lc($badbarycenter)) {
                $badbarycenter =~ s/AB/A-B/;
                $objects{$systemindex . $csort{$curprim}}{'Parent'} = $badbarycenter;
                $objects{$systemindex . $csort{$cursec}}{'Parent'} = $badbarycenter;
                my $badnames = $objects{$systemindex . $csort{$currentmult}}{'Names'};
                for (my $i = 0; $i < scalar @$badnames; $i++) {
                    $$badnames[$i] =~ s/AB/A-B/ if ($$badnames[$i] !~ /^(ADS|CCDM)/);
                }
            }
        }

        # Keep going up until all of the component masses are filled in.
        # Currently, only one iteration is needed and it applies to only one system (XI Tau)
        # but it may be better to use a while loop for this.
        # First, check if there's a barycenter "higher up" by using the "Parent" value.
        if ($objects{$systemindex . $csort{$currentmult}}{'Parent'}) {
            my $parentname = $objects{$systemindex . $csort{$currentmult}}{'Parent'};
            
            # Go through each of the keys and find out which one it is, so we can extract
            # the "mult" parameter.
            foreach my $key (sort keys %objects) {
                my $namearray = $objects{$key}{'Names'};
                if (@$namearray[0] eq $parentname) {
                    $currentmult = $objects{$key}{'mult'};
                    ($curprim, $cursec) = GetComponents($currentmult, $currentmult);
                    $objects{$systemindex . $csort{$currentmult}}{'AxisNote'} =~ s/\*/$MassRatioInsert/;
                    $objects{$systemindex . $csort{$curprim}}{'AxisNote'} =~ s/\*/$MassRatioInsert/;
                    $objects{$systemindex . $csort{$cursec}}{'AxisNote'} =~ s/\*/$MassRatioInsert/;
                    
                    last;
                }
            }
        }
    }

    # Next, the primary:
    $objects{$systemindex . $csort{$primary}} = {
        'Names' => [@names_A],
        'Barycenter' => $isabarycenter,
        'RA' => $RA_A,
        'Dec' => $Dec_A,
        'HIP' => $HIP_A,
        'Parent' => $names_AB[0],
        'Distance' => $dist,
        'SpectralType' => $SpT_A,
        'SpectralTypeNote' => $SpT_ANote,
        'AppMag' => $AppMagA,
        'AbsMag' => $AbsMagA,
        'MagNote' => $Mag_ANote,
        'Radius' => $RadiusA,
        'Mass' => $M_A,
        'MassNote' => $M_ANote,
        'AxisNote' => $AxisNote,
        'InclinationNote' => $InclinationNote,
        'NodeNote' => $NodeNote,
        'ArgNote' => $ArgNote,
        'EpochNote' => $EpochNote,
        'FullySpecified' => $FullySpecified,
        'Period' => $period,
        'SemiMajorAxis' => $A1,
        'Eccentricity' => $e,
        'Inclination' => $Inclination,
        'AscendingNode' => $AscendingNode,
        'ArgOfPericenter' => $Arg1,
        'MeanAnomaly' => $MeanAnomaly,
        'Epoch' => $Epoch,
        'RotationPeriod' => $RotationPeriod,
        'mult' => $mult,
        'Duplicate' => $duplicate,
    };
    
    # Finally, the secondary:
    $objects{$systemindex . $csort{$secondary}} = {
        'Names' => [@names_B],
        'Barycenter' => $isbbarycenter,
        'RA' => $RA_B,
        'Dec' => $Dec_B,
        'HIP' => $HIP_B,
        'Parent' => $names_AB[0],
        'Distance' => $dist,
        'SpectralType' => $SpT_B,
        'SpectralTypeNote' => $SpT_BNote,
        'AppMag' => $AppMagB,
        'AbsMag' => $AbsMagB,
        'MagNote' => $Mag_BNote,
        'Radius' => $RadiusB,
        'Mass' => $M_B,
        'MassNote' => $M_BNote,
        'AxisNote' => $AxisNote,
        'InclinationNote' => $InclinationNote,
        'NodeNote' => $NodeNote,
        'ArgNote' => $ArgNote,
        'EpochNote' => $EpochNote,
        'FullySpecified' => $FullySpecified,
        'Period' => $period,
        'SemiMajorAxis' => $A2,
        'Eccentricity' => $e,
        'Inclination' => $Inclination,
        'AscendingNode' => $AscendingNode,
        'ArgOfPericenter' => $Arg2,
        'MeanAnomaly' => $MeanAnomaly,
        'Epoch' => $Epoch,
        'mult' => $mult,
        'Duplicate' => $duplicate,
    };

    # Finally, pass name onto these variables so it can be compared with the next entry
    $previousHIP = $HIP;
    if ($newsystem) {
        $systemmult = $mult;
        # Increment number of systems
        $systemscount++ if (!$duplicate);
    }
}

# ==================================== OUTPUT THE STC FILE =================================== #

open(OUTPUT, '>', $output) or die $!;

print "Outputting .stc file...\n" if $verbose;

# Print header
print OUTPUT <<END;
# Catalogue of named multiple stars for Celestia
#
# This file was made using a Perl script and a table containing the relevant data for
# these systems. The two were custom-made for each other. The sources for the data are
# in a separate file, called references.txt.
#
# Data sources and tools used for parameter calculation:
#
# Grant Hutchison's Star Orbit Translation: starorbs.xls
# https://www.classe.cornell.edu/~seb/celestia/hutchison/spreadsheets.html#2
#
# Mamajek, E. "A Modern Mean Dwarf Stellar Color and Effective Temperature Sequence"
# http://www.pas.rochester.edu/~emamajek/EEM_dwarf_UBVIJHK_colors_Teff.txt
#
# Reed (1998), JRASC Vol.92, p.36
# "The Composite Observational-Theoretical HR Diagram"
#
# Straizys & Kuriliene (1981), Ap&SS Vol.80, p.353
# "Fundamental stellar parameters derived from the evolutionary tracks"
#
# Because information on star systems is often lacking, various relations are employed
# to fill in missing data. Missing spectral types have been estimated from temperatures,
# masses, or V-magnitudes (in that order) using Mamajek's stellar parameter table.
# Missing magnitudes have been estimated from luminosities and bolometric corrections,
# of which the latter are estimated from effective temperatures using the formula by
# Reed (1998). Otherwise, they have been estimated from the spectral type or mass using
# Mamajek's table.
#
# Stellar masses are used to calculate the semimajor axes for the stars if an angular or
# true value is missing. They are estimated from Mamajek's table, or if lacking, Table 6
# from Straizys & Kuriliene (1981) for other stars.
#
# Extrasolar orbits are typically given with the plane-of-sky used as a reference frame.
# These must be transformed into Celestia's reference frame, which is the J2000 ecliptic.
# This is done using the equations from Grant Hutchison's star orbit spreadsheet. For
# purely visual orbits, the true orientation of the orbital plane is unknown, because it
# is not possible to tell whether the periastron is pointing towards or away from the
# observer. This can only be found using radial velocity or eclipses. Systems where the
# true orientation is known are marked with the phrase "Fully specified orbit".
#
# In this file, unresolved spectroscopic binaries are not included, unless they are a
# component of a system with a visual orbit. In some cases, it is possible to estimate an
# inclination using the velocity semi-amplitude (K1). If unknown, the plane-of-sky
# inclination and ascending node have been set to 0, and are marked with the phrase
# "Plane-of-sky parameter unknown".
#
# Some binary pairs are assumed to be tidally locked. The (somewhat arbitrary) criteria
# are an orbital period less than a year, and an eccentricity less than 0.01.
#
# HIP indices are generally applied to the barycenter. This is so that their components
# receive their respective TYC designations. However, some binary stars have two HD
# designations, for the primary and the secondary, but one of them is mapped to the HIP
# designation. In this case, the HIP index is applied to the primary, so that the HD
# index is applied to the correct star.
#
# Although the IAU states that the name of a star system only applies to the brightest
# component, many sources treat the names as if they were designations (e.g. "Sirius A"
# and "Sirius B"), and this file does the same.
END

my $starscount = 0;

# TODO: clean up this code
foreach my $index (sort keys %objects)
{
    my $namearray = $objects{$index}{'Names'};
    my $Barycenter = $objects{$index}{'Barycenter'};
    my $HIP = $objects{$index}{'HIP'};
    my $Parent = $objects{$index}{'Parent'};
    my $RA = $objects{$index}{'RA'};
    my $Dec = $objects{$index}{'Dec'};
    my $Distance = $objects{$index}{'Distance'};
    my $SpT = $objects{$index}{'SpectralType'};
    my $SpTNote = $objects{$index}{'SpectralTypeNote'};
    my $AppMag = $objects{$index}{'AppMag'};
    my $AbsMag = $objects{$index}{'AbsMag'};
    my $MagNote = $objects{$index}{'MagNote'};
    my $Radius = $objects{$index}{'Radius'};
    my $Period = $objects{$index}{'Period'};
    my $SemiMajorAxis = $objects{$index}{'SemiMajorAxis'};
    my $AxisNote = $objects{$index}{'AxisNote'};
    my $InclinationNote = $objects{$index}{'InclinationNote'};
    my $NodeNote = $objects{$index}{'NodeNote'};
    my $ArgNote = $objects{$index}{'ArgNote'};
    my $EpochNote = $objects{$index}{'EpochNote'};
    my $FullySpecified = $objects{$index}{'FullySpecified'};
    my $Eccentricity = $objects{$index}{'Eccentricity'};
    my $Inclination = $objects{$index}{'Inclination'};
    my $AscendingNode = $objects{$index}{'AscendingNode'};
    my $ArgOfPericenter = $objects{$index}{'ArgOfPericenter'};
    my $MeanAnomaly = $objects{$index}{'MeanAnomaly'};
    my $Epoch = $objects{$index}{'Epoch'};
    my $RotationPeriod = $objects{$index}{'RotationPeriod'};
    my $Duplicate = $objects{$index}{'Duplicate'};
    
    # Check if there's a trailing period, and remove it
    if ($SemiMajorAxis && $SemiMajorAxis =~ /\.$/) {
        $SemiMajorAxis = substr($SemiMajorAxis, 0, -1);
    }
    
    print OUTPUT "\n";
    # print OUTPUT "# $index\n";
    if ($Barycenter) {
        print OUTPUT "Barycenter ";
    } elsif ($Duplicate) {
        print OUTPUT "Replace ";
    }
    print OUTPUT "$HIP " if ($HIP);
    
    my $names = join(":", @$namearray);
    print OUTPUT "\"$names\"\n";
    print OUTPUT "{\n";
    
    if ($Parent && $Period) {
        print OUTPUT "\tOrbitBarycenter \"$Parent\"\n";
    } else {
        print OUTPUT sprintf("\tRA %.8f\n", $RA);
        print OUTPUT sprintf("\tDec %.8f%s\n", $Dec, $AxisNote);
        print OUTPUT sprintf("\tDistance %.3f\n", $Distance);
    }
    
    if (!$Barycenter) {
        print OUTPUT "\tSpectralType \"$SpT\"$SpTNote\n";
        if ($AppMag) {
            print OUTPUT "\tAppMag $AppMag$MagNote\n";
        }
        elsif ($AbsMag) {
            if (length($AbsMag - int($AbsMag)) > 3) {
                print OUTPUT sprintf("\tAbsMag %.2f%s\n", $AbsMag, $MagNote);
            } else {
                print OUTPUT sprintf("\tAbsMag %s%s\n", $AbsMag, $MagNote);
            }
        } else {
            print OUTPUT "\tAbsMag $MagNote\n";
        }
        print OUTPUT "\tRadius $Radius\n" if ($Radius);
        $starscount++;
    }
    
    if ($Parent && $Period) {
        print OUTPUT "\n";
        print OUTPUT "\tEllipticalOrbit {$FullySpecified\n";
        print OUTPUT "\t\tPeriod          $Period\n";
        print OUTPUT "\t\tSemiMajorAxis   $SemiMajorAxis$AxisNote\n";
        print OUTPUT "\t\tEccentricity    $Eccentricity\n" if ($Eccentricity);
        print OUTPUT "\t\tInclination     $Inclination$InclinationNote\n" if ($Inclination);
        print OUTPUT "\t\tAscendingNode   $AscendingNode$NodeNote\n" if ($AscendingNode);
        print OUTPUT "\t\tArgOfPericenter $ArgOfPericenter$ArgNote\n" if ($ArgOfPericenter);
        print OUTPUT "\t\tEpoch           $Epoch$EpochNote\n" if ($Epoch);
        print OUTPUT "\t\tMeanAnomaly     $MeanAnomaly\n" if ($MeanAnomaly);
        print OUTPUT "\t}\n";
    }
    
    # TODO: make this code less ugly
    my $RotationNote = "";
    # If eccentricity < 0.01 and the orbital period is less than a year,
    # assume bodies are tidally locked
    if ((!$Eccentricity || $Eccentricity < 0.01) && ($Period && $Period < 1) && !$RotationPeriod) {
        $RotationPeriod = $Period * 365.25 * 24;
        $RotationPeriod = FormatSigFigs($RotationPeriod, CountSigFigs($Period));
        $RotationNote = "# Assume tidal locking";
    }
    if ($RotationNote && !$Barycenter && $Inclination) {
        print OUTPUT "\n";
        print OUTPUT "\t$RotationNote\n";
        print OUTPUT "\tUniformRotation {\n";
        print OUTPUT "\t\tPeriod          $RotationPeriod\n" if ($RotationPeriod);
        print OUTPUT "\t\tInclination     $Inclination\n" if ($Inclination);
        print OUTPUT "\t\tAscendingNode   $AscendingNode\n" if ($AscendingNode);
        print OUTPUT "\t}\n";
    } elsif ($RotationPeriod) {
        $RotationNote = "# " . $RotationPeriod / 24 . " d" if (!$RotationNote);
        print OUTPUT "\n";
        print OUTPUT "\tRotationPeriod $RotationPeriod $RotationNote\n";
    }
    
    print OUTPUT "}\n";
}

print "$starscount stars in $systemscount systems newly added.\n" if $verbose;

# =========================== FUNDAMENTAL MATHEMATICAL FUNCTIONS ============================= #

# Log, base 10.
sub log10
{
    my $value = shift;
    return log($value)/log(10);
}

# Mod function, that allows negative and decimal inputs. Note that this is actually quite
# inefficient, but it has allowed me to catch errors such as accidentally typing "B" instead
# of "JD" as the units of time.
sub mod
{
    my ($a, $b) = @_;
    if ($a < 0) {
        while ($a < 0) {
            $a += $b;
        }
    } elsif ($a >= $b) {
        while ($a >= $b) {
            $a -= $b;
        }
    }
    return $a;
}

# ======================= FUNCTIONS FOR ESTIMATING STELLAR PARAMETERS ======================== #

sub Teff_to_SpT
{
    my $teff = shift;
    my $st = '';
    my $dist = 9999;
    foreach my $test (sort {$SpTeff{$a} <=> $SpTeff{$b}} keys %SpTeff)
    {
        # Check each spectral type to see which one has the closest temperature
        if (abs($SpTeff{$test} - $teff) < $dist)
        {
            $st = $test;
            $dist = abs($SpTeff{$test} - $teff);
        }
    }
    return $st;
}

sub Mass_to_SpT
{
    my $mass = shift;
    my $st = '';
    my $dist = 9999;
    foreach my $test (sort {$SpMass{$a} <=> $SpMass{$b}} keys %SpMass)
    {
        # Check each spectral type that is main-sequence
        # and see which one has the closest mass
        if ($test =~ /\dV/ && abs($SpMass{$test} - $mass) < $dist)
        {
            $st = $test;
            $dist = abs($SpMass{$test} - $mass);
        }
    }
    return $st;
}

sub AbsMag_to_SpT
{
    my $absmag = shift;
    my $st = '';
    my $dist = 9999;
    foreach my $test (sort {$SpAbsMag{$a} <=> $SpAbsMag{$b}} keys %SpAbsMag)
    {
        # Check each spectral type to see which one has the closest absolute magnitude
        if (abs($SpAbsMag{$test} - $absmag) < $dist)
        {
            $st = $test;
            $dist = abs($SpAbsMag{$test} - $absmag);
        }
    }
    return $st;
}

sub App_to_AbsMag
{
    my $appmag = $_[0];
    my $dist_ly = $_[1];
    my $absmag = $appmag - 5 * log10($dist_ly / (10 * $LY_TO_PARSEC));
    return $absmag;
}

# Calculate bolometric correction from temperature
# Equation from Reed (1998), JRASC Vol.92, p.36
# "The Composite Observational-Theoretical HR Diagram"
# https://ui.adsabs.harvard.edu/abs/1998JRASC..92...36R/abstract
sub Teff_to_BC
{
    my $Teff = shift;
    my $BC = -8.499 * (log10($Teff) - 4)**4
    + 13.421 * (log10($Teff) - 4)**3
    - 8.131 * (log10($Teff) - 4)**2
    - 3.901 * (log10($Teff) - 4)
    - 0.438;
    return $BC;
}

# Stefan-Boltzmann law: calculate luminosity from radius and temperature
sub SBL_to_Lum
{
    my $radius = $_[0];
    my $Teff = $_[1];
    my $Lum = (4 * pi * (($radius * $SOLAR_RADIUS * 1000)**2) * 5.67036713e-8 * ($Teff**4)) / 3.828e+26;
    return $Lum;
}

sub SpT_to_Mass
{
    my $SpT = shift;
    my $SpClass = "";
    my $SpSubClass = "";
    my $SpLum = "";
    my $SpTKey = "";
    if ($SpT =~ /([OBAFGKM])([0-9])+[\/.-]?[0-9]?([-\/IabV*]+)/) {
        ($SpClass, $SpSubClass, $SpLum) = (substr($1, 0, 1), substr($2, 0, 1), $3);
    } else {
        return "";
    }
    if ($SpLum =~ /^(V|IV|III|II|Ib|Iab|Ia)$/) {
        # These values don't need to be changed
    } elsif ($SpLum !~ /(I|V)/) {
        return "";
    } elsif ($SpLum =~ /^V/ or $SpLum eq "IV-V") {
        $SpLum = "V";
    } elsif ($SpLum =~ /V/) {
        $SpLum = "IV";
    } elsif ($SpLum =~ /III/) {
        $SpLum = "III";
    } elsif ($SpLum =~ /II/) {
        $SpLum = "II";
    } elsif ($SpLum =~ /Ib/) {
        $SpLum = "Ib";
    } elsif ($SpLum =~ /Iab/ or $SpLum eq "I") {
        $SpLum = "Iab";
    } elsif ($SpLum =~ /Ia/) {
        $SpLum = "Ia";
    }
    
    # Mamajek's table has all the spectral types from O3 to M9, but
    # Straizys' table does not.
    # First, "round" missing spectral types in Straizys' table,
    # to the nearest one
    if ($SpLum ne "V") {
        if ($SpClass eq "O" && $SpSubClass < 5) {
            $SpSubClass = 4;
        } elsif (($SpClass eq "B" || $SpClass eq "A") && $SpSubClass == 4) {
            $SpSubClass = 5;
        } elsif ($SpClass eq "A" && $SpSubClass > 5) {
            $SpSubClass = 7;
        } elsif ($SpClass eq "F" || $SpClass eq "G") {
            if ($SpSubClass == 1 || $SpSubClass == 3) {
                $SpSubClass = 2;
            } elsif ($SpSubClass == 4 || $SpSubClass == 6) {
                $SpSubClass = 5;
            } elsif ($SpSubClass > 6) {
                $SpSubClass = 8;
            }
        } elsif ($SpClass eq "K" && $SpSubClass > 5) {
            $SpSubClass = 5;
        } elsif ($SpClass eq "M" && $SpSubClass > 6) {
            $SpSubClass = 6;
        }
    }
    # Next, clean up subgiant stars.
    if ($SpLum eq "IV") {
        if ($SpClass eq "K" && $SpSubClass > 1) {
            $SpSubClass = 2;
        }
        # M-type subgiants should not exist. Therefore, this subroutine doesn't
        # attempt to "round" them off to the nearest value.
    # Clean up giant stars.
    } elsif ($SpLum eq "III") {
        if ($SpClass eq "F" && $SpSubClass == 8) {
            $SpSubClass = 5;
        } elsif ($SpClass eq "G" && $SpSubClass == 0) {
            $SpSubClass = 2;
        }
    } elsif ($SpLum =~ /^(II|Ib|Iab|Ia)$/) {
        if ($SpClass eq "O" && $SpSubClass == 5 && $SpLum eq "Ia") {
            $SpSubClass = 6;
        } elsif ($SpClass eq "K" && $SpSubClass == 4) {
            $SpSubClass = 5;
        } elsif ($SpClass eq "M" && $SpSubClass > 3) {
            $SpSubClass = 3;
        }
    }
    $SpTKey = $SpClass . $SpSubClass . $SpLum;
    return $SpMass{$SpTKey};
}

# Converts angular separation (in arcsecs) to physical separation (in AU), using distance (in ly).
sub Arcsec_to_AU
{
    my $arcsec = $_[0];
    my $dist_ly = $_[1];
    my $sep_ly = 2 * $dist_ly * tan(deg2rad($arcsec/3600) / 2);
    return $sep_ly * 63241.1;
}

sub OrbitTransform
{
    my $RA = $_[0];
    my $Dec = $_[1];
    my $period = $_[2]; # in years, to match Grant Hutchison's spreadsheet
    my $inclination = $_[3];
    my $node = $_[4];
    my $periarg = $_[5];
    
    # If missing, set default values
    $RA = 0 if (!$RA);
    $Dec = 0.00001 if (!$Dec); # Prevent divide-by-zero errors
    $inclination = 0 if (!$inclination);
    $node = 0 if (!$node);
    $periarg = 0 if (!$periarg);
    
    # Prevent divide-by-zero errors
    $inclination = 89.99999 if ($inclination == 90);
    
    my $alpha0 = deg2rad($RA) - pi;
    my $delta0 = -deg2rad($Dec);
    my $l270 = deg2rad($node);
    my $b = deg2rad(90-$inclination);
    my $alpha = atan(cos($b)*cos($l270)/(sin($b)*cos($delta0)-cos($b)*sin($delta0)*sin($l270)))+$alpha0;
    $alpha = $alpha + pi if (sin($b)*cos($delta0)-cos($b)*sin($delta0)*sin($l270) < 0);
    my $delta = asin(cos($b)*cos($delta0)*sin($l270)+sin($b)*sin($delta0));
    my $epsilon = deg2rad($OBLIQUITY);
    my $lambda = atan((sin($alpha)*cos($epsilon)+tan($delta)*sin($epsilon))/cos($alpha));
    $lambda = $lambda + pi if (cos($alpha) < 0);
    my $beta = asin(sin($delta)*cos($epsilon)-cos($delta)*sin($epsilon)*sin($alpha));
    my $om270 = ($node-270)/180*pi;
    my $deltaom = asin(cos($delta0)*sin($om270));
    my $alphaom = atan(cos($om270)/-sin($delta0)/sin($om270))+$alpha0;
    $alphaom = $alphaom + pi if (-sin($delta0)*sin($om270) < 0);
    my $lambdaom = atan((sin($alphaom)*cos($epsilon)+tan($deltaom)*sin($epsilon))/cos($alphaom));
    $lambdaom = $lambdaom + pi if (cos($alphaom) < 0);
    my $betaom = asin(sin($deltaom)*cos($epsilon)-cos($deltaom)*sin($epsilon)*sin($alphaom));
    my $d = acos(cos($betaom)*cos($lambdaom-$lambda - pi/2));
    $d = -$d if ($betaom < 0);
    my $finalinclination = 90-$beta/pi*180;
    my $finalnode = mod($lambda/pi*180+90, 360);
    my $finalperi = mod($periarg+$d/pi*180, 360);
    return $finalinclination, $finalnode, $finalperi;
}

# ================================ STRING-HANDLING FUNCTIONS ================================= #

sub CountDecimalDigits
{
    my $num = shift;
    my $numdigits = 0;
    if ($num =~ /\d+\.(\d+)/) {
        $numdigits = length($1);
    }
    return $numdigits;
}

# Figure out what the "children" of a specific component designation would be.
# Note that three-letter components are typically ambiguous (i.e. not possible to tell if
# a system is AB-C or A-BC), so they have to be fixed by the next row in the catalog.
sub GetComponents
{
    my ($designation, $mult) = @_;
    my ($a, $b) = ("", "");
    if (length($designation) == 1) {
        $a = $designation . "a";
        $b = $designation . "b";
    } elsif (length($designation) == 2) {
        if ($designation =~ /[A-Z][A-Z]/) {
            $a = substr($designation, 0, 1);
            $b = substr($designation, -1);
        } elsif ($designation =~ /[A-Z][a-z]/) {
            $a = $designation . "1";
            $b = $designation . "2";
        }
    } elsif ($designation =~ /^(A[C-Z])B/) {
        $a = $1;
        $b = "B";
    } elsif ($designation =~ /^A([C-Z][C-Z])/) {
        $a = "A";
        $b = $1;
    # Use $mult to help figure out components
    } elsif ($mult =~ /(.+)-(.+)/) {
        $a = $1;
        $b = $2;
    }
    return $a, $b;
}

sub StripComponents
{
    my @array = @_;
    my @newarray = ();
    for (my $i = 0; $i < scalar @array; $i++) {
        if ($array[$i] =~ /^(ADS |CCDM J)(\d+|\d+[-\+]\d+)\s*([A-Za-z]+)$/) {
            push @newarray, $1 . $2;
        } else {
            push @newarray, $array[$i];
        }
    }
    return @newarray;
}

sub GetSibling
{
    my ($parent, $child) = @_;
    my $sibling = $parent;
    if ($parent =~ /$child/) {
        $sibling =~ s/$child//;
        $sibling =~ s/\-//;
    }
    return $sibling;
}

sub TrimSpType
{
    my $sptype = shift;
    if ($sptype =~ /([OBAFGKM]\d)[-\/][OBAFGKM]*\d([IV]+)/) {
        $sptype = $1 . $2;
    }
    return $sptype;
}
