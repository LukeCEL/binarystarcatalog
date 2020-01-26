# binarystarcatalog
Creates database of binary and multiple stars for Celestia.

## About
This catalogue was intended to increase the number of stars in Celestia, an open-source space simulator. In the default version, only a few hundred systems are rendered as binary, even though a significant fraction of stars are in star systems. Many famous star systems, such as Algol, appear as single stars. This catalogue was intended to increase this fraction.

## Prerequisites
This script was written in Perl. In addition to Perl, you also need these modules:

* Math::SigFigs
* Data::Dumper
* Array::Utils

It is recommended that you install both of these modules using `cpanm`. You can find instructions for installing `cpanm` here at https://www.cpan.org/modules/INSTALL.html.

## Usage
1. Install Perl and the prerequisite modules.

2. Download the Binary Star Catalogue (the latest version can be found [here](https://celestia.space/forum/viewtopic.php?f=23&t=20498).

3. In the command-line, `cd` (change the current directory) to the one that has the source catalog (catalog.txt).

4. Move the .stc file visualbins.stc (which can be found in the default Celestia package) to the same directory as catalog.txt.

5. In the command-line, type: `perl readcatalog.pl`. Instead of `readcatalog.pl`, you may also add the location of the file relative to the directory from which you are running.

The Perl script should automatically create an SSC file called `binarycatalog.stc`.

This script has one command-line argument, `-v` or `--verbose`. If specified, the script will output its current status, as well a count of the number of systems and stars added.

## Function
This Perl script reads data from a tab-separated text file, called catalog.txt, which contains all the data that is needed to build the .stc file. The two are designed to work with each other, and there are several idiosyncrasies that are specific to the way the handles systems.

First, the script reads each line of catalog.txt and puts it into a hash called `%binaries`. Then, the script goes through each entry into `%binaries`, and adds all the data for each _object_ (either a barycenter or star) into `%objects`. Finally, the script outputs the .stc file from `%objects`.

### Nomenclature
For each star, names are added in the following order: 1) the `ProperNames` column, 2) the `Names1` or `Names2` column if not HD, SAO, or TYC, 3) the `BAY` column, 4) the `FLA` column, 5) the `ADS` column, 6) the `CCDM` column, and 7) the remaining names of the `Names1` or `Names2` column.

The `mult` column has the component designations for the pair. For example, if the `mult` is "A", then the primary's parameters refer to component Aa and the secondary's, to component Ab. The `ADS` and `CCDM` have their own component designations, because they sometimes differ from the typical component designations.

If a name appears in `ProperNames` or `OtherNames`, but does not in the first line of the system, then the component designations for that name start at A and B. Therefore, Castor Ca = YY Gem A, and Castor Cb = YY Gem B.

Some systems are already in default Celestia, so they need to be replaced to avoid duplicates. For this reason, the script reads visualbins.stc to find which stars are already in Celestia.

### Hierarchies
Multiple star systems are generally hierarchical, with stars orbiting in pairs, and those pairs orbiting in more pairs, and so on. This hierarchical structure is reflected in how this script builds systems.

Suppose this script reads a line where the `mult` is AB. The script creates three objects, a barycenter entry with a component designation of AB and two entry whose component designations are A and B, respectively. Barycenters should be marked by an asterisk (*) in the spectral type columns (`SpT A` and `SpT B`). If the next line is A, then instead of creating a new barycenter entry with a component designation of A, it simply adds information onto the existing A entry.

In Celestia, barycenters must be declared before the stars orbiting said barycenter, but besides that there are no restrictions on the order of stars in a catalog file. That said, this script orders star systems, first alphabetically by the Latin transcription of the Greek letter of a Bayer designation, then by the Flamsteed designation if missing, then by the (Latin-letter) Bayer designation. Within a system, stars and objects follow an 
"ABACABA pattern": for example AB, A, Aa, Ab, B, Ba, and Bb.

### Parameters
Knowledge on binary and multiple stars is very incomplete, and the parameters necessary for Celestia (e.g. magnitudes, spectral types) often have to be estimated from known parameters. Therefore, the following sources are employed to fill in the missing data:

Mamajek (2019) "[A Modern Mean Dwarf Stellar Color and Effective Temperature Sequence](http://www.pas.rochester.edu/~emamajek/EEM_dwarf_UBVIJHK_colors_Teff.txt)".

Reed (1998), JRASC Vol.92, p.36
"[The Composite Observational-Theoretical HR Diagram](https://ui.adsabs.harvard.edu/abs/1998JRASC..92...36R/abstract)".

Straizys and Kuriliene (1981), Ap&SS 80, p.353
"[Fundamental stellar parameters derived from the evolutionary tracks](https://ui.adsabs.harvard.edu/abs/1981Ap%26SS..80..353S/abstract)".

### Orbits
Orbits in Celestia must be in the default J2000 ecliptic reference frame. Since orbits in the scientific literature are given relative to the plane-of-sky, the orbits must be transformed for Celestia. This script uses the equations from Grant Hutchison's spreadsheet:

Grant Hutchison's Star Orbit Translation: starorbs.xls
https://www.classe.cornell.edu/~seb/celestia/hutchison/spreadsheets.html#2

Additionally, orbital elements are often given in slightly different formats that must be taken into account. The main varieties are 1) a photocentric orbit, which is the orbit of the photocenter around the center of motion (this can usually be interpreted as the orbit of the primary relative to the barycenter); 2) for spectroscopic orbits, the argument of pericenter ω can refer to the primary, instead of the secondary; 3) for eclipsing binaries, T<sub>min</sub> (the epoch of the primary minimum) is sometimes given instead of the periastron epoch.

Ambiguity over whether a left-handed or right-handed coordinate system is used means that while a binary system's apparent orbit in the sky will be correct, it may be inconsistent with radial velocity observations (i.e. receding or approaching the observer at the wrong times). To fix this, this manually has to be checked in Celestia, and if the conventions are reversed relative to what Hutchison's spreadsheet uses, then the `flip` column says yes, otherwise it says no.

## Acknowledgements
This add-on has made extensive use of the SIMBAD database, operated at CDS, Strasbourg, France.

Thank you to Chris Laurel and everyone who helped make Celestia great in the first place. Additionally, thank you to Grant Hutchison for the orbital transformation, without which this add-on would not be possible.