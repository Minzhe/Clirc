# ===================================== #
#           Clirc_library.pl            #
# ===================================== #
# This perl script takes user input ro build a pseudo reference genome containing the regular reference genome together with user-supplied circRNA.
# arguments: see below help message.
# demo: perl bin/Clirc_library.pl -coord data/test/dm3_chr4_circRNAs.txt -genome data/test/chr4.fa.masked -library data/test/library
#!/usr/bin/perl

use strict;
use warnings;
use Cwd 'abs_path';
use Getopt::Long;
use File::Spec;


my ($coord, $genome, $library, $expand, $help, $message, $path);
$expand = 50;
$help = 0;
$path = abs_path($0);
$path =~ s/Clirc_library\.pl$//;

GetOptions(
    "coord=s" => \$coord,
    "genome=s" => \$genome,
    "expand=i" => \$expand,
    "library=s" => \$library,
    "help" => \$help
) or die("Error in command line arguments!\n");

$coord = File::Spec->rel2abs($coord);
$genome = File::Spec->rel2abs($genome);
$library = File::Spec->rel2abs($library);
if (!-d $library) {
    mkdir $library
}

if ($help == 1) {
    $message = << 'HELP';

    Clirc_library: build alignment library for identifying RBP-bound circRNAs from CLIP-Seq data
Options:
    -coord          coordinates of known circRNAs, a tab-delimited file with 3 columns, chromosomes, start coordinates and end coordinates
    -genome         file name (including path) of the reference genome fasta file
    -expand         if the circRNA is larger than this length, only expand nt of RNA sequences up- and down-stream of the junction site will be used to build library (default=50)
    -library        path to the directory where intermediate and gsnap index files will be stored. Specify an empty folder for use only by the Clirc software
    -help           print this help message

HELP
    print $message;
    exit;
}

0 == system("perl " . $path . "/circRNA2region.pl " . $expand . " " . $library . "/circRNA_coordinates.txt " . "$coord")
    or die("Extracting circRNA coordinates failed!\n");
0 == system("perl " . $path . "/index.pl " . $genome . " " . $library . "/index.txt")
    or die("Indexing genome failed!\n");
0 == system("perl " . $path . "query.pl " . $library . "/index.txt " . $genome . " " . $library . "/circRNA_coordinates.txt " . $library . "/halves.txt 1")
    or die("Extracting circRNA sequences failed!\n");
0 == system("perl " . $path . "half2library.pl " . $library . "/halves.txt " . $library . "/circRNA.fa")
    or die("Assembling circRNA fastq file failed!\n");
print "Calling gmap_build, it could be slow ...\n";
0 == system("gmap_build -D " . $library . " -d circRNA " . $genome . " " . $library . "/circRNA.fa >/dev/null 2>&1")
    or die("Calling gmap_build failed!\n");



