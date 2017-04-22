# ===================================== #
#           circRNA2region.pl           #
# ===================================== #
# This perl script take the strandless circRNA coordinates and extract the left and right boundaries for sequence queries.
# $file: the length of the RNA sequences to expand from the junction site
# @items: the absolute path of output file (a teb-delimited file with 5 columns, chromosomes, +, start coordinates, end coordinates and id)
# $id: the absolute path of input file (coordinates of known circRNAs, a tab-delimited file with 3 columns, chromosomes, start coordinates and end coordinates)
# test: perl bin/circRNA2region.pl 50 /Users/minzhe/Documents/Project/Clirc/data/test/library/circRNA_coordinates.txt /Users/minzhe/Documents/Project/Clirc/data/test/dm3_chr4_circRNAs.txt
#!/usr/bin/perl

use List::Util qw(min max);
use strict;
use warnings;

my ($file, @items, $id);
my ($expand, $output_file, @input_files) = @ARGV;

open(FILE_OUT, "> " . $output_file);

foreach $file (@input_files) {
    open(FILE_IN, $file) or die("Open file $file failed: $!");

    while (<FILE_IN>) {
        @items = split("\t", $_);
        $items[2] =~ /([0-9]+)/;
        $items[2] = $1;
        $id = $items[0] . "_" . $items[1] . "_" . $items[2] . ":" . $file;
        print FILE_OUT $items[0] . "\t+\t" . $items[1] . "\t" . min($items[2], $items[1]+$expand-1) . "\t" . $id . "\n";
        print FILE_OUT $items[0] . "\t+\t" . max($items[1], $items[2]-$expand+1) . "\t" . $items[2] . "\t" . $id . "\n";
    }

    close(FILE_IN);
}

close(FILE_OUT);

print "Done extracting circRNA coordinates.\n";
exit;