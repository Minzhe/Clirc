# ===================================== #
#            half2library.pl            #
# ===================================== #
# This perl script is to extract circRNA half sequence and assembly into a fasta file for alignment index building
# $input_fil: input file
# $output_file: output file
# demo: perl bin/half2library.pl /Users/minzhe/Documents/Project/Clirc/data/test/library/halves.txt /Users/minzhe/Documents/Project/Clirc/data/test/library/circRNA.fa
#!/usr/bin/perl

use strict;
use warnings;

my ($input_file, $output_file) = @ARGV;
my ($N, $header, $left, $right, %circ, $coord, $source);
my $N_limit = 0.4;

open(FILE_IN, $input_file) or die("Open file $input_file failed: $!");
open(FILE_OUT, "> " . $output_file);

while ($header = <FILE_IN>) {
    $left = <FILE_IN>;
    $right = <FILE_IN>;
    $right = <FILE_IN>;
    $header =~ s/\n//;
    ($coord, $source) = split(":", $header);
    $left =~ s/\n//;
    $right =~ s/\n//;
    if (! exists $circ{$right . $left}) {
        $circ{$right . $left} = $coord . "_circ_" . length($right) . "_" . length($left);
    }
}

foreach (keys %circ) {
    $N = () = $_ =~ /N/gi;        # too many Ns
    if ($N > $N_limit * length($_)) { next; }
    print FILE_OUT $circ{$_} . "\n" . $_ . "\n";
}

close(FILE_IN);
close(FILE_OUT);

print "Done assembling circRNA fastq file.\n";
exit;