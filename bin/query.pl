# ===================================== #
#                query.pl               #
# ===================================== #
# This perl script take the strandless circRNA coordinates and extract the left and right boundaries for sequence queries.
# @ARGV[0]: the index file
# @ARGV[1]: the genome file
# @ARGV[2]: the coordinates file (tab-delimited with 5 columns, chr, strand, start, end, name)
# @ARGV[3]: the output sequence file
# @ARGV[4]: specify the format of output (0: raw, 1: fasta, 2: fastq)
# test：perl bin/query.pl /Users/minzhe/Documents/Project/Clirc/data/test/library/index.txt /Users/minzhe/Documents/Project/Clirc/data/test/chr4.fa.masked /Users/minzhe/Documents/Project/Clirc/data/test/library/circRNA_coordinates.txt /Users/minzhe/Documents/Project/Clirc/data/test/library/halves.txt 1
#!/usr/bin/perl

use strict;
use warnings;
use POSIX 'floor';

my ($index, $genome, $coords, $seq, $format) = @ARGV;
my ($name, %index, @items, $chr, $strand, $start, $end, $s, $str, $line_s, $line_e, $sequence, $region);
my $jump = 5000;

##########    read index file    ##########
open(FILE_IN, $index) or die("Open file $index failed: $!");

while (<FILE_IN>) {
    $_ =~ s/\n//;
    @items = split("\t", $_);
    $index{$items[0]} = [$items[1], $items[2]];
}

close(FILE_IN);

##########    query sequence    ##########
open(FILE_IN1, $coords) or die("Open file $coords failed: $!");
open(FILE_IN2, $genome) or die("Open file $genome failed: $!");
open(FILE_OUT, "> " . $seq);

while ($region = <FILE_IN1>) {

    # locate start and end
    ($chr, $strand, $start, $end, $name) = split("\t", $region);
    $s = floor($start/$jump) * $jump + 1;
    if ($s == $start+1) { $s -= $jump; }    # boundary case
    unless (exists $index{$chr . "_" . $s}) {
        print "The following region dose not exist: ". $region;
    }
    $sequence = "";

    ### the first segment
    @items = @{$index{$chr . "_" . $s}};
    seek FILE_IN2, $items[1], 0;            # go to the line
    $line_s = $s - $items[0] + 1;           # the coordinate of the first base on that line

    while ($line_s <= $start) {             # untile the line in which $strat resides is read
        $str = <FILE_IN2>;
        unless (defined $str && $str !~ />/) {
            print "The following region dose not exist: " . $region;
        }
        $line_s += (length($str) - 1);      # the coordinate of the first base on each line
    }

    $str =~ s/\n//;
    $line_e = $line_s - 1;
    $line_s -= length($str);

    if ($end > $line_e) {                   # $start and $end are not on the same line
        $sequence = substr $str, ($start - $line_s);
    } else {                                # $start and $end are on the same line
        $sequence = substr $str, ($start - $line_s), ($end - $start + 1);
        print_seq($sequence, $strand, $end-$start+1, $region, $format, $name);
        next;
    }


    ### the middle segments
    $line_s += length($str);
    while ($line_s <= $end) {
        $str = <FILE_IN2>;
        unless (defined $str && $str !~ />/) {
            print "The following region dose not exist: " . $region;
        }
        $str =~ s/\n//;
        $line_s += length($str);
        $sequence .= $str
    }


    ### the last segment
    if ($line_s-$end-1 != 0) {
        $sequence = substr $sequence, 0, -($line_s - $end - 1);
    }
    print_seq($sequence, $strand, $end-$start+1, $region, $format, $name);
}

close(FILE_IN1);
close(FILE_IN2);
close(FILE_OUT);

print "Done extracting circRNA sequences.\n";
exit;

##########    subroutine    ##########
sub print_seq {
    my ($sequence, $strand, $len, $region, $format, $name) = @_;

    if ($len != length($sequence)) {
        print "The following region dose not exist: ". $region;
    }

    if ($strand eq "-") {
        $sequence =~ tr/ABCDGHMNRSTUVWXYabcdghmnrstuvwxy/TVGHCDKNYSAABWXRtvghcdknysaabwxr/;
        $sequence = reverse($sequence);
    }

    $sequence = uc $sequence;

    if ($format >= 1) { print FILE_OUT ">".$name; }
    print FILE_OUT $sequence . "\n";
    if ($format == 2) { print FILE_OUT "+\n-\n"; }
}