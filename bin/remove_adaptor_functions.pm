# ===================================================== #
#             remove_adaptor_functions.pl               #
# ===================================================== #
# This perl module provide subroutines for to trim adaptors.

package remove_adaptor_functions;
use strict;
use warnings;
use List::Util qw(max min);
use Exporter 'import';

our @EXPORT = qw(reverse_read prepare_adaptor trim_read);

##########       subroutine       ##########
### reverse reads if the adaptor is at 5' end
sub reverse_read {
    my ($read_ref, $format) = @_;
    $$read_ref[1] =~ s/\n//;
    $$read_ref[1] = (scalar reverse $$read_ref[1]) . "\n";
    if ($format eq "fastq") {
        $$read_ref[3] =~ s/\n//;
        $$read_ref[3] = (scalar reverse $$read_ref[3]) . "\n";
    }
}

### prepare adaptor hash from the adaptor sequence
sub prepare_adaptor {
    my ($adaptor, $regex_ref) = @_;
    my $ad_len = 0;
    my $adaptor_len_max = 20;

    $adaptor = substr $adaptor, 0, $adaptor_len_max;
    $adaptor = uc $adaptor;
    $adaptor =~ s/ //g;

    while ($adaptor =~ /([ATCGN])/g) {
        $regex_ref -> {$ad_len} = {N => 1, $1 => 1, n => 1, lc $1 => 1};
        if ($1 eq "N") {
            $regex_ref -> {$ad_len} = {A => 1, T => 1, C => 1, G => 1, a => 1, t => 1, c => 1, g => 1};
        }
        $ad_len++;
    }
}

### trim adaptor from read
# AGCTAGCTAGCTAGCTAGCT\n
# 123456789...
#              i----- (1-based)
#              -----j (0-based)
sub trim_read {
    my ($read_ref, $regex_ref, $format, $mismatch) = @_;
    my ($i, $j, $len2compare, $count, $mismatch_count);
    my $sequence = $$read_ref[1];
    my @char = split("", $sequence);

    foreach $i (1 .. length($sequence)-1) {                                     # each position in the sequence read
        $len2compare = min(length($sequence)-$i, scalar(keys %$regex_ref));
        $mismatch_count = $len2compare * $mismatch;                             # limit of mismatch
        $count = 0;

        foreach $j (0 .. $len2compare-1) {
            if (! exists $regex_ref -> {$j} -> {$char[$i-1+$j]}) { $count++; }      # count mismatch
            if ($count > $mismatch_count) { last; }                                 # too many mismatches
        }

        if ($count > $mismatch_count) { next; }
        # a match is found
        $$read_ref[1] = substr $$read_ref[1], 0, $i-1;
        $$read_ref[1] .= "\n";
        if ($format eq "fastq") {
            $$read_ref[3] = substr $$read_ref[3], 0, $i-1;
            $$read_ref[3] .= "\n";
        }
    last;
    }
}

1;