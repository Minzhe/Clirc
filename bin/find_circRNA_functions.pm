# ===================================================== #
#             remove_adaptor_functions.pl               #
# ===================================================== #
# This perl module provide subroutines for to search circRNA.

package find_circRNA_functions;
use strict;
use warnings;
use List::Util qw(max min);
use Exporter qw(import);

our @EXPORT = qw(find_circ read_alignment);

##########       subroutine       ##########
### find evidence of circRNAs from a gread group
sub find_circ {
	my ($reads_ref, $mismatch, $overhang, $string) = @_;
	my ($circ_chr, $circ_strand, $circ_stat);
	my $linear_qual = -1;
	my $linear_len = -1;
	my $circ_qual = -1;
	my $circ_len = -1;
	my $circ_hash = -1;

	foreach (keys %$reads_ref) {
		is_circ($reads_ref->{$_}, $mismatch, $overhang);                # decide the type of this read
		
		if ($reads_ref->{$_}->{"type"} eq "circRNA") {                  # mapping to circRNA and indeed across the splice junction
   			if ($circ_qual>$reads_ref->{$_}->{"qual"}) { next; }        # unless it reaches a better mapping quality, skip
   			if ($circ_qual == $reads_ref->{$_}->{"qual"} && $circ_len > $reads_ref->{$_}->{"len"}) { next; }    # always keep the longer alignment
   			if ($circ_qual == $reads_ref->{$_}->{"qual"} && $circ_len == $reads_ref->{$_}->{"len"} && $circ_hash > $reads_ref->{$_}->{"hash"}) { next; } 

			$circ_qual = $reads_ref->{$_}->{"qual"};
   			$circ_len = $reads_ref->{$_}->{"len"};
			$circ_chr = $reads_ref->{$_}->{"chr"};
			$circ_strand = $reads_ref->{$_}->{"strand"};
			$circ_stat = $reads_ref->{$_}->{"stat"};
			$circ_hash = $reads_ref->{$_}->{"hash"};
	 	} elsif ($reads_ref->{$_}->{"type"} eq "linear") {              # mapping to linear form
			# record the best mapping quality and longest length to the linear form
			$linear_qual = max($reads_ref->{$_}->{"qual"}, $linear_qual);
			$linear_len = max($reads_ref->{$_}->{"len"}, $linear_len);
	 	}
	}

	if ($circ_qual <= $linear_qual || $circ_len <= $linear_len) { return "Not a circRNA\n"; }
	if ($string == 1 && $linear_qual > -1) { return "Not a circRNA\n"; } # more stringent, if any linear mapping is found, it's not a circRNA mapping
	return join("\t", ($circ_chr, $$reads_ref{0}->{"name"}, $circ_strand, $circ_stat))."\n";
}

### read alignment from each line of SAM file
# $line: each SAM file line
# $save: 0 or 1, whether to save results
# $ref: hash reference
# return name of the read
sub read_alignment {
	my ($line, $save, $hash_ref) = @_;
	my $i = scalar(keys %$hash_ref);
	my @items_read = split("\t", $line);
	if ($save == 0) { return $items_read[0]; }

	if ($items_read[1] == 256) {$items_read[1] = 0;}
	if ($items_read[1] == 272) {$items_read[1] = 16;}
	if ($items_read[1] != 0 && $items_read[1] != 16) { return $items_read[0]; } # unmapped
	if ($items_read[1] == 0) {$items_read[1] = "+";} else { $items_read[1] = "-"; } # strand to +/-

	$hash_ref->{$i} = {
        "name" => $items_read[0],
        "strand" => $items_read[1],
        "chr" => $items_read[2],
        "seq" => $items_read[9],
        "start" => $items_read[3],
        "qual" => $items_read[4],
        "CIGAR" => $items_read[5]
    };
	if ($line =~ /MD\:Z\:(.*?)(\t|$)/) { $hash_ref->{$i}->{"MD"} = $1; } else { $hash_ref->{$i}->{"MD"} = "na"; }
	return $items_read[0];
}

### whether one alignment corresponds to a circRNA
sub is_circ {
	my ($each_ref, $mismatch, $overhang) = @_;
	my ($circ_start, $m_start, $m_end, $md, $md_element);
	my ($start, $end, $cigar, $overhang_start, $overhang_end);
	
	# initialize
	$m_start = 0;
	$m_end = 0;
	$each_ref->{"type"} = "linear";

	# length of alignment
	$start = $each_ref->{"start"};
	$end = $start-1;
	$cigar = $each_ref->{"CIGAR"};
	while ($cigar =~ /([0-9]+)([MD])/g) { $end += $1; }
	$each_ref->{"len"} = $end - $start + 1;

	# if mapping is within a normal chromosome
	if ($each_ref->{"chr"}!~/circ_([0-9]+)_/) { return 1; } # mapped to linear form
	$circ_start = $1; # split circRNA annotation
	if ($cigar =~ /N/) { return 2; }	
	$overhang_start = $circ_start-$start+1;
	$overhang_end = $end-$circ_start;
	if (min($overhang_start, $overhang_end) < $overhang) { return 3; }; # not long enough or doesn't extend across junctions

	# Mismatches from insertions
	$end = $start-1;
	while ($cigar =~ /([0-9]+)([MDI])/g) {
		if ($2 eq "M" || $2 eq "D") { $end+=$1; }
		if ($2 eq "I") {
            if ($end <= $circ_start) { $m_start += $1; } else { $m_end+=$1; }
        }
	}

	# mismatches from deletions and substitutions
	$md = $each_ref->{"MD"};
	if ($md ne "na") {
		$end = $start-1;
		$md = uc $md;
		while ($md =~ /([ACTGN]|[0-9]+)/g) {
			$md_element = $1;

			if ($md_element =~ /[ACTGN]/) {
				$end++;
				if ($end <= $circ_start) { $m_start++; } else { $m_end++; }
			} elsif ($md_element =~ /[0-9]+/) {
				$end += $md_element;
			}
		}
	}

	$each_ref->{"type"} = "circRNA";
	$each_ref->{"hash"} = hash_string($each_ref->{"chr"});
	$each_ref->{"stat"} = join("\t", ($overhang_start, $overhang_end, $m_start, $m_end, $cigar, $md, $each_ref->{"seq"}));
	if ($m_start >= $overhang_start * $mismatch || $m_end >= $overhang_end * $mismatch) { $each_ref->{"type"} = "invalid_circRNA"; } # invalid mapping to circRNA
	return 4;
}

### hash a string by suming the ASCII code of each character
sub hash_string {
	my $string = $_[0];
	my $hash = 0;
	while ($string =~ /(.)/g) {$hash += ord($1);}
	return $hash;
}

1;

