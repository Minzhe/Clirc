# ============================================= #
#              filter_circRNA.pl                #
# ============================================= #
# this perl script combines and filters output from several files generated by find_circRNA.pl
# min: minimum number of overlapping reads to call a true circRNA binding site
# ratio: minimum number of unique start/end positions divided by all start/end positions
# len: minimum length of the CLIP cluster on each circRNA
# output_file: output file
# library: library folder where gsnap index files are stored
# circ_files: input files separated by space
#!/usr/bin/perl
use strict;
use warnings;
use List::MoreUtils qw(uniq);
use List::Util qw(max);

my ($min, $ratio, $len, $output_file, $library, @circ_files) = @ARGV;
my (@items, $circ, %circ, $circ_file, $line, $header);
my ($n_unique, @starts, @ends, $total);
my (%circ_seq);

######  read each file generated by find_circRNA.pl  ###############

foreach $circ_file (@circ_files) {
	open(FILE_IN, $circ_file) or die "Cannot open file ".$circ_file."\n";
	$header = <FILE_IN>;
	$header =~ /Total mapped\: ([0-9]+)/;
	$total += $1;
	$header = <FILE_IN>;

	while ($line = <FILE_IN>) {
		@items = split("\t", $line);
		$circ = join("\t", @items[(0, 1, 2, 4)]); # spell the circRNA
		unless (exists $circ{$circ}) {$circ{$circ} = {};}
		$circ{$circ}->{$items[5]."_".$items[6]} = $line;
	}

	close(FILE_IN);
}

###############  write results  ##########################

open(FILE_OUT, "> " . $output_file) or die "Cannot write to file " . $output_file . "\n";
print FILE_OUT "#Total mapped: " . $total . "\n";
print FILE_OUT $header;

# read circRNA sequences
open(FILE_SEQ, $library."/circRNA.fa");
while ($line = <FILE_SEQ>) {
	$line =~ /\>(.*)_([0-9]+)_([0-9]+)_circ_[0-9]+_[0-9]+/;
	$circ = $1 . "\t" . $2 . "\t" . $3;
	if (exists $circ_seq{$circ}) { print $circ . " already exists!\n"; }
	$circ_seq{$circ} = <FILE_SEQ>;
}
close(FILE_SEQ);

foreach $circ (keys %circ) {
	# there must be enough number of reads
	$n_unique = scalar(values %{$circ{$circ}});
	if ($n_unique<$min) { next; }
	
	# there must be enough unique start/end positions
	@starts = ();
	@ends = ();
	foreach (keys %{$circ{$circ}}) {
		$_ =~ /^([0-9]+?)_([0-9]+?)$/;
		push @starts, $1;
		push @ends, $2;
	}
	if (scalar(uniq @starts) <= $ratio*$n_unique) { next; } 
	if (scalar(uniq @ends) <= $ratio*$n_unique) { next; }

	# the length of the CLIP cluster is too small
	if (max(@starts)+max(@ends) <= $len) { next; }

	# print 
	print FILE_OUT values %{$circ{$circ}}; # result file
}

close(FILE_OUT);

print "Done removing adaptors.\n";
exit;

##########    subroutine    ##########

sub extract_sequence {
	my ($circ, $circ_seq_ref, $max_start, $max_end) = @_;
	my @items = split("\t", $circ);
	my $seq = uc $$circ_seq_ref{join("\t", @items[0 .. 2])}; # get sequence
	my $rc = "";

	$seq = ~s/\n//;
	$seq = substr $seq, length($seq)/2-$max_start, $max_start+$max_end; # regions that are bound by RBP at least once
	if ($items[3] eq "+") { return $seq."\n"; } # positive strand
	$seq = scalar reverse $seq; # reverse complement

	while ($seq =~ /([ATCGN])/g) {
		if ($1 eq "A") { $rc .= "T"; }
		if ($1 eq "T") { $rc .= "A"; }
		if ($1 eq "C") { $rc .= "G"; }
		if ($1 eq "G") { $rc .= "C"; }
		if ($1 eq "N") { $rc .= "N"; }
	}

	return $rc . "\n";
}


