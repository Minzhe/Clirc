# =========================================== #
#              find_circRNA.pl                #
# =========================================== #
# This perl script finds evidence of circRNAs from alignment data
# (1) no duplicate collapsing will be done within this script;
# (2) assume that each sequencing read has a unique name;
# $sam_file: SAM format file
# $circ_file: file to write found circRNA to
# $overhang: cutoff for overhang
# $mismatch: max allowed mismatch rate
# $process_n: number of processes to spawn
# $string: 0 or 1, 0=low stringency, 1=high stringency
#!/usr/bin/perl

use strict;
use warnings;
use Parallel::ForkManager;
use IO::Pipe;
use IO::Select;
use List::Util qw(max min);
use Cwd 'abs_path';
my $path;
BEGIN {
	$path = abs_path($0);                 
	$path =~ s/find_circRNA\.pl//;              
}               
use lib $path;              
use find_circRNA_functions;               

# $SIG{__WARN__} = sub { die "Alignment file format error!\n" };              
$SIG{CHLD} = 'IGNORE';

my ($sam_file, $circ_file, $overhang, $mismatch, $process_n, $string) = @ARGV;
my ($result, $line, @items, $count, %reads, %circ, $name);
my ($save, $circRNA, $circRNA1);
my ($select, @select, $process, $pid, %pipes, $pm, @ready);
my $n_mapped = 0;

#############  prepare for forking  ##############

$pm = Parallel::ForkManager->new($process_n);                   # prepare for fork
$select = IO::Select->new();                                    # communication back to parent
foreach (1 .. $process_n) { $pipes{$_} = IO::Pipe->new(); } 

# -------------  child process starts  ---------------- #

foreach $process (1..$process_n) {      
	# establish communication with parent process
	$pid = $pm->start and next;                                 # do the fork
	$pipes{$process}->writer();
	$pipes{$process}->autoflush(1);
	
	open(FILE_IN, $sam_file) or die "Cannot open SAM file!\n";
	while ($line = <FILE_IN>) { if ($line !~ /^@/) { last; } }  # skip headers
	$count = 0;

	while ($line) {
		# read the process-th read group
        $count++;
		if ($count==$process_n) { $count=0; }
		if (($count-$process) % $process_n != 0) { $save = 0; } else { $save = 1; }

		%reads = ();
		$name = read_alignment($line, $save, \%reads);          # the first of each read group

		while ($line = <FILE_IN>) {                             # read the other ones and the first of the next read group
			if ($name ne read_alignment($line, $save, \%reads)) {
				map {
                    if ($name ne $reads{$_}->{"name"}) { delete $reads{$_}; }
                } keys %reads;
				last;
			}
		}


		# find evidence of circRNA
		if (scalar(keys %reads) == 0) { next; }                 # not the responsibility of this child or no valid reads found
		$n_mapped++;
		$result = find_circ(\%reads, $mismatch, $overhang, $string);
		if ($result !~ /Not a circRNA/) { print {$pipes{$process}} $result; }
	}

	close(FILE_IN);
	print {$pipes{$process}} "child finished " . $n_mapped . "\n";
	$pm->finish; 
}
				
# ----------------  child process ends  ------------------ #

# ---------------  parent process starts  ----------------- #
foreach $process (1 .. $process_n) {
	$pipes{$process}->reader();
	$select->add($pipes{$process});
}

while ($select->count() > 0) {            # all child processes that stil haven't finished yet
	@ready = $select->can_read(1);
	foreach (@ready) {
		$line = readline($_);
		if ($line =~ /child finished ([0-9]+)\n/) {
			$n_mapped += $1;
			$select->remove($_);
			next;
		}
		
		@items = split("\t", $line);                                    # items[0] is the name of the circRNA
		if (! exists $circ{$items[0]}) { $circ{$items[0]} = []; }
		push @{$circ{$items[0]}} , join("\t", @items[1 .. $#items]);    # print the other elements as they are 
	}					
}

open(FILE_OUT, "> " . $circ_file) or die "Cannot write to result file!\n";
print FILE_OUT "#Total mapped: " . $n_mapped . "\n";                    # total number of mapped reads
print FILE_OUT "chr\tstart\tend\tname\tstrand\toverhang_start\toverhang_end\tmutation_start\tmutation_end\tCIGAR\tMD\tSequence\n";
foreach $circRNA (keys %circ) {            # count how many reads each circRNA has
	$circRNA =~ /^(.*)_([0-9]+)_([0-9]+)_circ_/;
	$circRNA1 = $1 . "\t" . $2 . "\t" . $3;
	foreach (@{$circ{$circRNA}}) { print FILE_OUT $circRNA1 . "\t" . $_; }	
}
close(FILE_OUT);

# ---------------  parent process ends  ----------------- #

print "Done searching circRNA.\n";
exit;