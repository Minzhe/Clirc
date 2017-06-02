# =========================================== #
#             remove_adaptor.pl               #
# =========================================== #
# This perl script is to remove adaptor sequence.
# $format: fastq or fasta
# $file: sequencing file
# $adaptor: adaptor sequence
# $len_limit: if length of the sequence after adaptor removal is smaller than this limit, don't keep it
# $mismatch: fraction of mismatch allowed (0-1)
# $process_n: number of concorrent processes to run
# $unique: unique or all, whether to keep unique reads only or keep all reads
# parent outputs processed reads, child process 1, 2, ... trims adaptor
# test: perl bin/remove_adaptor.pl fastq data/test/test1.fastq 5TGGC 18 0.2 1 unique
# test: s
#!/usr/bin/perl

use strict;
use warnings;
use Parallel::ForkManager;
use IO::Pipe;
use IO::Select;
use Cwd 'abs_path';
my $path;
BEGIN {
	$path = abs_path($0);
	$path =~ s/remove_adaptor\.pl//;
}
use lib $path; 
use remove_adaptor_functions;
$SIG{__WARN__} = sub { die "Sequencing file format error!\n" };
$SIG{CHLD} = 'IGNORE';

my ($format, $file, $adaptor, $len_limit, $mismatch, $process_n, $unique) = @ARGV;
my ($select, @select, $process, $pid, %pipes, $pm);
my (@read, %regex, $line, $count, @ready, %unique);
my $five_prime = 0;

# -----------------  prepare for forking  ---------------- #

if ($adaptor =~ s/^5//) {                                       # adaptor is at 5' end
	$five_prime = 1;
	$adaptor = scalar reverse $adaptor;
} 
prepare_adaptor($adaptor, \%regex);                             # prepare adaptor hash
$pm = Parallel::ForkManager->new($process_n);                   # prepare for fork
$select = IO::Select->new();                                    # communication back to parent
foreach (1 .. $process_n) { $pipes{$_} = IO::Pipe->new(); } 


# -----------------  child process starts  ------------------- #
foreach $process (1 .. $process_n) {	
	$pid = $pm->start and next; # do the fork

	# establish communication with parent process
	$pipes{$process}->writer();
	$pipes{$process}->autoflush(1);

	open(FILE_IN, $file) or die "Cannot open file!\n";
	$count = 0;

	while ($read[0] = <FILE_IN>) {
		# read each read
		$count++;
		$read[1] = <FILE_IN>;
		if ($format eq "fastq") {
			$read[2] = <FILE_IN>;
			$read[3] = <FILE_IN>;
		}

		# process each sequencing read
		if (($count-1) % $process_n == $process-1) {
			if ($five_prime == 1) { reverse_read(\@read, $format); } # if adaptor is at 5'end
			trim_read(\@read, \%regex, $format, $mismatch);
			if ($five_prime == 1) { reverse_read(\@read, $format); }
			if (length($read[1]) > $len_limit) { print {$pipes{$process}} @read; }
		}
		if ($count == $process_n) { $count = 0; }
	}

	close(FILE_IN);
	print {$pipes{$process}} "Last read sent\n";
	$pm->finish; 
}
# ----------------  child process ends  ------------------ #

# ---------------  parent process starts  ------------------ #
foreach $process (1 .. $process_n) {
	$pipes{$process}->reader();
	$select->add($pipes{$process});
}

open(FILE_OUT, "> " . $file . ".removed") or die "Cannot write to file!\n";

while ($select->count()>0) {                        # all child processes that have processed output
	@ready = $select->can_read(1);

	foreach (@ready) {
		# read reads
		$read[0] = readline($_);
		if ($read[0] =~ /^Last read sent/) {
			$select->remove($_);
			next;
		}

		$read[1] = readline($_);
		if ($format eq "fastq") {
			$read[2] = readline($_);
			$read[3] = readline($_);
		}

		# print reads
		if ((! exists $unique{$read[1]}) || $unique eq "all") {
			print FILE_OUT @read;
			$unique{$read[1]} = 1;
		}
	}
}

close(FILE_OUT);

#-------  parent process ends  --------------------#

print "Done removing adaptors.\n";
exit;