# ===================================== #
#           Clirc_library.pl            #
# ===================================== #
# This perl script combine the results if multiple fastq files are used and to refine the search to generate final output.
# arguments: see below help message.
# demo: perl bin/Clirc_filter.pl -out data/test/test_results.txt -input data/test/test1_results.txt,data/test/test2_results.txt -library data/test/library -cleanup
#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;
use Cwd 'abs_path';
use File::Spec;
use File::Path;

my ($out, $overlap, $ratio, $len, $input, $library, $cleanup);
my (@input, $help, $message, $path);

$path = abs_path($0);
$path =~ s/Clirc_filter\.pl$//;
$help = 0;
$overlap = 2;
$ratio = 0.34;
$len = 20;
$cleanup = 0;

GetOptions(
	"out=s" => \$out,
	"overlap=i" => \$overlap,
	"ratio=f" => \$ratio,
	"len=i" => \$len,
	"input=s" => \$input,
	"library=s" => \$library,
	"cleanup" => \$cleanup,
	"help" => \$help,
) or die("Error in command line arguments\n");

$out = File::Spec->rel2abs($out);
$library = File::Spec->rel2abs($library);

if ($help == 1) {
	$message = << 'HELP';
Clirc_filter: filter Clirc_search results to generate final list of RBP-bound circRNAs.

Options:
	-out		output file name for final analysis results
	-overlap	minimum number of overlapping reads to call a true CLIP cluster on circRNA (default=2)
	-ratio		minimum number of unique start/end positions divided by all start/end positions (default=0.34)
	-len		minimum length of the CLIP cluster on circRNA (default=20)
	-input		comma separated list of analysis output generated by Clirc_search.pl (the -results parameter)
	-library	path to the directory where intermediate and gsnap index files are stored (same as the one used in Clirc_library)
	-cleanup	whether to automatically delete intermediate files after analysis is complete	
	-help		print this help message
HELP

	print $message;
	exit;
}

print "* Start combining and filtering RBP-bound circRNA in alignment file. *\n";

@input = split(",", $input);
foreach (@input) {
	if ($_ =~ s/^~//) { $_ = $ENV{"HOME"} . $_; }
	$_ = File::Spec->rel2abs($_);
}
$input = join(" ", @input);
0 == system("perl " . $path . "/filter_circRNA.pl " . $overlap . " " . $ratio . " " . $len . " " . $out . " " . $library . " " . $input)
    or die("filtering for circRNAs failed!\n");

if ($cleanup == 1) {
	map { unlink($_) } split(" ", $input);
	if (-e $library) { rmtree($library); }	
}

print "* Done filtering. *\n";