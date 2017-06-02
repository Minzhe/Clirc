# ===================================== #
#            Clirc_search.pl            #
# ===================================== #
# This perl script search CLIP-Seq reads spanning circRNAs on one or multiple fastq files. Each Clirc_search command processes one fastq file, which internally calls gsnap and can be very slow. The user can submit many jobs in parallel, each containing one Clirc_search command, to cut down the processing time. One thing to note is that each fastq file must have unique read names. Duplicate read names may arise when the user catenate two fastq files into one.
# arguments: see below help message.
# demo: perl bin/Clirc_search.pl -adaptor 5TGGC,TGAGATCGGAAGAGCGGTTCAGC -fastq data/test/test1.fastq -library data/test/library -results data/test/test1_results.txt -cleanup
# demo: perl bin/Clirc_search.pl -adaptor 5TGGC,TGAGATCGGAAGAGCGGTTCAGC -fastq data/test/test2.fastq -library data/test/library -results data/test/test2_results.txt -cleanup
#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;
use Cwd 'abs_path';
use File::Spec;

my ($adaptor, $library, $fastq, $sam, $overhang, $mismatch, $process_n);
my ($thread, $string, $cleanup, $help, $message, $results, $path);
my $flag = 0;

$help = 0;
$adaptor = "";
$cleanup = 0;
$process_n = 1;
$overhang = 5;
$mismatch = 0.15;
$string = 1;
$thread = 1;
$path = abs_path($0);
$path =~ s/Clirc_search\.pl$//;

GetOptions(
    "adaptor=s" => \$adaptor,
    "fastq=s" => \$fastq,
    "library=s" => \$library,
    "cleanup" => \$cleanup,
    "process_n=i" => \$process_n,
    "overhang=i" => \$overhang,
    "mismatch=f" => \$mismatch,
    "string=i" => \$string,
    "thread=i" => \$thread,
    "results=s" => \$results,
    "help" => \$help
) or die("Error in command line arguments.\n");

$fastq = File::Spec->rel2abs($fastq);
$library = File::Spec->rel2abs($library);
$results = File::Spec->rel2abs($results);

if ($help == 1) {
    $message = << 'HELP';
Clirc_search: search RBP-bound circRNAs in each SAM alignment file.

Options:
    -adaptor        adaptor sequence(s) to be trimmed from the fastq file. This is optional. Multiple adaptor sequences can be specified as a comma separated string. Simply specify a character string to be an adaptor sequnece to be trimmed from 3' end; specify a character string with "5" in the leading character to be an adaptor sequence to be trimmed from 5' end. Example "-adaptor TCGATCGA,5ATCGATCG" will trim the read ATCGATCGTTTAAAGGGCCCTCGATCGA
    -fastq          file name (including path) of the fastq file. Note: some intermediate files will be generated in the same folder of the fastq file when Clirc works
    -library        path to the directory where intermediate and gsnap index files are stored (same as the one used in Clirc_library)
    -process_n      number of parallel processes to start (default=1)
    -thread         number of threads for gsnap (default=1)
    -overhang       minimum overhang on both sides of the junction site for a CLIP-Seq read to be considered as a potential circRNA (default=5)
    -mismatch       maximum allowed mismatch rate between the CLIP-Seq read and the circRNA reference (default=0.15)
    -string         stringency of the analysis, 1 is more stringent (default), 0 is less stringent
    -results        output file name for analysis results
    -cleanup        whether to automatically delete intermediate files after analysis is complete
    -help           print this help message

HELP

    print $message;
    exit;
}

print "* Start searching RBP-bound circRNA in alignment file. *\n";

if ($adaptor ne "") {
    foreach (split(",", $adaptor)) {
        0 == system("perl " . $path . "/remove_adaptor.pl fastq " . $fastq . " " . $_ . " 18 0.2 " . $process_n . " unique\n")
            or die("Removing adaptor failed!\n");
        if ($flag == 1) { unlink($fastq); }
        $fastq .= ".removed";
        $flag = 1;
    }
}

$sam = $fastq . ".sam";
print("Calling gsnap, it could be slow ...\n");
0 == system("gsnap --no-sam-headers -B 4 -A sam -Q -n 5 --query-unk-mismatch=1 --merge-distant-samechr -N 1 -t " . $thread . " -m 4 --nofails --sam-multiple-primaries -D " . $library . " -d circRNA " . $fastq . " > " . $sam . " 2>&-")
    or die("gsnap alignment failed!\n");
0 == system("perl " . $path . "/find_circRNA.pl " . $sam . " " . $results . " " . $overhang . " " . $mismatch . " " . $process_n . " " . $string)
    or die("circRNA search failed!\n");

if ($cleanup == 1) {
    if ($flag == 1) { unlink($fastq); }
    unlink($sam);
}

print "* Done searching RBP-bound circRNA! *\n";