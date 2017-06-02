# ===================================== #
#                index.pl               #
# ===================================== #
# This perl script is to index the reference genome
# $input_file: the absolute path of reference genome fasta file
# $output_file: the absolute path of output index file
# test: perl /Users/minzhe/Documents/Project/Clirc/bin/index.pl /Users/minzhe/Documents/Project/Clirc/data/test/chr4.fa.masked /Users/minzhe/Documents/Project/Clirc/data/test/library/index.txt
#!/usr/bin/perl

my ($chr, $i, @items, $len, $pos, $next, $pre_len);
my $jump = 5000;
my ($input_file, $output_file) = @ARGV;

open(FILE_IN, $input_file) or die("Open file $input_file failed: $!");
$pos = tell FILE_IN;
open(FILE_OUT, "> " . $output_file);

while (<FILE_IN>) {
    if ($_ =~ />(.+)\n/) {
        @items = split("\t", $1);
        $chr = $items[0];
        $len = 0;
        $next = 0;
    } else {
        $pre_len = $len;
        $len += (length($_) - 1);

        if ($len - $jump * $next > 0) {
            print FILE_OUT $chr . "_" . ($jump * $next + 1) . "\t" . ($jump * $next + 1 - $pre_len) . "\t" . $pos . "\n";
            $next++;
        }
    }

    $pos = tell FILE_IN
}

close(FILE_IN);
close(FILE_OUT);

print "Done indexing genome.\n";
exit;