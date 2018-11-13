![Clirc_logo](QBRC.jpg)

# Clirc

Clirc: a user-friendly bioinformatics software to identify RBP-bound circRNAs through analysis of CLIP-Seq data

## Description

The Clirc software is designed to identify RBP-bound circRNAs in CLIP-Seq data. The Clirc software internally calls the gsnap aligner for alignment of CLIP-Seq reads. The CLIP-Seq technology comprises of three variants: HITS-CLIP, PAR-CLIP, and iCLIP. The Clirc software works on all three versions of them. It can also work on any type of RNA-Seq data in general, although non-polyA selected RNA-Seq is preferred. The Clirc software is designed to search for known circRNAs in high-throughput sequencing data whose coordinates are provided by users. The users can collect these coordinates from previous publication or their own search using RNA-Seq data. This avoids the common problem of short sequncing reads in CLIP-Seq data. But this approach cannot identify novel circRNAs in CLIP-Seq data.
 
For the user's convenience, we also provided comprehensive circRNA coordinates for hg19, dm3 and mm9 in our data folder. These coordinates are the ones used in our own publication.

Clirc is designed to be used on Unix or Unix-like systems only.

Cite our paper:


## Install

1. A recent version of **Perl** is needed.

The user may also need to install `List::MoreUtils` and `Parallel::ForkManager` if these two modules are not installed yet. A txt file called MODULE is provided in the doc folder if the user doesn't know how to install Perl modules.

2. A recent version of **gsnap** is needed. 

## Usage

The whole Clirc software contains 3 commands in the bin directory to run for identifying RBP-bound circRNAs from CLIP-Seq data.
```{}
Clirc_library	build alignment library for identifying RBP-bound circRNAs from CLIP-Seq data
Clirc_search	search RBP-bound circRNAs in each SAM alignment file
Clirc_filter	filter Clirc_search results to generate final list of RBP-bound circRNAs
```
Please use `perl Clirc_library.pl -help`, `perl Clirc_search.pl -help`, and `perl Clirc_filter.pl -help` to show the detailed instructions for running each command. Here I will describe the general workflow of Clirc

1. Run `Clirc_library` to build a pseudo reference genome containing the regular reference genome together with user-supplied circRNAs
2. Run `Clirc_search` on one or multiple fastq files to identify CLIP-Seq reads spanning circRNAs. Each Clirc_search command processes one fastq file, which internally calls gsnap and can be very slow. The user can submit many jobs in parallel, each containing one `Clirc_search` command, to cut down the processing time. One thing to note is that each fastq file must have unique read names. Duplicate read names may arise when the user catenate two fastq files into one.
3. Run `Clirc_filter` to combine the results if multiple fastq files are used and to refine the search to generate final output

## DEMO

To demonstrate the usage of Clirc, we put two test fastq files in the data/test/ folder. And here are the demo commands to analyze them:
```{perl}
perl bin/Clirc_library.pl -coord data/test/dm3_chr4_circRNAs.txt -genome data/test/chr4.fa.masked -library data/test/library
perl bin/Clirc_search.pl -adaptor 5TGGC,TGAGATCGGAAGAGCGGTTCAGC -fastq data/test/test1.fastq -library data/test/library -results data/test/test1_results.txt -cleanup
perl bin/Clirc_search.pl -adaptor 5TGGC,TGAGATCGGAAGAGCGGTTCAGC -fastq data/test/test2.fastq -library data/test/library -results data/test/test2_results.txt -cleanup
perl bin/Clirc_filter.pl -out data/test/test_results.txt -input data/test/test1_results.txt,data/test/test2_results.txt -library data/test/library -cleanup
```

## Contact

Minzhe Zhang <<minzhe.zhang@utsouthwestern.edu>>

## Version update
0.1.0: First release. (11/12/2018)
