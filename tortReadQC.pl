#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use Cwd;
use Parallel::ForkManager;


my $help = 0;
my $readsDir;
my $adaptersDir;
my $outDir;
my $logFile;
my $threadsMax = 4;
my $trim;
my $map;

GetOptions  ("reads=s"         => \$readsDir,
             "adapters=s"      => \$adaptersDir,
             "out=s"           => \$outDir,
             "log=s"           => \$logFile,
             "threads=i"       => \$threadsMax,
             "trim"            => \$trim,
             "map"             => \$map,
             "help|man" => \$help) || pod2usage(2);

if (!$readsDir or !$adaptersDir or !$outDir or $help) {
    pod2usage(-exitval => 0, -verbose => 2, -noperldoc => 1);
}


my $startingDir = getcwd();
unless (-d $outDir) {
    mkdir $outDir or die "Couldn't make output directory $outDir: $!\n";
}
open(my $logFH, ">", $logFile) or die "Couldn't open log file $logFile for writing: $!\n";


# First gather up all the file names of the files in the specified reads directory
my %readFilesHash;
my %sampleNamesHash;
chdir($readsDir);
opendir my $readsDirFH, "./";
my @readsFiles = readdir $readsDirFH or die "Couldn't readdir $readsDirFH: $!\n";
closedir $readsDirFH;
chdir $startingDir;

my $trimmomaticDir = $outDir . "/trimmomatic";
unless (-d $trimmomaticDir) {
    mkdir $trimmomaticDir or die "Couldn't make directory $trimmomaticDir: $!\n";
}

# Now do sequence QC using Trimmomatic
{
my @trimmomaticCommands;
print $logFH "Generating trimmomatic commands on all read files.\n";
foreach my $file (@readsFiles) {
    if ($file =~ /(.*)_(.*)_R1.fastq.gz/) { # Only looking for the R1s here, because we only want one trimmomatic run per library (or in the case of too many reads for a single fastq.gz file, do a trimmomatic run for each pair of read files)
        # $1 = sample name
        # $2 = index sequence
        print $logFH "--------------------------------------------------\n";
        print $logFH "Sequence file found: $1\_$2\_R1.fastq.gz\n";
        my $R1File = $readsDir . "$1\_$2\_R1.fastq.gz";
        print $logFH "R1 file: $readsDir" . "$1\_$2\_R1.fastq.gz\n";
        my $R2File = $readsDir . "$1\_$2\_R2.fastq.gz";
        print $logFH "R2 file: $readsDir" . "$1\_$2\_R2.fastq.gz\n";
        my $readGroupName = $1 . "_" . $2;
        
        my $adaptersFile = $adaptersDir . "$1\.adapters";
        print $logFH "Adapters file to be used for the sequence group: " . $adaptersFile . "\n";
        
        my $R1OutFilePaired = "$trimmomaticDir/$1\_$2\_R1_paired_trimmed.fastq.gz";
        print $logFH "R1 paired trimmomatic output file: $trimmomaticDir/$1\_$2\_R1_paired_trimmed.fastq.gz\n";
        my $R1OutFileSingles = "$trimmomaticDir/$1\_$2\_R1_singles_trimmed.fastq.gz";
        print $logFH "R1 singles trimmomatic output file: $trimmomaticDir/$1\_$2\_R1_singles_trimmed.fastq.gz\n";
        my $R2OutFilePaired = "$trimmomaticDir/$1\_$2\_R2_paired_trimmed.fastq.gz";
        print $logFH "R2 paired trimmomatic output file: $trimmomaticDir/$1\_$2\_R2_paired_trimmed.fastq.gz\n";
        my $R2OutFileSingles = "$trimmomaticDir/$1\_$2\_R2_singles_trimmed.fastq.gz";
        print $logFH "R2 singles trimmomatic output file: $trimmomaticDir/$1\_$2\_R2_singles_trimmed.fastq.gz\n";
        $sampleNamesHash{$readGroupName}{'R1_paired_trimmed'} = $R1OutFilePaired;
        $sampleNamesHash{$readGroupName}{'R1_singles_trimmed'} = $R1OutFileSingles;
        $sampleNamesHash{$readGroupName}{'R2_paired_trimmed'} = $R2OutFilePaired;
        $sampleNamesHash{$readGroupName}{'R2_singles_trimmed'} = $R2OutFileSingles;        
        push (@trimmomaticCommands, "java -Xmx8G -jar ~/bin/Trimmomatic-0.32/trimmomatic-0.32.jar PE -threads 2 -phred33 $R1File $R2File $R1OutFilePaired $R1OutFileSingles $R2OutFilePaired $R2OutFileSingles ILLUMINACLIP:$adaptersFile:2:30:10 LEADING:5 TRAILING:15 SLIDINGWINDOW:4:20 MINLEN:40");
    }
}
print $logFH "--------------------------------------------------\n";
print $logFH "Finished generating trimmomatic commands on all read files.\n\n\n";


if ($trim) {
    print $logFH "Running all trimmomatic commands\n";
    my $counter = 0;
    my $forkManager = new Parallel::ForkManager($threadsMax);
    foreach my $trimCommand (@trimmomaticCommands) {
        $counter++;
        print $logFH "--------------------------------------------------\n";
        print $logFH "Trimmomatic command $counter: \n\t"; # Indent the next line to make it easier to find the commands in the text
        print $logFH $trimCommand . "\n";
        sleep 10;
        print "\n";
        $forkManager->start and next;
        print "\n";
        system("$trimCommand");
        print "Finished running the following:\n\t$trimCommand\n\n";
        $forkManager->finish;
    }
    $forkManager->wait_all_children;
    print $logFH "--------------------------------------------------\n";
    print $logFH "Finished running all trimmomatic commands\n";
    print $logFH "--------------------------------------------------\n\n";
    }
    
    my $fqjDir = $outDir . "/fastq-join";
    unless (-d $fqjDir) {
        mkdir $fqjDir;
    }
    
    # print $logFH "Running fastq-join to merge overlapping paired end reads\n";
    foreach my $readGroup (sort keys %sampleNamesHash) {
        sleep 10;
        print $logFH "\n\nStarting to process $readGroup through the fastq-join and assembly stages\n\n";
        my $combinedTrimmomaticSingles = $fqjDir . "/" . $readGroup . "_R1andR2trimmomaticSingles.fastq.gz";
        my $R1pairedTrimmedUnzipped;
        my $R2pairedTrimmedUnzipped;
        if ($sampleNamesHash{$readGroup}{'R1_paired_trimmed'} =~ /(.*).gz$/) {
            $R1pairedTrimmedUnzipped = $1;
            system("gunzip $sampleNamesHash{$readGroup}{'R1_paired_trimmed'}");
        }
        if ($sampleNamesHash{$readGroup}{'R2_paired_trimmed'} =~ /(.*).gz$/) {
            $R2pairedTrimmedUnzipped = $1;
            system("gunzip $sampleNamesHash{$readGroup}{'R2_paired_trimmed'}");
        }
        my $fqjOutputPrefix = $fqjDir . "/" . $readGroup . "_trimmed_fqj.%.fastq";
        system("fastq-join -v ' ' $R1pairedTrimmedUnzipped $R2pairedTrimmedUnzipped -o $fqjOutputPrefix");
        system("cat $sampleNamesHash{$readGroup}{'R1_singles_trimmed'} $sampleNamesHash{$readGroup}{'R2_singles_trimmed'} > $combinedTrimmomaticSingles");
        system("gunzip $combinedTrimmomaticSingles");
        my $joinedReads = $fqjDir . "/" . $readGroup . "_trimmed_fqj.join.fastq";
        my $joinedAndSingles = $fqjDir . "/" . $readGroup . "_combinedJoinedAndSingles.fastq";
        # Need to update this file name since we unzipped the file above
        if ($combinedTrimmomaticSingles =~ /(.*).gz$/) {
            $combinedTrimmomaticSingles = $1;
        }
        
        system("cat $joinedReads $combinedTrimmomaticSingles > $joinedAndSingles");
    }
    
    print $logFH "--------------------------------------------------\n";
    print $logFH "Finished running fastq-join on all samples\n";
    print $logFH "--------------------------------------------------\n\n";
}

if ($map) {
    foreach my $readGroup (sort keys %sampleNamesHash) {
        
        
        
        
    }
}












#Documentation
__END__

=head1 NAME

tortReadQC.pl

=head1 SYNOPSIS 

perl tortReadQC.pl --reads <file> --adapters <file>

 Options:
   -reads=s         Directory with raw reads in gzipped fastq format
   -adapters=s      Directory with adapters fasta files
   -out             Name of output directory (it will be created)
   -log             Name of logfile to print output (you will probably also want
                    to capture STDERR manually)
   -threads         Threads to use for multithreading (default=4)
   -help|man        Prints out documentation


=head1 DESCRIPTION

This was written for the purposes of QCing HiSeq reads for a desert tortoise
project

=cut