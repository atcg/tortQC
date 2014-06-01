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
my $reference = "/mnt/Data4/genomes/Galapagos.fasta"; # Default to galapagos tortoise

GetOptions  ("reads=s"         => \$readsDir,
             "adapters=s"      => \$adaptersDir,
             "out=s"           => \$outDir,
             "log=s"           => \$logFile,
             "threads=i"       => \$threadsMax,
             "trim"            => \$trim,
             "map"             => \$map,
             "reference=s"     => \$reference,
             "help|man" => \$help) || pod2usage(2);

if (!$readsDir or !$adaptersDir or !$outDir or $help) {
    pod2usage(-exitval => 0, -verbose => 2, -noperldoc => 1);
}


my $startingDir = getcwd();
unless (-d $outDir) {
    mkdir $outDir or die "Couldn't make output directory $outDir: $!\n";
}
open(my $logFH, ">", $logFile) or die "Couldn't open log file $logFile for writing: $!\n";

my $fqjDir = $outDir . "/fastq-join";
unless (-d $fqjDir) {
    mkdir $fqjDir;
}

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
    

    print $logFH "Running fastq-join to merge overlapping paired end reads\n";
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
    my $mappingDir = $startingDir . "/mapping";
    unless (-d $mappingDir) {
        mkdir $mappingDir;
    }
    foreach my $readGroup (sort keys %sampleNamesHash) {
        my $singlesSamFile = $mappingDir . "/" . $readGroup . ".singlesAndJoined.sam";
        my $singlesBamFile = $mappingDir . "/" . $readGroup . ".singlesAndJoined.bam";
        my $pairedSamFile = $mappingDir . "/" . $readGroup . ".paired.sam";
        my $pairedBamFile = $mappingDir . "/" . $readGroup . ".paired.bam";
        my $mergedBamFile = $mappingDir . "/" . $readGroup . "_merged.bam";
        my $reads1 = $fqjDir . "/" . $readGroup . "_trimmed_fqj.un1.fastq";
        my $reads2 = $fqjDir . "/" . $readGroup . "_trimmed_fqj.un2.fastq";
        my $readsSingles = $fqjDir . "/" . $readGroup . "_R1andR2trimmomaticSingles.fastq";
        system("bwa mem -t 12 $reference $readsSingles > $singlesSamFile");
        system("bwa mem -t 12 $reference $reads1 $reads2 > $pairedSamFile");
        
        system("samtools view -b -S $singlesSamFile > $singlesBamFile");
        system("samtools view -b -S $pairedSamFile > $pairedBamFile");
        
        system("samtools merge $mergedBamFile $singlesBamFile $pairedBamFile");
        
        # Mark duplicates and use mpileup
        my $cleanedBam = $mappingDir . "/" . $readGroup . ".merged.cleaned.bam";
        my $sortedBam = $mappingDir . "/" . $readGroup . ".merged.cleaned.sorted.bam";
        my $markDupsBam = $mappingDir . "/" . $readGroup . ".merged.cleaned.sorted.markDups.bam";
        my $markDupsMetrics = $mappingDir . "/" . $readGroup . ".merged.sorted.cleaned.markDups.metrics";
        my $pileupFile = $mappingDir . "/" . $readGroup . ".mpileup";
        system("java -jar ~/bin/picard/picard-tools/CleanSam.jar I=$mergedBamFile O=$cleanedBam");    
        system("java -jar ~/bin/picard/picard-tools/AddOrReplaceReadGroups.jar I=$cleanedBam O=$sortedBam SORT_ORDER=coordinate RGPL=illumina RGPU=Test RGLB=Lib1 RGID=$readGroup RGSM=$readGroup VALIDATION_STRINGENCY=LENIENT");
        system("java -jar ~/bin/picard/picard-tools/MarkDuplicates.jar I=$sortedBam O=$markDupsBam METRICS_FILE=$markDupsMetrics MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=250 ASSUME_SORTED=true REMOVE_DUPLICATES=false");
        system("samtools mpileup $markDupsBam > $pileupFile");
    
        print "Stats for $readGroup:\n";
        system("samtools flagstat $markDupsBam");
        print "\n\n\n";
        # All we need to keep is the $markDupsBam file, so get rid of the other sams and bams
        unlink ($singlesSamFile, $pairedSamFile, $singlesBamFile, $pairedBamFile, $mergedBamFile, $cleanedBam, $sortedBam,);
        
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
   -trim            Perform read trimming using Trimmomatic and join with fastq-join
   -map             Perform mapping to reference genome
   -threads         Threads to use for multithreading (default=4)
   -reference       Reference to use for mapping (default to /mnt/Data4/genomes/Galapagos.fasta)
   -help|man        Prints out documentation


=head1 DESCRIPTION

This was written for the purposes of QCing HiSeq reads for a desert tortoise
project

=cut










