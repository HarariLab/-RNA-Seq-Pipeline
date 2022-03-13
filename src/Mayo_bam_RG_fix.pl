#!/usr/bin/env perl
# Mayo_bam_RG_fix.pl
# 2019_01_24
use strict;
use warnings;

# This file will take all the Mayo bams in a folder and parse their @RG header lines to create the picard AddOrReplaceReadGroups commands

my ($bam, $RG, $RGID, $RGLB, $RGPL, $RGPU, $RGSM, $extra, $preRGPU, $bamout, $command, $outfile, $time);
my @bams = ();

# First, save the names of the bams to an array
print "ls *.bam > km.bam.tmp\n";
system 'ls *.bam > km.bam.tmp';
open(IN, "km.bam.tmp") or die "Couldn't open your list of bams $!";
while (my $line = <IN>) {
  next if ($line =~/^\n/);
	chomp ($line);
  $bam = $line;
  push(@bams, $bam);
}

# Now run samtools view -H | grep '@RG' on each bam and parse the results, printing a log file and running the command
my $x=0;
my $count = scalar(@bams);
while ($count > $x) {
  $bam = $bams["$x"];
  $outfile = $bam; $outfile =~s/\.bam/_RG.fix.log/;
  $bamout = $bam; $bamout =~s/\.bam/_RG.fix.bam/;
  open (OUT, ">$outfile") or die "Couldn't create the file to hold your bam RG fix info $!";
  print OUT "Bam\tRGID\tRGLB\tRGPL\tRGPU\tRGSM\tcommand\n";
  print "samtools view -H $bam | grep \'\@RG\'\n";
  $RG = `samtools view -H $bam | grep \'\@RG\'`;
  #@RG	ID:1000_CER.FCC7KUNACXX_L5IATCACG	PL:Illumina	PU:pu	LB:lb	SM:sm
  ($extra, $RGID, $RGPL, $RGPU, $RGLB, $RGSM) = split(/\t/, $RG, 6);
  $RGID =~s/^ID://; print "\$RGID = $RGID\n";
  $RGPL =~s/^PL://;  print "\$RGPL = $RGPL\n";
  ($extra, $preRGPU) = split (/\./, $RGID, 2); print "\$preRGPU = $preRGPU\n";
  ($RGPU, $extra) =  split (/\_/, $preRGPU, 2); print "\$RGPU = $RGPU\n";
  ($RGSM, $extra) = split(/\./, $RGID, 2); print "\$RGSM = $RGSM\n";
  $RGLB = $RGSM;  print "\$RGLB = $RGLB\n";
  $command = "java -jar /opt/picard-tools/picard.jar AddOrReplaceReadGroups I=$bam O=$bamout RGID=$RGID RGLB=$RGLB RGPL=$RGPL RGPU=$RGPU RGSM=$RGSM";
  print OUT "$bam\t$RGID\t$RGLB\t$RGPL\t$RGPU\t$RGSM\t$command\n";
  $time = scalar(localtime(time + 0));
  print OUT "Start: $time\n";
  `$command`;
  $x++;
}
close(IN);
$time = scalar(localtime(time + 0));
print OUT "Stop: $time\n";
close(OUT);
`rm -rf km.bam.tmp`;

print "\nSapiens qui prospicit: Wise is he who looks ahead.  Your fixed bam is $bamout\n";
