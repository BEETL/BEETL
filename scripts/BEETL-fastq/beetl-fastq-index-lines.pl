#!/usr/bin/perl -w

use strict;

my ($rzFile,$rzFileIndex);

my $lineNum=0;
my $linesPerEntry=4;
my $bitsPerBlock=10; # 1024
my $entriesPerBlock=(1<<$bitsPerBlock);
my $blockMask=$entriesPerBlock-1;
my $sizeSoFar=0;
my $thisLine;
my @index;
my ($blockNum, $blockPos);

my $usage = "Usage:\n$0 myfile.fastq.rz - index myfile.fastq.rz\n$0 myfile.fastq.rz i_1 [i_2...i_n] - extract entries i_1 [i_2..i_n] from fastq";

if (@ARGV==0)
{
    die $usage;
}

$rzFile=shift(@ARGV);
$rzFileIndex=$rzFile.".idx";

if ($rzFile eq "-h")
{
  print ${usage}, "\n";
  exit 0;
}

if (@ARGV==0)
{
    # index;

    print STDERR "indexing $rzFile\n";

    open (IDX, ">$rzFileIndex") || die "Failed to open $rzFileIndex: $!";
    open (ZIP, "cat $rzFile | zcat |") || die "Failed to open $rzFile: $!";
    while ($thisLine=<ZIP>)
    {
   	if (($lineNum&$blockMask)==0)
	{
	    print IDX "$lineNum $sizeSoFar\n";
	    ++$blockNum;
	}
#	    $thisLine=<ZIP>;
#	    $thisLine=<ZIP>;
#	    $thisLine=<ZIP>;
	$lineNum++;
	$sizeSoFar=tell(ZIP);
    }
    print STDERR "$0: Indexed $lineNum entries with $blockNum index points\n";
    close (IDX);
    close (ZIP);
}
else
{
    open (IDX, "$rzFileIndex") || die "Failed to open $rzFileIndex: $!";
    while ($thisLine=<IDX>)
    {
#	print $thisLine;
	chomp ($thisLine);
	die unless ($thisLine=~/\d+\ (\d+)/);
	push @index, $1;
    }
#    print STDERR "$0: Found ", scalar(@index), " entries in index\n";

    close (IDX);
    for my $entryNum (@ARGV)
    {
#	print "$entryNum\n";
	$blockNum=($entryNum>>$bitsPerBlock);
	$blockPos=($entryNum&$blockMask);
#	print STDERR "$0: Want entry $blockPos in block $blockNum\n";
#	print STDERR "$0: Execute razip $rzFile -d -c -b $index[$blockNum]\n"; 
	my $childPid = open (ZIP, "razip $rzFile -d -c -b $index[$blockNum] |") 
	    || die "Failed to open $rzFile: $!";

	$lineNum=0;
	while($lineNum!=$blockPos)
	{   
	    $thisLine=<ZIP>; die $lineNum if (!defined($thisLine));
#	    $thisLine=<ZIP>;
#	    $thisLine=<ZIP>;
#	    $thisLine=<ZIP>;
	    $lineNum++;
	}
	$thisLine=<ZIP>; print $thisLine;
#	$thisLine=<ZIP>; print $thisLine;
#	$thisLine=<ZIP>; print $thisLine;
#	$thisLine=<ZIP>; print $thisLine;
    kill(9,$childPid);
	close(ZIP);
	next;
    }
}



