#!/usr/bin/perl -w

use strict;
use Carp;

if (scalar @ARGV != 2) {  die "Usage: $0 <read_list> <kmer_list>\n"; }

my $readFile = $ARGV[0];
my $kmerFile = $ARGV[1];

my %kmer2reads  = ();
my %read2offset = ();
my @unmatched   = ();

open(KIN,$kmerFile) or croak "Cannot read from <$kmerFile>.\n";
while(<KIN>)
{
	chomp;                  # no newline
  s/#.*//;                # no comments
  s/^\s+//;               # no leading white
  s/\s+$//;               # no trailing white
  next unless length;     # anything left?
	
	$kmer2reads{$_} = ();
}
close(KIN);

open(RIN,$readFile) or croak "Cannot read from <$readFile>.\n";
while(<RIN>)
{
	chomp;									# no newline
  s/#.*//;                # no comments
  s/^\s+//;               # no leading white
  s/\s+$//;               # no trailing white
  next unless length;     # anything left?
	my $read = $_;
	my $match = 0;
	foreach my $kmer (keys %kmer2reads)
	{	
		my $offset = index($read,$kmer);		
		if ($offset != -1)
		{
			++$match;
			#print "kmer = $kmer\n";
			#print "$read has $kmer at $offset\n";
			push @{$kmer2reads{$kmer} },$read;
			$read2offset{$read} = $offset;
		}
	}
	if ($match == 0)
	{
		# read has no matches, dodgy
		push @unmatched,$read;
	}
	elsif ($match > 1)
	{
		croak "Read $read has matches to multiple kmers. Not sure what do.\n";
	}
}
close(RIN);

foreach my $k (keys %kmer2reads)
{
	#print "$k :\n";
	my $maxOffset = 0;
	foreach my $rr (@{$kmer2reads{$k}})
	{
		if ($read2offset{$rr} > $maxOffset)
		{
			$maxOffset = $read2offset{$rr};
		}
	}
	#print "maxoffset is $maxOffset\n";
	
	my $kmerPadding = pack("A$maxOffset","");
	print $kmerPadding,$k,"\n";

	foreach my $r (@{$kmer2reads{$k}})
	{
		my $readOffset = $maxOffset - $read2offset{$r};
		my $readPadding = pack("A$readOffset","");
		print $readPadding,$r,"\n";
	}
	print "\n";
}
print "\n";
print "Reads without kmer matches:\n";
foreach my $r (@unmatched)
{
	print $r,"\n"; 
}

