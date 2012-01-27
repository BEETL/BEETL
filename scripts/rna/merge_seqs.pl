#!/usr/bin/perl -w

use strict;
use Carp;
use Getopt::Long;

my $infile     = "count.txt";
my $outfile    = "consensus.fastq";
my $outfile_r1 = "consensus_read1.fastq";
my $outfile_r2 = "consensus_read2.fastq";
my $verbose    = 0;
my $n_args     = @ARGV;

my $ret = GetOptions("infile=s"     => \$infile,      
										 "outfile=s"    => \$outfile,      
										 "outfile_r1=s" => \$outfile_r1,      
										 "outfile_r2=s" => \$outfile_r2,      
                     "verbose"      => \$verbose);  


#print "$ret $n_args\n"; 
if (!$ret || !$n_args) { croak "Usage: $0 infile=<reads_bkpt_assigned.txt> outfile=<consensus.fastq> outfile_r1=<consensus_r1.fastq> outfile_r2=<consensus_r2.fastq>\n"; }

print STDERR "Reading from <$infile>.\n";
print STDERR "Writing joint consensus sequence to <$outfile>.\n";
print STDERR "Writing R1 consensus sequence to <$outfile_r1>.\n";
print STDERR "Writing R2 consensus sequence to <$outfile_r2>.\n";

# stitch sequence of longest breakpoint and read
sub stitchBreakpoint($$);

my $prev_bseq  = "NA";
my $prev_bext  = "NA";
my @bkpt_reads = (); # stores the reads associated with the current breakpoint 

my $ct=0; # counter for consensus sequences

open(OUT,">$outfile") or croak "Cannot write to <$outfile>.\n";
open(R1_OUT,">$outfile_r1") or croak "Cannot write to <$outfile_r1>.\n";
open(R2_OUT,">$outfile_r2") or croak "Cannot write to <$outfile_r2>.\n";
open(IN,$infile) or croak "Cannot read from <$infile>.\n";
while(<IN>)
{
	next if (m/^NA/); # skip reads without associated breakpoint
	chomp;
	# Example line
	# AAAAAAAAAAAAAAAAAAAAAAAAATGTTTGATC BKPT CTAGTTTGTAAAAAAAAAAAAAAAAAAAAAAAAA 1:9:0:0:0:0 1:0:0:0:0:11 -> AAAAAAAAAAAAAAAAAAAAAAAAATGTTTGATCAAAGTAA READ AATGAAACTAGTTTGTAAAAAAAAAAAAAAAAAAAAAAAAA 2:0:0:5:0:0
	my ($bseq_rev,$btag,$bseq,$bext_tumour,$bext_normal,$arrow,$read_rev,$rtag,$read_fwd,$read_ext)=split;

	if ($rtag ne "READ" || $btag ne "BKPT") { carp "Weird line: $_\n"; }

	# check if previous breakpoint sequence is substring of the current one
	my $l=length($bseq);
	my $tmp = substr($bseq,0,$l);

	#if ($prev_bseq eq $tmp)
	#{
		#print "GOTCHA!!\n";
		#print "prev_seq=$prev_bseq\n";
		#print "tmp=$tmp\n";
		#print "bseq=$bseq\n";
		#push @bkpt_reads, $read_fwd;
	#}	
	if ($prev_bseq ne "NA")
	{
		# new breakpoint
		print "$prev_bseq $prev_bext\n";
		#my $bseq_ext = extendBreakpoint($prev_bseq,$prev_bext);
		#print "EXTENDED_BREAKPOINT : $bseq_ext\n";

		# print reads sorted by length
		my $longest="NA";
		foreach my $r ( sort { length($a) <=> length($b) }  @bkpt_reads)
		{
			print "\t",$r,"\t","\n";
			$longest=$r;
		}
		my ($prefix,$suffix) = stitchBreakpoint($longest,$prev_bseq);
		print "CONSENSUS : $longest $prefix $suffix $prev_bseq\n";

		my $q    =  ('a') x length($longest);
		my $q_r1 =  ('a') x length($prefix);
		my $q_r2 =  ('a') x length($suffix);

		# write joint consensus sequence (this might not longer be needed but just in case)
		print OUT "\@$ct\n";
		print OUT "$longest\n";
		print OUT "+\n";
		print OUT $q,"\n";

		# write breakpoint prefix (read 1)
		print R1_OUT "\@$ct\n";
		print R1_OUT "$prefix\n";
		print R1_OUT "+\n";
		print R1_OUT $q_r1,"\n";

		# write breakpoint suffix (read 2)
		print R2_OUT "\@$ct\n";
		print R2_OUT "$suffix\n";
		print R2_OUT "+\n";
		print R2_OUT $q_r2,"\n";
	
		print "\n";
		# clear read set
		$#bkpt_reads = -1;
		++$ct;
	}
	
	push @bkpt_reads, $read_fwd;

	$prev_bseq = $bseq;
	$prev_bext = $bext_tumour; 
}
close(IN);
close(OUT);
close(R1_OUT);
close(R2_OUT);

sub stitchBreakpoint($$)
{
	my ($read,$breakpoint) = @_;
	#print "StitchBreakpoint : $read $breakpoint\n";	

	my $pref = "NA"; # sequence before the breakpoint
	my $suf  = "NA"; # sequence before the breakpoint

	my $res = index($read,$breakpoint);
	#print "res=$res\n";	
	if ($res != -1)
	{
		# breakpoint sequence is contained in read, all good
		$pref = substr($read,0,$res);
		$suf  = substr($read,$res);
	}
	else
	{
		# should not happen
		print "Read = $read\n";
		print "Breakpoint = $breakpoint\n";
		croak "Fatal error: breakpoint sequence not contained in read...\n";
	}
	return ($pref,$suf);
}

