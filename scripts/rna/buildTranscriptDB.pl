#!/usr/bin/perl -w

use strict;
use Carp;

# Example GFF entry:
# Chr1_supercontig_000000000      miRNA   exon    524299  524395  .       +       .        gene_id "ENSSHAG00000020030"; transcript_id "ENSSHAT00000023646"; exon_number "1";

my $reffile = "/illumina/scratch/denovo/results/taz/phusion2/assembly7.2_stfix/tdevil_supercontig-v7.2.fasta";

print STDERR "Reading genome file from $reffile.\n";
my %genome=();
open(FIN,$reffile) or die "Cannot read from <$reffile>.\n";
my $header="NA";
my $seq="";
my $ct=0;
while(<FIN>)
{	
	++$ct;
	chomp;
  if (m/^>/)
  {
  	if ($header ne "NA")
    {
				$genome{$header} = $seq;
    }
    s/^>//;
    $header=$_;
   	$seq="";
	}
  else
  {
  	$seq.=$_;
	}
}
close(FIN);
$genome{$header} = $seq;
print "Read " . (scalar keys %genome) . " scaffolds.\n";

my %strands      = ();
my %transcripts  = ();
my %splice_sites = ();

while(<STDIN>)
{
	chomp;
	my ($ctg,$source,$feature,$start,$end,$score,$strand,$frame,$attribute)=split(/\t/,$_);
	#print $attribute,"\n";

	next if ($feature ne "CDS"); # skip annotations which are not CDS

	my @p=split(/\s+/,$attribute);
	my $t_id = $p[4]; # transcript id
	$t_id =~ s/\"//g;
	$t_id =~ s/;//g;
	#print $t_id,"\n";

	if (! exists($transcripts{$t_id}))
	{
		#print "New transcript $t_id\n";
		$strands{$t_id} = $strand;
		if ( !exists $genome{$ctg} ) 
		{
			die "$ctg of annotation $t_id not fine in genome file.\n";
		}
		my $len = ($end-$start);
		#print "Retrieving exon of length $len\n";
		$transcripts{$t_id} = substr ($genome{$ctg},$start,$len);
	}
	else
	{
		if ($strand ne $strands{$t_id})
		{
			die "Strand of exon in $t_id is $strand but expected ",$transcripts{$t_id},"\n";
		}
		my $len = ($end-$start);
		#print "Retrieving exon of length $len\n";
		if ( exists $splice_sites{$t_id})
		{
			push @{ $splice_sites{$t_id} }, length($transcripts{$t_id});
		}
		else
		{
			$splice_sites{$t_id} = ();
			push @{ $splice_sites{$t_id} }, length($transcripts{$t_id}); 
		}
		$transcripts{$t_id} .= substr ($genome{$ctg},$start,$len);
	}
}
my $tfile = "transcriptDB.fa";
open(OUT,">",$tfile) or croak "Cannot write to <$tfile>.\n";
foreach my $t (sort {$a cmp $b} keys %transcripts)
{
	print OUT ">",$t,"_",$strands{$t},"\n";
	print OUT $transcripts{$t},"\n";
}
close(OUT);

my $spsite_file = "transcript_splice_sites.txt";
open(OUT,">",$spsite_file) or croak "Cannot write to <$spsite_file>.\n";
foreach my $v (sort {$a cmp $b} keys %splice_sites)
{
	print OUT $v,"_",$strands{$v}," : ";
	if (scalar @{$splice_sites{$v} } > 0) 
	{
		print OUT join(" ",@{$splice_sites{$v}});
	}
	print OUT "\n";
}
close(OUT);

print "Keys: " . (scalar keys %splice_sites) . " " . (scalar keys %transcripts) . " " . (scalar keys %strands) . "\n";

