#!/usr/bin/perl -w

# Author: Tony Cox 17th February 2012

# implicit_sort.pl
# This takes BWT data as output by the Issue002 branch version of
# BCRext. Needs to be run in ASCII mode
# /path/Beetl ext -i ../test.txt -a
# A lower case character indicates the character is in the same
# SAP-interval as the previous character. A 'D' implies similarly
# for $ characters (but this isn't used so they are converted
# straight back to $ characters) 


# TBD #1: Buffer symbols, don't read them all into RAM
# TBD #2: Do a minimal permutation - only permute when number of runs in
# interval exceeds number of characters present
# TBD #3: look more closely at making the first run in an SAP-interval match
# the last run of a preceding interval (a preliminary look at this
# indicated it doesn't help compression very much)


$numberPermuted=0;

while(<STDIN>)
{
    $a.=$_
}

$l=length($a);
$a.="N"; # this triggers printing out the final interval
for ($i=0;$i<=$l;$i++)
{
    $c=substr($a,$i,1);
    $c='$' if ($c eq 'D');
    if ($c=~/[acgtn]/)
    {
	$c=uc($c);
	$h{$c}++;
    }
    else
    {
	# process any preceding interval
	for $k (keys %h)
	{
	    $numberPermuted+=$h{$k}*(scalar(keys %h)>1);
	    print $k x $h{$k};
	}    
	%h=();
	$h{$c}++;
    }
}
print STDERR "Permuted at most $numberPermuted out of $l\n";
