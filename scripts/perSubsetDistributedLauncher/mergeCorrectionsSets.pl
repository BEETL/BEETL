use warnings FATAL => 'all';
use strict;
use Cwd qw(abs_path);
use POSIX qw(strftime);
use IO::File;
use Carp;
use Cwd;
use Switch;
use List::Util qw( min max );
use Getopt::Long;

open (my $corrections, "<", $ARGV[0]);

my $firstLine = 1;
my $newLongestWitness;
my $newShortestWitness;
my $newRead;
my $newCorrectionString;
my $newPosition;
my $newStrand;
my $highestLongestWitness;
my $lowestShortestWitness;
my $longestCorrectionString = "";
my $lastRead = -1;
my $lastStrand = -1;
my $lastPosition = -1;
my $correctorStart;

print "read position reverse_strand correction corrector_start shortest_witness longest_witness" . "\n";

my %colNames = ();

foreach my $line (<$corrections>) {
    $line =~ s/\s+$//;

    my @parts = split(/ /,$line);
    if ($firstLine) {
        $firstLine = 0;

        for (my $index = 0; $index < scalar @parts; $index++) {
            $colNames{$parts[$index]} = $index;
        }
    } 
    else {
        $newRead = $parts[$colNames{"read"}];
        $newStrand = $parts[$colNames{"reverse_strand"}];
        $newPosition = $parts[$colNames{"position"}];

        if ($newRead != $lastRead or $newStrand != $lastStrand or $newPosition != $lastPosition) {
            
            if ($lastStrand == 0) {
                $correctorStart = $lastPosition - length($longestCorrectionString) + 1;
            } else {
                $correctorStart = $lastPosition;
            }

            if ($lastRead >= 0) {
                print $lastRead . " ";
                print $lastPosition . " ";
                print $lastStrand . " ";
                print $longestCorrectionString . " ";
                print $correctorStart . " ";
                print $lowestShortestWitness . " ";
                print $highestLongestWitness . "\n";
            }

            $highestLongestWitness = $parts[$colNames{"longest_witness"}];
            $lowestShortestWitness = $parts[$colNames{"shortest_witness"}];
            $longestCorrectionString = $parts[$colNames{"correction"}];

        }
        else
        {
            if ($parts[$colNames{"longest_witness"}] > $highestLongestWitness)
            {
                $highestLongestWitness = $parts[$colNames{"longest_witness"}];
            }
            if ($parts[$colNames{"shortest_witness"}] < $lowestShortestWitness)
            {
                $lowestShortestWitness = $parts[$colNames{"shortest_witness"}];
            }
            if (length($parts[$colNames{"correction"}]) > length($longestCorrectionString))
            {
                $longestCorrectionString = $parts[$colNames{"correction"}];
            }
        }

        $lastRead = $newRead;
        $lastStrand = $newStrand;
        $lastPosition = $newPosition;
    }
}

print $lastRead . " ";
print $lastPosition . " ";
print $lastStrand . " ";
print $longestCorrectionString . " ";
print $correctorStart . " ";
print $lowestShortestWitness . " ";
print $highestLongestWitness . "\n";
