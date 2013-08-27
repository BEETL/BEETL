#!/usr/bin/env perl

=head1 LICENSE

Copyright (c) 2011 Illumina, Inc.


This software is covered by the "Illumina Non-Commercial Use Software
and Source Code License Agreement" and any user of this software or
source file is bound by the terms therein (see accompanying file
Illumina_Non-Commercial_Use_Software_and_Source_Code_License_Agreement.pdf)

This file is part of the BEETL software package.

Citation: Markus J. Bauer, Anthony J. Cox and Giovanna Rosone
Lightweight BWT Construction for Very Large String Collections.
Proceedings of CPM 2011, pp.219-231


=head1 NAME

readEndPosFile.pl

=head1 DIAGNOSTICS

=head2 Exit status

0: successful completion
1: abnormal completion
2: fatal error

=head2 Errors

All error messages are prefixed with "ERROR: ".

=head2 Warnings

All warning messages are prefixed with "WARNING: ".

=head1 CONFIGURATION AND ENVIRONMENT

=back

=head1 BUGS AND LIMITATIONS

There are no known bugs in this module.

All documented features are fully implemented.

Please report problems to Illumina Technical Support (support@illumina.com)

Patches are welcome.

=head1 AUTHOR

Lilian Janin

=cut

use warnings FATAL => 'all';
use strict;
use Cwd qw(abs_path);
use POSIX qw(strftime);
use IO::File;
use Carp;

use Pod::Usage;
use Getopt::Long;


my $VERSION = '@BEETL_VERSION_FULL@';

my $programName = (File::Spec->splitpath(abs_path($0)))[2];
my $programPath = (File::Spec->splitpath(abs_path($0)))[1];
my $Version_text =
    "$programName $VERSION\n"
  . "Copyright (c) 2011 Illumina\n"
  . "This source file is covered by the \"Illumina Non-Commercial Use Software and\n"
  . "Source Code License Agreement\" and bound by the terms therein.\n";

my $usage =
    "Usage: $programName [options]\n"
  . "\t-i, --input=PATH             - {BWT}-end-pos file\n"
  . "\t-j, --input2=PATH            - ranges file\n"
  . "\t-k, --input3=PATH            - BWT prefix (temporary until we change endPos format)\n"
  . "\t-o, --output=PATH            - destination for output\n"

  . "\t--help                       - prints usage guide\n"
  . "\t--version                    - prints version information\n"

.<<'EXAMPLES_END';

EXAMPLES:
    (none)

EXAMPLES_END

my $help             = 'nohelp';
my $isVersion        = 0;
my %PARAMS           = ();

my $argvStr = join ' ', @ARGV;

$PARAMS{verbose} = 0;

$PARAMS{input} = "";
$PARAMS{input2} = "";
#$PARAMS{input3} = "";
$PARAMS{output} = "-";

my $result = GetOptions(
    "input|i=s"             => \$PARAMS{input},
    "input2|j=s"            => \$PARAMS{input2},
#    "input3|k=s"            => \$PARAMS{input3},
    "output|o=s"            => \$PARAMS{output},

    "version"               => \$isVersion,
    "help"                  => \$help
);

# display the version info
if ($isVersion) {
    print $Version_text;
    exit(0);
}

# display the help text when no output directory or other required options are given
if ( ( $result == 0 || !$PARAMS{input} || !$PARAMS{output} ) && 'nohelp' eq $help) {
    die "$usage";
}

die("ERROR: Unrecognized command-line argument(s): @ARGV")  if (0 < @ARGV);


# Check that we won't overwrite any existing file
(! -e "$PARAMS{output}") or die "$PARAMS{output} already exists. Aborting.";

# pile size
#my @pileStartPos;
#$pileStartPos[0] = 0;
#for (my $i=0; $i<6; $i++) {
#  my $filename = "$PARAMS{input3}-B0${i}";
#  my $pileSize = -s "$filename" or die "$filename doesn't exist";
#  $pileStartPos[$i+1] = $pileStartPos[$i] + $pileSize;
#  print "$pileStartPos[$i+1] = " . $pileStartPos[$i+1] . "\n";
#}

# Merging fragments.stats : sum of each int32
my $myInt32 = "";
my $myInt64 = "";
my $myInt8 = "";
open INF1, "<$PARAMS{input}" or die "Can't open $PARAMS{input}";
open OUTF, ">$PARAMS{output}" or die "Can't open $PARAMS{output} for writing";
binmode INF1;
binmode OUTF;
openRangeFile();

read (INF1, $myInt32, 4) or die "Can't read first 4 bytes in $PARAMS{input}";
my $entriesCount = unpack('L',$myInt32);
for (my $i=0; $i<$entriesCount; $i++) {
  read (INF1, $myInt32, 4) or die "Error reading 4 bytes in $PARAMS{input}";
  read (INF1, $myInt64, 8) or die "Error reading 8 bytes in $PARAMS{input}";
  read (INF1, $myInt8, 1) or die "Error reading 1 byte in $PARAMS{input}";
  my $seqN = unpack('L',$myInt32);
  my $posN = unpack('Q',$myInt64);
  my $pileN = unpack('C',$myInt8);
#  print OUTF "1 $posN\t$seqN\t$pileN\n";
#  $posN += $pileStartPos[$pileN];

#  print "posN=$posN \t seqN=$seqN\n";
  readAndProcessRangesUntil( $posN, $seqN, $pileN );

}
close INF1;
close OUTF;

print "EndPosFile successfully processed\n";



my $rangeFileAvailable;
my $currentRangePile;
my $currentRangePos;
my $currentRangeSize;

sub openRangeFile {
  if ($PARAMS{input2}) {
    open INF2, "<$PARAMS{input2}" or die "Can't open $PARAMS{input2}";
    $rangeFileAvailable = 1;
    $currentRangePile = 0;
    $currentRangePos = 0;
    $currentRangeSize = 0;
#    print OUTF "#seqN\n";
  }
  else {
    print OUTF "#posN\tseqN\n";
    $rangeFileAvailable = 0;
  }
}

sub readAndProcessRangesUntil {
  my ( $posN, $seqN, $pileN ) = @_;
  if ($rangeFileAvailable) {
    do {
      #print "$posN\t$seqN\t\tvs\t$currentRangePos\t$currentRangeSize\n";

      if ($currentRangePile < $pileN or ($currentRangePile == $pileN and $currentRangePos+$currentRangeSize <= $posN) ) {
        # Current range is before current position from endPos file => read next range
        #print " => read next range\n";
        my $line;
        if (not ($line = <INF2>)) {
          #print "No more ranges\n";
          exit 0;
        }
        if ($line =~ /READNUM\(1\) ([0-9]) ([0-9]*) ([0-9]*) ([0-9]*)/) {}
        $currentRangePile = $1;
        $currentRangePos = $2; # + $pileStartPos[$pile];
        $currentRangeSize = $3;
        my $dollarCount = $4;
        #print "  = $currentRangePile\t$currentRangePos\t$currentRangeSize\t$dollarCount\n";
      }
      else {
        # Current range ends after position from endPos file
        if ($currentRangePile == $pileN and $currentRangePos <= $posN) {
          # Current range includes position from endPos file
          #print " => GOOD!\t$pileN\t$seqN\t$posN\n";
          print OUTF "$seqN\n";
          return;
        }
        else {
          # Current range is after current position from endPos file => return from subfunction to read next endPos
          return;
        }
      }
    } while (1);
  }
  else {
    print OUTF "$posN\t$seqN\t$pileN\n";
  }
}
