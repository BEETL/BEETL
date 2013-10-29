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

readIntervalFile.pl

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
#    "input2|j=s"            => \$PARAMS{input2},
#    "input3|k=s"            => \$PARAMS{input3},
#    "output|o=s"            => \$PARAMS{output},

    "version"               => \$isVersion,
    "help"                  => \$help
);

# display the version info
if ($isVersion) {
    print $Version_text;
    exit(0);
}

# display the help text when no output directory or other required options are given
if ( ( $result == 0 || !$PARAMS{input} ) && 'nohelp' eq $help) {
    die "$usage";
}

die("ERROR: Unrecognized command-line argument(s): @ARGV")  if (0 < @ARGV);


# Check that we won't overwrite any existing file
#(! -e "$PARAMS{output}") or die "$PARAMS{output} already exists. Aborting.";

# Merging fragments.stats : sum of each int32
my $myInt32 = "";
my $myInt64 = "";
my $myInt8 = "";
open INF1, "<$PARAMS{input}" or die "Can't open $PARAMS{input}";
#open OUTF, ">$PARAMS{output}" or die "Can't open $PARAMS{output} for writing";
binmode INF1;
#binmode OUTF;
#openRangeFile();

my $pos = 0;
while (1) {
  my $num1 = readCompressedNum();
  my $offset = int($num1/4);
  print "offset=${offset}\n";
  $pos += $offset;
  print "pos=${pos}\n";

  my $len = readCompressedNum();
  print "len=${len}\n";
}

close INF1;
#close OUTF;

print "IntervalFile successfully processed\n";



sub readCompressedNum {
  read (INF1, $myInt8, 1) or die "Error reading 1 byte in $PARAMS{input}";
  my $val = unpack('C',${myInt8});
#  print "val=${val}\n";
  my $count = $val & 15;
  $val = $val >> 4;
#  print "count=${count}\n";
#  print "val=${val}\n";

#  my $num2 = 0;
  my $shift = 4;
  for (my $i=0; $i<$count; $i++) {
    read (INF1, $myInt8, 1) or die "Error reading 1 byte in $PARAMS{input}";
    my $byte = unpack('C',${myInt8});
    $val = ($byte << $shift) | $val;
    $shift += 8;
  }

  return $val;
}
