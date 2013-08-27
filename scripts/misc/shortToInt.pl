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

shortToInt.pl

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
  . "This source file is covered by the \"Illumina Non-Commercial Use Software and\"\n"
  . "Source Code License Agreement and bound by the terms therein.\n";

my $usage =
    "Usage: $programName [options]\n"
  . "\t-i, --input=PATH             - binary input file containing unsigned short numbers\n"
  . "\t-o, --output=PATH            - destination for output with unsigned int numbers\n"

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

my $myInt16 = "";
my $myInt32 = "";
open INF1, "<$PARAMS{input}" or die "Can't open $PARAMS{input}";
open OUTF, ">$PARAMS{output}" or die "Can't open $PARAMS{output} for writing";
binmode INF1;
binmode OUTF;

while (read (INF1, $myInt16, 2)) {
  my $number = unpack('S',$myInt16);
  $myInt32 = pack('L',$number) or die "Odd number of bytes in input file. Ignoring last byte.";
  print OUTF $myInt32;
}
close INF1;
close OUTF;

print "Conversion complete\n";
