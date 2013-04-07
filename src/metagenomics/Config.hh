/**
 	** Copyright (c) 2011 Illumina, Inc.
 	**
 	**
 	** This software is covered by the "Illumina Non-Commercial Use Software
 	** and Source Code License Agreement" and any user of this software or
 	** source file is bound by the terms therein (see accompanying file
 	** Illumina_Non-Commercial_Use_Software_and_Source_Code_License_Agreement.pdf)
 	**
 	** This file is part of the BEETL software package.
 	**
 	**/

#ifndef DEFINED_CONFIG_HH
#define DEFINED_CONFIG_HH

/* Control flags - modify and recompile to change program behaviour */

//#define DEBUG 1
//#define USE_POSIX_FILE_OPTIMIZATIONS 1
#define USE_STANDARD_LEXICOGRAPHIC_ORDER 1
#define USE_ZLIB 1

// USE_4_BITS_PER_BASE: if set, convert the input sequences from ASCII into
// a 4-bits-per-base format and use that from then on
#define USE_4_BITS_PER_BASE 1


// USE_PREFIX_ONLY: if set, don't copy the whole of the sequences but
// only the prefix that's yet to be processed
// NB - USE_PREFIX_ONLY must only be used with USE_4_BITS_PER_BASE
#define USE_PREFIX_ONLY 1

// TRACK_SEQUENCE_NUMBER: if set, store an unsigned integer with each
// sequence that reports its originating position
//#define TRACK_SEQUENCE_NUMBER 1

// COMPRESS_BWT: if set, store the temporary BWT files in simple run-length
// encoded form
#define COMPRESS_BWT 1

// REPORT_COMPRESSION_RATIO: if set, output the compression achieved in the
// partial BWTs (only has an effect if COMPRESS_BWT is set
#define REPORT_COMPRESSION_RATIO 1

// REMOVE_TEMPORARY_FILES: if set, remove all files except BWT once creation
// is complete
#define REMOVE_TEMPORARY_FILES 1


/* Control parameters - modify and recompile to change program behaviour */

const int maxSeqSize( 1023 );
const int bwtBufferSize( 16384 ); // 1<<20=1048576

// Read this many sequences into RAM at once
static const unsigned int ReadBufferSize( 1024 );


#endif
