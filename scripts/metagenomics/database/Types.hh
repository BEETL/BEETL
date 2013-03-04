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
 	** Citation: Markus J. Bauer, Anthony J. Cox and Giovanna Rosone
 	** Lightweight BWT Construction for Very Large String Collections.
 	** Proceedings of CPM 2011, pp.219-231
 	**
 	**/
 	
 	#ifndef DEFINED_TYPES_HH
 	#define DEFINED_TYPES_HH
 	
 	
 	/* Standard data types */
 	
 	typedef unsigned int uint;
 	typedef unsigned char uchar;
 	
 	// below limits to 4 billion reads max - change to unsigned long for more
 	typedef unsigned int SequenceNumberType;
 	
 	// Should work for BWT of up to 2^64 characters in size
 	typedef unsigned long long LetterCountType;
 	
 	const LetterCountType maxLetterCountType(static_cast<LetterCountType>(-1));
 	
 	#endif
