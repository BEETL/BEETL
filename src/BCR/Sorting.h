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

#ifndef SORTED_INCLUDED
#define SORTED_INCLUDED
		
//#include <stdlib.h>
#include <vector>
#include "Tools.h"

#if USE_ATTRIBUTE_PACKED == 1
struct sortElement {
  sortElement() {};
  sortElement( dataTypedimAlpha z, dataTypeNChar x,dataTypeNSeq y ) { pileN = z; posN = x; seqN = y; };
  ~sortElement() {};
  dataTypedimAlpha pileN;
  dataTypeNChar posN;
  dataTypeNSeq seqN;
}__attribute__ ((packed));
#else
struct sortElement {
  sortElement() {};
  sortElement( dataTypedimAlpha z, dataTypeNChar x,dataTypeNSeq y ) { pileN = z; posN = x; seqN = y; };
  ~sortElement() {};
  dataTypedimAlpha pileN;
  dataTypeNChar posN;
  dataTypeNSeq seqN;
};
#endif

void quickSort(std::vector< sortElement > &v);


#endif
