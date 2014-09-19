File formats used in BEETL
==========================

BWT
---

## RLE_v3 (Run Length Encoded) BWT -- current default format

This format allows different alphabet letters to be associated to different numbers of bytecodes in order to represent longer run lengths for A,C,G,T and shorter ones for $,N.  
It also includes a "continuation symbol" '+', which changes the usual additive process of consecutive symbols of the same alphabet letter, into a multiplicative process.

- Bytes 0-5: Magic header "BWT\13\10\26"  
    - "BWT" identifier
    - "\13\10" = "\r\n" sequence to check for invalid dos/unix format conversions  
    - "\26" = Ctrl-Z, makes some text viewers stop here and is non-ASCII to avoid confusion with text files
- Bytes 6-7: Format version number, on 2 bytes to help identify endianness problems.  
    - Value: 3 as "\3\0"
- Bytes 8-?: Conversion table, to match bytecodes to {letter, run length}
    - Made of ranges following the format:
        - base: 1 char
        - range length: 1 byte
        - first run length: 2 bytes
    - Must contain enough ranges to cover 256 entries  
    - The current BWT writer uses these ranges:  
                { 'A', 58, 1 } // bytecode 0 = {letter='A',runLength=1}, bytecode 1 = {'A',2} ... until runLength 58  
                { 'C', 58, 1 } // same 58 run lengths for C, G and T  
                { 'G', 58, 1 }  
                { 'T', 58, 1 }  
                { 'N', 4, 1 }  // 'N' and '$' only get run lengths from 1 to 4  
                { '$', 4, 1 }  
                { '+', 16, 0 } // the continuation symbol can take values from 0 to 15
- Following bytes: data


Example:  
When a continuation symbol is encountered, the run length is calculated as in:  
  {'A',5},{'+',3},{+,11} => runLength = 5 + 3*58 + 11*58*16 (little-endian in base 58 for first digit and base 16 for other digits, the base values being as defined in the bytecode conversion table)


## ASCII BWT

Easy: 1 ASCII char per BWT letter


## RLE BWT (for backward compatibility)

Each byte is made of:
- 4 low bits to encode the letter ($=0, A=1, C=2, G=3, N=4, T=5 by default)
- 4 high bits to encode the run length (number of times the character appears)


## RLE53 BWT

Each byte is made of:
- 5 high bits to encode the run length (number of times the character appears)
- 3 low bits to encode the letter ($=0, A=1, C=2, G=3, N=4, T=5 by default)



## end-pos file

This file maps the BWT dollar signs to their sequence numbers.

Header:
- `sequenceGroupCount`: 4 bytes (type SequenceNumber): number of sequences (counting main read, paired read and reverse-complemented entries together as one)
- `sequenceCountInGroup`: 1 byte (type uint8\_t): 2 if the BWT contains paired reads, otherwise 1 (future extension: used for more than 2 reads per sequence group)
- `hasRevComp`: 1 byte (type bool): true if the BWT includes reverse-complemented reads

The total number of sequences in the BWT (number of '$' signs) is `dollarSignCount = sequenceGroupCount * sequenceCountInGroup * (hasRevComp+1)`.

Body, for each of these `dollarSignCount` entries, containing information about the k^th dollar sign present in the BWT:
- `sequenceGroupNum`: 4 bytes (type SequenceNumber): 0-based sequence number in the original input file (e.g. fasta) modulo sequenceGroupCount
- `positionInGroup`: 1 byte (type uint8\_t): 0-based position in group (e.g. if sequenceCountInGroup==2 and hasRevComp==true: 0=main read, 1=paired read, 2=reverse-complemented main read, 3=rev-comp paired read)

The sequence number in the overall input file is `sequenceNum = sequenceGroupNum + positionInGroup * sequenceGroupCount`.


## index file (from beetl-index)




Meta-BEETL
----------

Key identifiers used across all the files:

1. GenBank id
2. metaBeetl genome num (with different numbers for _rev version of each genome)
3. metaBeetl file num
4. taxId



## headerFile.csv

//GenBank id (1), chr name, metaBeetl genome num (2)
metaBeetl genome num (2 without_rev), chr name


## filecounter.csv

metaBeetl file num (3), metaBeetl genome num (2)


## gi_taxid_nucl.dmp

GenBank id (1), taxId (4)


## nodes.dmp

taxId (4) | parent Tax Id | taxLevel | ...


## names.dmp

taxId (4) | name


## ncbiFileNumToTaxTree

metaBeetl file num (3), Rank 0 taxId, Rank 1 taxId, ..., taxId (4) [, then filled with taxId=0 until Rank 11]
