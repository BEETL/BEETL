#!/usr/bin/python
import os
import sys
import struct
import csv


def getCommandlineOptions() :
    """
    get user options form the command-line and validate required arguments
    """

    from optparse import OptionParser

    parser = OptionParser()
    parser.add_option('-e', '--endpos', dest='endpos', type='string', help='path and name of the endpos index file')	
    parser.add_option('-i', '--fastqidx', dest='fastqidx', type='string', help='path and name of the fastq reads index file')
    parser.add_option('-f', '--infastq', dest='infastq', type='string', help='path and name of the input fastq file')
    parser.add_option('-o', '--outdir', dest='outdir', type='string', action='store', help='output directory path')
    parser.add_option('-p', '--paired-end read flag', dest='paired', type='string', action='store', default='1', help='is paired-end input fastq file flag (default 1)')
    (options, args) = parser.parse_args()

    # check for any errors or missing requirements from argument parsing:
    if len(args) > 0:
        parser.print_help()
        sys.exit(2)
    if (options.fastqidx is None) or (options.endpos is None) or (options.infastq is None) or (options.outdir is None) or (options.paired is None) :
        parser.print_help()
        sys.exit(2)

    options.endpos = os.path.abspath(options.endpos)
    options.fastqidx = os.path.abspath(options.fastqidx)
    options.infastq = os.path.abspath(options.infastq)
    return options

	
options=getCommandlineOptions()
options_outfastq1 = os.path.join(os.path.abspath(options.outdir),'lane1_read1.fastq')
fh_inidx = open(options.fastqidx, 'rt')
fh2_out = open(options_outfastq1, 'w+')
if options.paired == '1':
    options_outfastq2 = os.path.join(os.path.abspath(options.outdir),'lane1_read2.fastq')
    fh3_out = open(options_outfastq2, 'w+')
fh_endpos = open(options.endpos, 'rt')
line = fh_inidx.readline()

print "infastq ", options.infastq
print "fastqidx ", options.fastqidx
print "endpos file", options.endpos

line = []
tag = []
tag_sep = ' '
offset = 0
num_lines = struct.unpack('<i', fh_endpos.read(4))[0]
line1 = '1'
print "num_lines ", num_lines
fh1_in = open(options.infastq, 'rt')
offset = int((num_lines+1)/2)
print "offset ", offset

'''

'''

Reader = csv.reader(fh_inidx, delimiter=' ')
for row in Reader:
    if len(row) > 3:
        tag.extend([row[0]]*(len(row)-1))
        line.extend(row[3:len(row)])
line_sorted = [y for (y,x) in sorted(zip(line,tag))]	
tag_sorted  = [x for (y,x) in sorted(zip(line,tag))]	


'''
Print out read 1 
'''

if options.paired == '1':
    counter = 1
    counter_line_sorted = 0
    while  counter < offset:
        if  (counter  ==  line_sorted[counter_line_sorted]):
            fh2_out.write(fh1_in.readline()+ tag_sep + tag_sorted[counter_line_sorted])
            fh2_out.write(fh1_in.readline())
            fh2_out.write(fh1_in.readline())
            fh2_out.write(fh1_in.readline())
            counter_line_sorted += 1

        elif (counter  > line_sorted[counter_line_sorted]):
            while (counter  > line_sorted[counter_line_sorted]):
                counter_line_sorted += 1           
        else:
            line1 = fh1_in.readline()
            line1 = fh1_in.readline()
            line1 = fh1_in.readline()
            line1 = fh1_in.readline()
        counter += 1  
        if line == '' or line1 == '':
            print "weired index line ", line1, '\n ', counter
            break

    print "finished printing read 1 "

    line1 = fh1_in.readline()
    line1 = fh1_in.readline()
    line1 = fh1_in.readline()
    line1 = fh1_in.readline()
    counter_line_sorted = 0
    counter = 1
    while  counter < offset:

        if  (counter  ==  line_sorted[counter_line_sorted]):
            fh3_out.write(fh1_in.readline() + tag_sep + tag_sorted[counter_line_sorted])
            fh3_out.write(fh1_in.readline())
            fh3_out.write(fh1_in.readline())
            fh3_out.write(fh1_in.readline())
            counter_line_sorted += 1
        elif (counter  > line_sorted[counter_line_sorted]):
            while (counter  > line_sorted[counter_line_sorted]):	
                counter_line_sorted += 1                    
        else:
            line1 = fh1_in.readline()
            line1 = fh1_in.readline()
            line1 = fh1_in.readline()
            line1 = fh1_in.readline()
        counter += 1
        if line_sorted[counter_line_sorted] == '' or line1 == '':
            break
    fh1_in.close()

	
    counter_line_sorted = 0
    fh1_in = open(options.infastq, 'rt')
    counter_line_sorted += 1
    while (line_sorted[counter_line_sorted] < offset):
        counter_line_sorted += 1


    line1 = fh1_in.readline()
    line1 = fh1_in.readline()
    line1 = fh1_in.readline()
    line1 = fh1_in.readline()

    counter = 1
    while counter < offset:
        counter += 1
        line1 = fh1_in.readline()
        line1 = fh1_in.readline()
        line1 = fh1_in.readline()
        line1 = fh1_in.readline() 	
    
    counter = offset
    while counter < num_lines:
    
        if  (counter  ==  line_sorted[counter_line_sorted]):
            fh3_out.write(fh1_in.readline() + tag_sep + tag_sorted[counter_line_sorted])
            fh3_out.write(fh1_in.readline())
            fh3_out.write(fh1_in.readline())
            fh3_out.write(fh1_in.readline()) 
            counter_line_sorted += 1
        elif (counter  > line_sorted[counter_line_sorted]):
            while (counter  > line_sorted[counter_line_sorted]):	
                counter_line_sorted += 1            
        else:
            line1 = fh1_in.readline()
            line1 = fh1_in.readline()
            line1 = fh1_in.readline()
            line1 = fh1_in.readline()
        counter += 1
        if line_sorted[counter_line_sorted] == '' or line1 == '':
            break

else:
    counter = 0
    counter_line_sorted = 0
    while  counter < num_lines and counter_line_sorted < len(line_sorted):
        if  (counter  ==  int(line_sorted[counter_line_sorted])):
            fh2_out.write(fh1_in.readline().rstrip() + tag_sep + tag_sorted[counter_line_sorted] + '\n')
            fh2_out.write(fh1_in.readline())
            fh2_out.write(fh1_in.readline())
            fh2_out.write(fh1_in.readline())
            counter_line_sorted += 1
 		
        elif (counter  > int(line_sorted[counter_line_sorted])):
            while (counter  > int(line_sorted[counter_line_sorted])):
                counter_line_sorted += 1           
        else:
            line1 = fh1_in.readline()
            line1 = fh1_in.readline()
            line1 = fh1_in.readline()
            line1 = fh1_in.readline()
            counter += 1  
        if line == '' or line1 == '':
            print "weired index line ", line1, '\n ', counter
            break


    print "finished printing read 1 "
   
fh1_in.close()
fh_inidx.close()
