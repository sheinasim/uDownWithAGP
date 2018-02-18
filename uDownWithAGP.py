#!/usr/bin/env python

import argparse
import subprocess as sp
import sys
import re
import numpy as np

def readFasta(fastaToParse):
    name, seq = None, []
    for line in fastaToParse:
        line = line.rstrip()
        if line.startswith(">"):
            if name: yield name, ''.join(seq)
            name, seq = line, []
        else:
            seq.append(line)
    if name: yield (name, ''.join(seq))

def fasta2agp(inFasta, minGapLength, outFile, linkageEvidence):
    with open(inFasta) as fastaFile:
        outAGP = open(outFile, "w")
        for name, seq in readFasta(fastaFile):
	    numberInSequence = 1
            name = name.strip(">").split()
            listOfHeaders = []
	    listOfScaffolds = []
            NString = "N"*minGapLength
            listOfHeaders.append(name[0])
	    listOfScaffolds.append(seq)
            listOfContigs = []
            listOfGaps = []
            contigs = seq.split(NString)
            contigs = filter(None, contigs)
            strippedContigs = []
            for contig in contigs:
                contig = contig.strip("N")
                strippedContigs.append(contig)
            listOfContigs.append(strippedContigs)
            gaps = re.split("A|G|T|C", seq)
            gaps = filter(None, gaps)
            newGaps = [gap for gap in gaps if len(gap) > 10]
            listOfGaps.append(newGaps)
            for i, x in enumerate(listOfHeaders):
		scaffoldStart = "1"
		scaffoldEnd = str(len(listOfScaffolds[i]))
                for j in range(len(listOfContigs)):
                    if len(listOfContigs[j]) == 1:
                        listA = [] # ["object", "object_beg", "object_end", "part_number", "component_type", "component_ID or gap_length", "component_beg or gap_type", "component_end or linkage", "orientation or linkage evidence"]
                        listA.append(listOfHeaders[i]) # scaffold
                        contigStartPos = str(1)
                        listA.append(scaffoldStart) # scaffold_beg
                        contigEndPos = str(len(listOfContigs[j][0]))
                        listA.append(scaffoldEnd)
                        listA.append(str(numberInSequence)) # part_number
                        listA.append("W") # component_type
                        listA.append("scaff"+str(listOfHeaders[i])+".contig"+str(j)) # component_ID or gap_length
                        listA.append(str(contigStartPos)) # component_beg or gap_type
                        listA.append(str(contigEndPos)) # component_end or linkage
                        listA.append("+") # orientation or linkage_evidence
                        listA.append("\n")
                        numberInSequence += 1
                        outAGP.write('\t'.join(map(str,listA)))
                    else:
                        for k in range(len(listOfContigs[j])):
                            listB = []
                            listB.append(listOfHeaders[i])
                            contigStartPos = seq.find(listOfContigs[j][k])+1
                            listB.append(scaffoldStart)
                            contigEndPos = seq.find(listOfContigs[j][k])+len(listOfContigs[j][k])
                            listB.append(scaffoldEnd)
                            listB.append(str(numberInSequence))
                            listB.append("W")
                            listB.append("scaff"+str(listOfHeaders[i])+".contig"+str(k)) # component_ID or gap_length
                            listB.append(str(contigStartPos)) # component_beg or gap_type
                            listB.append(str(contigEndPos)) # component_end or linkage
                            listB.append("+") # orientation or linkage_evidence
                            listB.append("\n")
                            numberInSequence += 1
                            outAGP.write('\t'.join(map(str,listB)))
                            if k < len(listOfGaps[j]):
                                if len(listOfGaps[j][k]) == 100:
                                    listC = []
                                    listC.append(listOfHeaders[i])
                                    listC.append(scaffoldStart)
                                    listC.append(scaffoldEnd)
                                    listC.append(str(numberInSequence))
                                    listC.append("U")
                                    listC.append("100")
                                    listC.append("scaffold")
                                    listC.append("yes")
                                    listC.append(linkageEvidence) # orientation or linkage_evidence
                                    listC.append("\n")
                                    outAGP.write('\t'.join(map(str,listC)))
                                else:
                                    listD = []
                                    listD.append(listOfHeaders[i])
                                    listD.append(scaffoldStart)
                                    listD.append(scaffoldEnd)
                                    listD.append(str(numberInSequence))
                                    listD.append("N")
                                    listD.append(str(len(listOfGaps[j][k])))
                                    listD.append("scaffold")
                                    listD.append("yes")
                                    listD.append(linkageEvidence)
                                    listD.append("\n")
                                    outAGP.write('\t'.join(map(str,listD)))
                                numberInSequence += 1
                            else:
                                pass
        outAGP.close()

parser = argparse.ArgumentParser(description="Write NCBI .agp file for scaffold .fasta.")
parser.add_argument("-?", "--question", action="store_true")
parser.add_argument("-f", "--inFasta", type=argparse.FileType('r'), help=".fasta file containing all your scaffolds.")
parser.add_argument("-minGap", "--minimumGapLength", help="Number of Ns denoting gap of unknown length.", type=int, default=10)
parser.add_argument("-out", "--outAGP", help=".agp file containing all your reads and gaps.", type=str, default="out.agp")
parser.add_argument("-LE", "--linkageEvidence", help="Evidence used to scaffold assembly.", type=str, default="unspecified")
args = parser.parse_args()

if args.question and not args.inFasta:
    print "Yeah you know me. \nPlease supply -f <.fasta> "
elif args.question and args.inFasta:
    print "Yeah you know me. "
    fasta2agp(args.inFasta.name, args.minimumGapLength, args.outAGP, args.linkageEvidence)
elif not args.question and args.inFasta:
    fasta2agp(args.inFasta.name, args.minimumGapLength, args.outAGP, args.linkageEvidence)
else:
    print "Please supply -f <.fasta> "
