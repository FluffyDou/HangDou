#!/usr/bin/python
# Copyright (C) 2005-2008 Washington University in St Louis, Baylor College of Medicine.  All rights reserved
# Author:  Mike Marsh (michael.marsh@bcm.edu)
# Class:  Sheet
# Class Description: Class that models beta sheet, which act as containers for beta strand objects
#                    More info in: seq_model-doc.txt
#

from seq_model.Strand import Strand

class Sheet:
    def __init__(self, initialStrand=None):
        self.strandList = {}
        self.bonds      = [];
        
        if initialStrand is not None:
            self.strandList[initialStrand.strandNo]=initialStrand

    @classmethod
    def parsePDB(cls,line,chain,whichChainID):
        """
Given a line from a PDB and a chain, this creates a Sheet object and 
adds it to the chain.
        """
        # Parse Entry
        strandNo   = int(line[7:10].strip())
        sheetID    = line[11:14].strip()
        chainID    = line[21].strip()
        initSeqNum = int(line[22:26].strip())
        endSeqNum  = int(line[33:37].strip())
        label      = 'S' + sheetID

        # Check Chain
        if chainID != whichChainID:
            return

        # Create Sheet/Strands
        strand = Strand(chain, strandNo, label, initSeqNum, endSeqNum)
        if sheetID not in chain.sheets.keys():
            sheet = Sheet(strand)
            chain.addSheet(sheetID, sheet)
            chain.addStrand(strand, strandNo, sheetID)
        else:
            sheet = chain.sheets[sheetID]
            chain.addStrand(strand, strandNo, sheetID)
    
        # Calculate and Add Sheet Bonds
        sense = int(line[38:40].strip())
        if sense != 0:
            # Parse Entry
            curAtom     = line[41:45].strip()
            curChainId  = line[49]
            curResSeq   = int(line[50:54].strip())
            prevAtom    = line[56:60].strip()
            prevChainId = line[64]
            prevResSeq  = int(line[65:69].strip())
            
            # Given Bond
            sheet.bonds.append((curResSeq, prevResSeq))

            # Get Previous Strand
            currentStrand  = sheet.strandList[strandNo]
            previousStrand = sheet.strandList[strandNo-1]

            # Anti-Parallel
            if sense == -1:
                # Backward Iteration
                currentResidueIndexIterator  = curResSeq - 2
                previousResidueIndexIterator = prevResSeq + 2
                while currentStrand.startIndex <= currentResidueIndexIterator and currentResidueIndexIterator <= currentStrand.stopIndex and previousStrand.startIndex <= previousResidueIndexIterator and previousResidueIndexIterator <= previousStrand.stopIndex:
                    sheet.bonds.append((currentResidueIndexIterator, previousResidueIndexIterator))
                    
                    currentResidueIndexIterator  -= 2
                    previousResidueIndexIterator += 2

                # Forward Iteration
                currentResidueIndexIterator  = curResSeq + 2
                previousResidueIndexIterator = prevResSeq - 2
                while currentStrand.startIndex <= currentResidueIndexIterator and currentResidueIndexIterator <= currentStrand.stopIndex and previousStrand.startIndex <= previousResidueIndexIterator and previousResidueIndexIterator <= previousStrand.stopIndex:
                    sheet.bonds.append((currentResidueIndexIterator, previousResidueIndexIterator))
                    
                    currentResidueIndexIterator  += 2
                    previousResidueIndexIterator -= 2
            # Parallel
            elif sense == 1:
                for iteratorChange in (-2, 2):
                    currentResidueIndexIterator  = curResSeq
                    previousResidueIndexIterator = prevResSeq

                    if curAtom == 'O' and prevAtom == 'N':
                        if iteratorChange < 0:
                            firstStrandToIterate = "Previous"
                        elif iteratorChange > 0:
                            firstStrandToIterate = "Current"
                    elif curAtom == 'N' and prevAtom == 'O':
                        if iteratorChange < 0:
                            firstStrandToIterate = "Current"
                        elif iteratorChange > 0:
                            firstStrandToIterate = "Previous"

                    while True:
                        if firstStrandToIterate == "Previous":
                            previousResidueIndexIterator += iteratorChange
                            if previousStrand.startIndex <= previousResidueIndexIterator and previousResidueIndexIterator <= previousStrand.stopIndex:
                                sheet.bonds.append((currentResidueIndexIterator, previousResidueIndexIterator))
                            else:
                                break;
                            
                            currentResidueIndexIterator += iteratorChange
                            if currentStrand.startIndex <= currentResidueIndexIterator and currentResidueIndexIterator <= currentStrand.stopIndex:
                                sheet.bonds.append((currentResidueIndexIterator, previousResidueIndexIterator))
                            else:
                                break;
                        elif firstStrandToIterate == "Current":
                            currentResidueIndexIterator += iteratorChange
                            if currentStrand.startIndex <= currentResidueIndexIterator and currentResidueIndexIterator <= currentStrand.stopIndex:
                                sheet.bonds.append((currentResidueIndexIterator, previousResidueIndexIterator))
                            else:
                                break;
                           
                            previousResidueIndexIterator += iteratorChange
                            if previousStrand.startIndex <= previousResidueIndexIterator and previousResidueIndexIterator <= previousStrand.stopIndex:
                                sheet.bonds.append((currentResidueIndexIterator, previousResidueIndexIterator))
                            else:
                                break;

    def toPDB(self,sheetID):
        """
This generates a 'SHEET' line from this sheet for a PDB file.
        """
        s=''
        for strandNo in sorted(self.strandList.keys()):
            strand=self.strandList[strandNo]
            s=s+ strand.toPDB(sheetID,self)
        return s
