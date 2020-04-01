import os
import sys
import subprocess as s
from subprocess import Popen
import argparse
from typing import List, Tuple
import pandas as pd
import re
from itertools import groupby


class MRRI():

    def __init__(self, args, idx1Pos0, idx2Pos0):
        '''
            Constructor for the MRRI class with default value from args
        '''
        self.regex = "[ACGTUacgtu]+"
        self.args = args
        self.querySeq = self.parseFasta2dict( args.query, 'query' ) # dictionary of queryID 2 query sequence
        self.targetSeq = self.parseFasta2dict( args.target, 'target' ) # dictionary of targetID 2 target sequence
        self.b1 = None
        self.b2 = None
        self.b3 = None
        self.idx1Pos0 = idx1Pos0
        self.idx2Pos0 = idx2Pos0
    
    
    def csv2dict( self, intarnaCsvOut ):
        if (intarnaCsvOut == None):
            return {}
        keys = intarnaCsvOut.split("\n")[0].split(";")
        values = intarnaCsvOut.split("\n")[1].split(";")
        op = zip(keys, values)
        return dict(op)

    def read_fasta_file(self, fastaname):
        """
        given a fasta file. yield tuples of id, sequence
        """
        fh = open(fastaname)
        fastaiter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))
        for idx in fastaiter:
            # drop the ">"
            idx = idx.next()[1:].strip()
            # join all sequence lines to one.
            seq = "".join(s.strip() for s in fastaiter.next())
            yield idx, seq

      

    def parseFasta2dict(self, sequenceInput, seqname ):
       
        if os.path.isfile(sequenceInput) == True:
            fastaReturn = self.read_fasta_file(sequenceInput)
            return fastaReturn
        else:
            matches = re.finditer(self.regex, sequenceInput, re.MULTILINE | re.IGNORECASE)
            tempDict = {}
            tempDict[seqname] = ''
            for match in matches:
                if match.group().strip("") != "":
                    tempDict[seqname] += str(match.group())
            return tempDict
   
    def runIntaRNA(self, B1= None):
        '''
            Run IntaRNA on given query and target
            params have value of start and end numbers
        
        '''
        complete = str(self.args.intarnaBin) + ' -q ' + self.querySeq[B1['id2']] + ' -t ' + self.targetSeq[B1['id1']] +' --outMode=C -n 1 '#--energyNoDangles'
        # add parameterFile to call if given
        if self.args.parameterFile :
            complete += " --parameterFile="+self.args.parameterFile
        # set require CSV columns
        complete += " --outCsvCols=id1,start1,end1,id2,start2,end2,subseqDP,hybridDP,E,E_hybrid,ED1,ED2"
        # TODO (somewhen) parse parameterFile for outCsvCols and add user-requested csv-col ids not already within the list
        if B1.keys().__contains__('start1') and B1.keys().__contains__('start2'):
            complete += ' --tAccConstr="b:'+B1['start1']+'-'+B1['end1']+'" --qAccConstr="b:'+B1['start2']+'-'+B1['end2']+'" '
        return self.csv2dict(self.runCmdLine(complete)[0].replace("query",B1['id2']).replace("target",B1['id1']))
    


    def get0basedIndex( self, idx, origIdxPos0 ):
        #print('origIdxPos0:', origIdxPos0)
        idxNew = [x - origIdxPos0 for x in map(int, idx.copy())]
        if int(origIdxPos0) < 0:
            # shift positive range by -1 since 0 is skipped  
            for i in range(len(idxNew)-1):
                if int(idx[i]) > 0:
                    idxNew[i] =  int(idxNew[i])-1
        return idxNew
                

    def getEDunconstraint(self, B1 ):
        complete = str(self.args.intarnaBin) + ' -q ' + self.querySeq[B1['id2']] + ' -t ' + self.targetSeq[B1['id1']] +' --out=/dev/null -n 0 '#--energyNoDangles'
        # add parameterFile to call if given
        if self.args.parameterFile :
            complete += " --parameterFile="+self.args.parameterFile
        # set require ED output
        complete += " --out=tAcc:STDOUT --out=qAcc:STDERR"

        # first convert indices into 0-based indexing
        startEnd1 = self.get0basedIndex(idx=[B1['start1'], B1['end1']], origIdxPos0=int(self.idx1Pos0))
        startEnd2 = self.get0basedIndex(idx=[B1['start2'], B1['end2']], origIdxPos0=int(self.idx2Pos0))
        
        # confine output to region of interest
        l1 = startEnd1[1]-startEnd1[0]+1
        l2 = startEnd2[1]-startEnd2[0]+1

        complete += " --tIntLenMax=" + str(l1) + " --qIntLenMax=" +str(l2)
        # run intarna
        outputED = self.runCmdLine( complete )
        # parse ED1 from stdout == outputED[0]      
        ED1 = round(float(outputED[0].split('\n')[ startEnd1[1]+2 ].split()[l1]),2)
            # parse ED2 from stderr == outputED[1]
        ED2 = round(float(outputED[1].split('\n')[ startEnd2[1]+2 ].split()[l2]),2)
            # return combined constraints
        return [ED1,ED2]
        

    def runCmdLine(self,completeCall):
        #print(completeCall)
        ps = s.Popen(str(completeCall),  stdout = s.PIPE, stderr=s.PIPE, universal_newlines=True)
        (stdout, stderr) = ps.communicate()
        if ps.returncode:  # If IntaRNA exits with a returncode != 0, skip this iteration
            sys.stderr.write("IntaRNA failed ({}) for call {}\n".format(stderr,completeCall))
            exit(ps.returncode)
        #print('output of ps command:', stdout)
        return [stdout,stderr]

