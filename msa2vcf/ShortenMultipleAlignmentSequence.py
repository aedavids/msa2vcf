'''
Created on Oct 27, 2019

@author: Andrew Davidson (aedavids@ucsc.edu)
'''

import argparse
import logging
import os
#from setupLogging import setupLogging
import sys


class ShortenMultipleAlignmentSequence:
    logger = logging.getLogger(__name__)
    
    ################################################################################            
    def __init__(self,refSeq, MSARefSeq, sampleSeq):
        '''
        arguments:
            refSequece: the reference sequence
                        example : AAAAA
            MSARefSeq: the reference sequence in a multiple sequence alignment file
                        example:  A..AGGC
            sampleSeq: the sample sequence in a multiple sequence alignment file
                        example: *--G**
        '''
        self.refSeq = refSeq
        self.MSARefSeq = MSARefSeq
        self.sampleSeq = sampleSeq
    ################################################################################                
    def __repr__(self, *args, **kwargs):  
        fmt = "refSeq:{}\nMSARefSeq:{}\nsampleSeq:{}\n"  
        ret = fmt.format(self.refSeq, self.MSARefSeq, self.sampleSeq) 
        return ret   
        
    ################################################################################            
    def shorten(self, n):
        """
        shortens refSequence, MSARefseq, and sampleSeq such that they remain aligned
        
        for example:
            input:
                n = 2
                refSequece = ACGTA
                MSARefSeq =  A|CGT.A
                sampleSeq =  -|---G-
                
            output:
                refSequece = GTA
                MSARefSeq =  GT.A
                sampleSeq =  --G-            
        """
        if n == 0:
            return
        
        self.refSeq = self.refSeq[n:]
        position = 1
        i = 0
        metaSymbols = {'.', '|'}
        while position <= n:
            if self.MSARefSeq[i] not in  metaSymbols:
                position += 1
            i += 1
            
        self.MSARefSeq = self.MSARefSeq[i:]
        self.sampleSeq = self.sampleSeq[i:]
        
################################################################################                
def main():
    logConfigFilePath = setupLogging("logging.ini.json")
    logger = logging.getLogger(__name__)
    logger.info("using logging config file:{}".format(logConfigFilePath))
        
    description = "trims bases from reference sequence while preserving the alignment" \
                + " See logging.ini.json to change log level"
                        
    usage = "n inputFile outputFile [options]"
    parser = argparse.ArgumentParser(description=description, usage=usage)

    helpStr = "{} {} {} {}".format(
        "This file has 3 lines.",
        "The first line is the reference sequence.",
        "The second line is the reference sequence in multiple sequence alignment format.",
        "The third line is the sample sequence in multiple sequence alignment format."
        )
    parser.add_argument("n", type=int,
                        help="The number of nucleotides to trim from the front of the reference sequenc")    
    parser.add_argument("inputFile", type=str,
                        help=helpStr)
    parser.add_argument("outputFile", type=str,
                        help="The trimed version of the input file")
    
    parser.add_argument('--version', action='version', version='%(prog)s 1.0')
    
    
    # if no arguments print help
    if len(sys.argv)==1:
        parser.print_help(sys.stderr)
        sys.exit(1)
        
    options = parser.parse_args()  
    
    n = options.n
    # make sure input file exists and output does not  
    inputFilePath = options.inputFile
    if not os.path.isfile(inputFilePath) :
        logger.error("not found inputFilePath: {}".format(inputFilePath))
        return 1
    
    # make sure output file does not exist
    outputFilePath = options.outputFile
    if os.path.isfile(outputFilePath) :
        logger.error("outputFilePath can not exist: {}".format(outputFilePath))
        return 1  
    
    with open(inputFilePath) as fp:
        refSeq, MSARefSeq, sampleSeq = fp 
        
    smas = ShortenMultipleAlignmentSequence(refSeq, MSARefSeq, sampleSeq)
    smas.shorten(n)
    # "a+" create if it does not exist else append
    with open(outputFilePath, "a+") as fp:
        fp.write(smas.refSeq)
        fp.write(smas.MSARefSeq)
        fp.write(smas.sampleSeq)
        
    return 0

if __name__ == '__main__':
    main()
