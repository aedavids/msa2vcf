'''
Created on Oct 13, 2019

@author: Andrew Davidson (aedavids@ucsc.edu)
'''
import logging

class VCFSample(object):
    '''
    ref: https://samtools.github.io/hts-specs/VCFv4.3.pdf
    
    TODO: add doc
    '''

    logger = logging.getLogger(__name__)    
    
    @staticmethod
    def header(sampleNamesList):
        """
        TODO: add doc
        """
        manditory = "CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT"
        samples = "\t".join(sampleNamesList)
        hdr = "#" + manditory + "\t" + samples
        return hdr

    def __init__(self, msfs, logicalVector):
        '''
        argument
            msfs: a MSFSample object
            
            logicalVector: a list of zeros and 1, use to identify which samples are match the variation
                described on msfs
        '''
        self.refChromName = msfs.refChromName 
        self.position = msfs.position 
        self.ident = msfs.ident 
        self.ref = msfs.ref 
        self.alt = msfs.alt 
        self.qual = msfs.qual 
        self.filterArg = msfs.filterArg 
        self.info = msfs.info
        self.format = msfs.format
        self.logicalVector = logicalVector
        
    def __repr__(self, *args, **kwargs):
        fmt = "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}"
        logicalVectorStr = "\t".join(str(i) for i in self.logicalVector)
        ret = fmt.format(self.refChromName, 
                        self.position, 
                        self.ident, 
                        self.ref, 
                        self.alt, 
                        self.qual, 
                        self.filterArg, 
                        self.info, 
                        self.format,
                        logicalVectorStr)
        
        return ret
        