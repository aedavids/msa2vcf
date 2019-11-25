'''
Created on Oct 11, 2019

@author: Andrew Davidson (aedavids@ucsc.edu)
'''
from enum import Enum
import logging

class VariantType(Enum):
    """
    FIXME: add doc
    """
    unknown = 0
    polymorphism = 1
    insertion = 2
    deletion = 3
    
################################################################################                
class MSASample(object):
    """
    parsed representation of a sample from a MSF file plus information needed to convert to
    a VCF 
    """
    
    logger = logging.getLogger(__name__)    
    
    ################################################################################                
    def __init__(self, refChromName, sampleName, position, ident, ref, alt, qual, filterArg, info, formatArg, variantType):
        """
        constructor
        """
        self.refChromName = refChromName 
        self.sampleName = sampleName 
        self.position = position 
        self.ident = ident 
        self.ref = ref 
        self.alt = alt 
        self.qual = qual 
        self.filterArg = filterArg 
        self.info = info
        self.format = formatArg
        self.variantType = variantType
                
#     def __str__(self, *args, **kwargs):
#         return object.__str__(self, *args, **kwargs)
    
#     ################################################################################                
#     def __repr__(self, *args, **kwargs):
#         return object.__repr__(self, *args, **kwargs)
     
    ################################################################################                
#     def __str__(self):
    def __repr__(self, *args, **kwargs):
#     def __str__(self, *args, **kwargs):
        """
        create a string representation of MSFSample, to make testing easier
        """
        fmt = "CHROM:{} sample:{} POS:{} ID:{} REF:{} ALT:{} QUAL{} FILTER:{} INFO:{} FORMAT:{} {}"
        ret = fmt.format(self.refChromName, self.sampleName, self.position, self.ident, 
                        self.ref, self.alt, self.qual, self.filterArg, self.info, self.format, self.variantType)
        
        return ret  