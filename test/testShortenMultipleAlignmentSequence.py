'''
Created on Oct 27, 2019

@author: andrew davidson (aedavids@ucsc.edu)
'''
import logging
import msa2vcf.setupLogging as sl
from msa2vcf import ShortenMultipleAlignmentSequence as SMAS
import unittest


class TestShortenMultipleAlignmentSequence(unittest.TestCase):
    configFilePath = sl.setupLogging( default_path='logging.test.ini.json')
    logger = logging.getLogger(__name__)
    logger.info("using logging configuration file:{}".format(configFilePath))

    ################################################################################    
    def setUp(self):
        pass

    ################################################################################    
    def tearDown(self):
        pass

    ################################################################################ 
    @unittest.skip   
    def testTemplate(self):
        self.logger.info("BEGIN")
        self.logger.info("END\n")
        
    ################################################################################    
    def testSimple(self):
        self.logger.info("BEGIN")
        
        refSeq    = "ACGTACGT"
        MSARefSeq = "A|CGTA.CGT"
        sampleSeq = "-|-*--G---"
        self.logger.info("   refSeq:{}".format(refSeq))
        self.logger.info("MSARefSeq:{}".format(MSARefSeq))
        self.logger.info("sampleSeq:{}".format(sampleSeq))
        
        smas = SMAS.ShortenMultipleAlignmentSequence(refSeq, MSARefSeq, sampleSeq)
        n = 0
        self.logger.info("n = {}:\nsmas:\n{}".format(n, smas))
        
        e1 = "refSeq:GTACGT\n" + \
            "MSARefSeq:GTA.CGT\n" + \
            "sampleSeq:*--G---\n"
             
        n = 2
        smas.shorten(n)
        self.logger.info("n = {}:\nsmas:\n{}".format(n, smas))
        self.logger.info("e1:\n{}".format(e1))
        self.assertEqual(str(smas), e1)
        
        n = 4
        e2 = "refSeq:GT\n" + \
            "MSARefSeq:GT\n" + \
            "sampleSeq:--\n"        
        smas.shorten(n)
        self.logger.info("n = {}:\nsmas:\n{}".format(n, smas)) 
        self.assertEqual(str(smas), e2)
               
        self.logger.info("END\n")

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
