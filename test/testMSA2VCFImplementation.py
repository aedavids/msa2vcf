'''
Created on Oct 10, 2019

@author: andrewdavidson (aedavids@ucsc.edu)
'''
from msa2vcf import MSA2VCFImplementation as MSAVCF
import msa2vcf.setupLogging as sl
import logging
import os
import pprint
import unittest


class TestMSA2VCFImplementation(unittest.TestCase):
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
    @unittest.skip 
    def testLogging(self):
        print("AEDWIP pwd:{}".format(os.getcwd()))
        self.logger.info("pwd:{}".format(os.getcwd()))
        self.logger.info("BEGIN")
        self.logger.info("END\n")
        

    ################################################################################   
    def testParseMSASample_GoodSNP(self):
        self.logger.info("BEGIN")
        
        # data/A_gen.txt
        # snp at beginning
        ref  = "A*01:01:01:01    CAGGAGCAGA GGGGTCAGGG ATTGGGGAGT CCCAGCCTTG "
        seq1 = "A*01:01:01:02    T--------- ---------- ---------- ---------- "
        # CHROM:A*01:01:01:01 sample:A*01:01:01:02 POS:1 ID:. REF:C ALT:T QUAL. FILTER:PASS INFO: FORMAT:GT VariantType.polymorphism
                
        # snp in middle
        seq2 = "A*01:01:01:03    -----C---- ---------- ---------- ---------- "
        # result2:A*01:01:01:03    6    .    G    C    .    PASS    .
        # "CHROM:A*01:01:01:01 sample:A*01:01:01:03 POS:6 ID:. REF:G ALT:C QUAL. FILTER:PASS INFO:. FORMAT:GT VariantType.polymorphism"
        
        # snp at end of block
        seq3 = "A*01:01:01:04    ---------C ---------- ---------- ---------- "
        # result3:A*01:01:01:04    10    .    A    C    .    PASS    .
        # CHROM:A*01:01:01:01 sample:A*01:01:01:04 POS:10 ID:. REF:A ALT:C QUAL. FILTER:PASS INFO:. FORMAT:GT VariantType.polymorphism
        
        # multiple snps
        seq4 = "A*01:01:01:05    ---------C ---------- -AA------- ---------- "
        # A*01:01:01:05    10    .    A    C    .    PASS    .
        # "CHROM:A*01:01:01:01 sample:A*01:01:01:05 POS:10 ID:. REF:A ALT:C QUAL. FILTER:PASS INFO:. FORMAT:GT VariantType.polymorphism",
        # A*01:01:01:05    22    .    T    A    .    PASS
        #"CHROM:A*01:01:01:01 sample:A*01:01:01:05 POS:22 ID:. REF:T ALT:A QUAL. FILTER:PASS INFO:. FORMAT:GT VariantType.polymorphism",
        # A*01:01:01:05    23    .    T    A    .    PASS    .
        # "CHROM:A*01:01:01:01 sample:A*01:01:01:05 POS:23 ID:. REF:T ALT:A QUAL. FILTER:PASS INFO:. FORMAT:GT VariantType.polymorphism"        
        
        msa2vcf = MSAVCF.MSA2VCFImplementation()
        
        offset = 0
        result1List = msa2vcf._parseMSASample(offset, ref, seq1)
        self.logger.info("result1[0]:\n{}\t".format(result1List[0]))
        r1Expected = "CHROM:A*01:01:01:01 sample:A*01:01:01:02 POS:1 ID:. REF:C ALT:T QUAL. FILTER:PASS INFO:. FORMAT:GT VariantType.polymorphism"
        self.assertEqual(str(result1List[0]), r1Expected)
        
         
        result2List = msa2vcf._parseMSASample(offset, ref, seq2)
        self.logger.info("result2[0]:\n{}".format(result2List[0]))
        r2Expected = "CHROM:A*01:01:01:01 sample:A*01:01:01:03 POS:6 ID:. REF:G ALT:C QUAL. FILTER:PASS INFO:. FORMAT:GT VariantType.polymorphism"
        self.assertEqual(str(result2List[0]), r2Expected)
        
        result3List = msa2vcf._parseMSASample(offset, ref, seq3)
        self.logger.info("result3[0]:\n{}".format(result3List[0]))
        r3Expected = "CHROM:A*01:01:01:01 sample:A*01:01:01:04 POS:10 ID:. REF:A ALT:C QUAL. FILTER:PASS INFO:. FORMAT:GT VariantType.polymorphism"
        self.assertEqual(str(result3List[0]), r3Expected)        
        
        result4List = msa2vcf._parseMSASample(offset, ref, seq4)
        r4Expected = [
            "CHROM:A*01:01:01:01 sample:A*01:01:01:05 POS:10 ID:. REF:A ALT:C QUAL. FILTER:PASS INFO:. FORMAT:GT VariantType.polymorphism",
            "CHROM:A*01:01:01:01 sample:A*01:01:01:05 POS:22 ID:. REF:T ALT:A QUAL. FILTER:PASS INFO:. FORMAT:GT VariantType.polymorphism",
            "CHROM:A*01:01:01:01 sample:A*01:01:01:05 POS:23 ID:. REF:T ALT:A QUAL. FILTER:PASS INFO:. FORMAT:GT VariantType.polymorphism"
        ]
        for i in range(len(result4List)):
            msas = result4List[i]
            self.logger.info("msas{}:\n{}".format(i, msas))
            self.assertEqual(str(msas), r4Expected[i])
               
        self.logger.info("END\n")
        
        
    ################################################################################ 
    def testParseMSASample_DeleteSNP(self):
        self.logger.info("BEGIN")
        
        # data/A_gen.txt
        # SNP delete at beginning
        # true                   CAG
        # pos                    123
        ref  = "A*01:01:01:01    CAGGAGCAGA GGGGTCAGGG ATTGGGGAGT CCCAGCCTTG "
        seq1 = "A*01:01:01:02    .--------- ---------- ---------- ---------- "
#       result1:A*01:01:01:02    1    .    C    .    .    PASS    .
#       CHROM:A*01:01:01:01 sample:A*01:01:01:02 POS:1 ID:. REF:C ALT:. QUAL. FILTER:PASS INFO:. FORMAT:GT VariantType.deletion

        # SNP delete at middle of block
        #true                    CAGGAGCAGA
        #pos                     1234567890
        seq2 = "A*01:01:01:03    -----.---- ---------- ---------- ---------- "
#       result2:A*01:01:01:03    6    .    AG    A    .    PASS    .
#       CHROM:A*01:01:01:01 sample:A*01:01:01:03 POS:5 ID:. REF:AG ALT:A QUAL. FILTER:PASS INFO:. FORMAT:GT VariantType.deletion

        # SNP delete at end of block
        #true                    CAGGAGCAGA
        #pos                     1234567890
        seq3 = "A*01:01:01:04    ---------. ---------- ---------- ---------- "
#       result3:A*01:01:01:04    10    .    GA    G    .    PASS    .
#       CHROM:A*01:01:01:01 sample:A*01:01:01:04 POS:9 ID:. REF:GA ALT:G QUAL. FILTER:PASS INFO:. FORMAT:GT VariantType.deletion

        #true                    CAGGAGCAGA GGGGTCAGGG
        #pos                     1234567890 1234567890
        seq4 = "A*01:01:01:05    ---------- .--------- ---------- ---------- "
#       result4:A*01:01:01:05    11    .    AG    A    .    PASS    .
#       CHROM:A*01:01:01:01 sample:A*01:01:01:05 POS:10 ID:. REF:AG ALT:A QUAL. FILTER:PASS INFO:. FORMAT:GT VariantType.deletion

        # multiple SNP deletes
        #true                    CAGGAGCAGA GGGGTCAGGG ATTGGGGAGT 
        #pos                     1234567890 1234567890 1234567890
        seq5 = "A*01:01:01:06    ---------. ---------- -..------- ---------- "
#               A*01:01:01:06     9    .    GA     G    .    PASS    .
#               A*01:01:01:06    21    .    ATT    A    .    PASS    .

        msa2vcf = MSAVCF.MSA2VCFImplementation()
        offset = 0
        
        result1List = msa2vcf._parseMSASample(offset, ref, seq1)
        self.logger.info("result1[0]:\n{}".format(result1List[0]))
        r1Expected = "CHROM:A*01:01:01:01 sample:A*01:01:01:02 POS:1 ID:. REF:CA ALT:A QUAL. FILTER:PASS INFO:. FORMAT:GT VariantType.deletion"
        self.assertEqual(str(result1List[0]), r1Expected)

        result2List = msa2vcf._parseMSASample(offset, ref, seq2)
        self.logger.info("result2[0]:\n{}".format(result2List[0]))
        r2Expected = "CHROM:A*01:01:01:01 sample:A*01:01:01:03 POS:5 ID:. REF:AG ALT:A QUAL. FILTER:PASS INFO:. FORMAT:GT VariantType.deletion"
        self.assertEqual(str(result2List[0]), r2Expected)
        
        result3List = msa2vcf._parseMSASample(offset, ref, seq3)
        self.logger.info("result3[0]:\n{}".format(result3List[0]))
        r3Expected = "CHROM:A*01:01:01:01 sample:A*01:01:01:04 POS:9 ID:. REF:GA ALT:G QUAL. FILTER:PASS INFO:. FORMAT:GT VariantType.deletion"
        self.assertEqual(str(result3List[0]), r3Expected) 
        
        result4List = msa2vcf._parseMSASample(offset, ref, seq4)
        self.logger.info("result4[0]:\n{}".format(result4List[0]))
        r4Expected = "CHROM:A*01:01:01:01 sample:A*01:01:01:05 POS:10 ID:. REF:AG ALT:A QUAL. FILTER:PASS INFO:. FORMAT:GT VariantType.deletion"
        self.assertEqual(str(result4List[0]), r4Expected) 
        
        result5List = msa2vcf._parseMSASample(offset, ref, seq5)   
        r5Expected = [
                        "CHROM:A*01:01:01:01 sample:A*01:01:01:06 POS:9 ID:. REF:GA ALT:G QUAL. FILTER:PASS INFO:. FORMAT:GT VariantType.deletion",
                        "CHROM:A*01:01:01:01 sample:A*01:01:01:06 POS:21 ID:. REF:ATT ALT:A QUAL. FILTER:PASS INFO:. FORMAT:GT VariantType.deletion",
                    ]
        for i in range(len(result5List)):
            msas = result5List[i]
            self.logger.info("\nmsas:{}:{}".format(i, msas))
            self.assertEqual(str(msas), r5Expected[i])
        
        self.logger.info("END\n")
        
    ################################################################################ 
    def testParseMSASample_DeleteEndSNP(self):
        self.logger.info("BEGIN")

        # data/A_gen.txt
        # SNP delete at end 
        #true                    CAGGAGCAGA
        #pos                     1234567890
        ref  = "A*01:01:01:01    CAGGAGCAGA"        
        seq1 = "A*01:01:01:04    ---------."
#       result3:A*01:01:01:04    10    .    GA    G    .    PASS    .
#       CHROM:A*01:01:01:01 sample:A*01:01:01:04 POS:9 ID:. REF:GA ALT:G QUAL. FILTER:PASS INFO:. FORMAT:GT VariantType.deletion
          
        msa2vcf = MSAVCF.MSA2VCFImplementation()
        offset = 0
        
        result1List = msa2vcf._parseMSASample(offset, ref, seq1)
        self.logger.info("result1[0]:\n{}".format(result1List[0]))
        r1Expected = "CHROM:A*01:01:01:01 sample:A*01:01:01:04 POS:9 ID:. REF:GA ALT:G QUAL. FILTER:PASS INFO:. FORMAT:GT VariantType.deletion"
        self.assertEqual(str(result1List[0]), r1Expected)
        
                
        self.logger.info("END\n")
                
    ################################################################################ 
    def testDeleteRunAtBegining(self):
        self.logger.info("BEGIN")
        # true         CAG
        # pos          1234567890
        ref  = "ref    CAGGAGCAGA"
        seq1 = "s1     ...------T"
#       result1:A*01:01:01:02    1    .    C    .    .    PASS    .
#       CHROM:A*01:01:01:01 sample:A*01:01:01:02 POS:7 ID:. REF:CAGA ALT:C QUAL. FILTER:PASS INFO:. FORMAT:GT VariantType.deletion
        msa2vcf = MSAVCF.MSA2VCFImplementation()
        offset = 0        
        
        msaList = msa2vcf._parseMSASample(offset, ref, seq1)
        for msa in msaList:
            self.logger.info("msa:{}".format(msa))
            
        self.logger.info("msaList:\n{}".format(msaList))
        e1 = "CHROM:ref sample:s1 POS:1 ID:. REF:CAGG ALT:G QUAL. FILTER:PASS INFO:. FORMAT:GT VariantType.deletion"
        self.assertEqual(str(msaList[0]), e1)        
        
        e2 = "CHROM:ref sample:s1 POS:10 ID:. REF:A ALT:T QUAL. FILTER:PASS INFO:. FORMAT:GT VariantType.polymorphism"
        self.assertEqual(str(msaList[1]), e2) 
                   
        self.logger.info("END\n")
        
    ################################################################################ 
    def testDeleteRunAtEnd(self):
        self.logger.info("BEGIN")
        self.logger.info("BEGIN")
        
        # true         CAGGAGCAGA
        # pos          1234567890
        ref  = "ref    CAGGAGCAGA"
        seq1 = "s1    -------..."
#       result1:A*01:01:01:02    1    .    C    .    .    PASS    .
#       CHROM:ref sample:seq1 POS:7 ID:. REF:CAGA ALT:C QUAL. FILTER:PASS INFO:. FORMAT:GT VariantType.deletion
        msa2vcf = MSAVCF.MSA2VCFImplementation()
        offset = 0        
        
        msaList = msa2vcf._parseMSASample(offset, ref, seq1)
        self.logger.info("msaList:\n{}".format(msaList))
        e1 = "CHROM:ref sample:s1 POS:7 ID:. REF:CAGA ALT:C QUAL. FILTER:PASS INFO:. FORMAT:GT VariantType.deletion"
        self.assertEqual(str(msaList[0]), e1)        
        
        self.logger.info("END\n")
        
    ################################################################################         
    def testDeleteRun(self):
        self.logger.info("BEGIN")
        
        # true        GATT
        # pos         1234
        ref = "ref    GATT"
        s1  = "s1     -..-" # middle
        s2  = "s1     -..." # end
        s3  = "s1     ...-" # beginning

        msa2vcf = MSAVCF.MSA2VCFImplementation()
        offset = 0        
        
        msaList = msa2vcf._parseMSASample(offset, ref, s1)
        self.logger.info("msaList:\n{}".format(msaList))
        e1 = "CHROM:ref sample:s1 POS:1 ID:. REF:GAT ALT:G QUAL. FILTER:PASS INFO:. FORMAT:GT VariantType.deletion"
        self.assertEqual(str(msaList[0]), e1)
        
        msaList = msa2vcf._parseMSASample(offset, ref, s2)
        self.logger.info("msaList:\n{}".format(msaList))       
        e2 = "CHROM:ref sample:s1 POS:1 ID:. REF:GATT ALT:G QUAL. FILTER:PASS INFO:. FORMAT:GT VariantType.deletion" 
        self.assertEqual(str(msaList[0]), e2)
        
        msaList = msa2vcf._parseMSASample(offset, ref, s3)
        self.logger.info("msaList:\n{}".format(msaList))               
        e3 = "CHROM:ref sample:s1 POS:1 ID:. REF:GATT ALT:T QUAL. FILTER:PASS INFO:. FORMAT:GT VariantType.deletion"
        self.assertEqual(str(msaList[0]), e3)

        self.logger.info("END\n")
        
    ################################################################################         
    def testCreateAbstractSyntaxTree(self):
        self.logger.info("BEGIN")
        
        ref  = "r1    GAT"
        s1   = "s1    -T-"
        s2   = "s2    -C-"
        s3   = "s3    -T-"
        
        msa2vcf = MSAVCF.MSA2VCFImplementation()
        offset = 0
        
        msasList1  = msa2vcf._parseMSASample(offset, ref, s1)
        msasList2  = msa2vcf._parseMSASample(offset, ref, s2)
        msasList3  = msa2vcf._parseMSASample(offset, ref, s3)
        
        msasList = msasList1 + msasList3 + msasList2 # test sampleNamesList is sorted correctly
        for msas in msasList:
            self.logger.info("msas:{}".format(msas))
        sampleNamesList, sampleNameLookUp, AST = msa2vcf._createAbstractSyntaxTree(msasList)
        
        self.logger.info("sampleNamesList:{}".format(sampleNamesList))
        self.assertEqual(sampleNamesList, ['s1', 's2', 's3'])
        
        self.logger.info("sampleNameLookUp:{}".format(sampleNameLookUp))
        self.assertEqual(sampleNameLookUp, {'s1': 0, 's2': 1, 's3': 2})
        
        self.logger.info("ast:{}".format(AST))
        
        pp = pprint.PrettyPrinter()
        #pp.pprint(AST) # debug line comment this line out mess up unit test summary output
        ASTStr = pp.pformat(AST)
        
        deletionVCF = AST['r1POS:2'][MSAVCF.VariantType.deletion]
        self.logger.info("deletionVCF:\n{}".format(deletionVCF))
        self.assertEqual(deletionVCF, {})

        insertionsVCF = AST['r1POS:2'][MSAVCF.VariantType.insertion]
        self.logger.info("insertionsVCF:\n{}".format(insertionsVCF))
        self.assertEqual(insertionsVCF, {})

        unknownVCF = AST['r1POS:2'][MSAVCF.VariantType.unknown]
        self.logger.info("unknownVCF:\n{}".format(unknownVCF))
        self.assertEqual(unknownVCF, {})        
        
        polymorphismVCF = AST['r1POS:2'][MSAVCF.VariantType.polymorphism]
        self.logger.info("polymorphismVCF:\n{}".format(polymorphismVCF))
        expectedP = "{'T': [CHROM:r1 sample:s1 POS:2 ID:. REF:A ALT:T QUAL. FILTER:PASS INFO:. FORMAT:GT VariantType.polymorphism, CHROM:r1 sample:s3 POS:2 ID:. REF:A ALT:T QUAL. FILTER:PASS INFO:. FORMAT:GT VariantType.polymorphism], 'C': [CHROM:r1 sample:s2 POS:2 ID:. REF:A ALT:C QUAL. FILTER:PASS INFO:. FORMAT:GT VariantType.polymorphism]}"
        self.assertEqual(str(polymorphismVCF), expectedP)
        
        self.logger.info("END\n")        
        
    ################################################################################ 
    def testParsemsaSampleList(self):
        self.logger.info("BEGIN")
        
        ref  = "r1    GAT"
        s1   = "s1    -T-"
        s2   = "s2    -C-"
        s3   = "s3    -T-"
        s4   = "s4    -.-"
        s5   = "s5    -.-"
        
        msa2vcf = MSAVCF.MSA2VCFImplementation()
        offset = 0
        
        vcfList, sampleNamesList = msa2vcf._parseMSASampleList(offset, ref, [s1, s2, s3, s4, s5])
        expected = [
            "r1    1    .    GA    G    .    PASS    .    GT    0 0 0 1 1",                        
            "r1    2    .    A    T    .    PASS    .    GT    1 0 1 0 0",
            "r1    2    .    A    C    .    PASS    .    GT    0 1 0 0 0",
            ]
        
        for i in range(len(vcfList)):
            vcf = vcfList[i]
            self.logger.info("i:{} vcf:{}".format(i, vcf))
            self.logger.info("       e:{}\n".format(expected[i]))
            e = expected[i]
            cleanE = "".join(e.split()) # remove whitespace
            cleanV = "".join(str(vcf).split())
            self.assertEqual(cleanV, cleanE, " i:{}".format(i))
        
        self.logger.info("END\n")     
        
    ################################################################################ 
    def testUnknown(self):
        self.logger.info("BEGIN")
        # TRUE       GAT
        # pos        123
        ref = "r1    GAT"
        s1  = "s1    -*-"
        s2  = "s2    **-"
        s3  = "s3    ***"
        
        msa2vcf = MSAVCF.MSA2VCFImplementation()
        offset = 0
        
        msasList1  = msa2vcf._parseMSASample(offset, ref, s1)
        self.logger.info("msasList1:\n{}".format(msasList1))
        self.assertEqual(msasList1, [])
        
        msasList2  = msa2vcf._parseMSASample(offset, ref, s3)
        self.logger.info("msasList2:\n{}".format(msasList2))
        self.assertEqual(msasList2, [])

        msasList3  = msa2vcf._parseMSASample(offset, ref, s3)
        self.logger.info("msasList3:\n{}".format(msasList3))
        self.assertEqual(msasList3, [])        

        vcfList, sampleNamesList = msa2vcf._parseMSASampleList(offset, ref, [s1, s2, s3])
        self.assertEqual(vcfList, [])    
                
        self.logger.info("END\n")
   
    ################################################################################ 
    def testInsert(self):
        self.logger.info("BEGIN")
        # TRUE       G  C
        # pos        1  2
        ref = " r    G.C"
        s1  = "s1    -A-"
        s2  = "s2    -A-"
        s3  = "s3    -G-"
        
        expected = [
            "r 1 . G GA . PASS . GT 1 1 0",
            "r 1 . G GG . PASS . GT 0 0 1",
            ]
        
        msa2vcf = MSAVCF.MSA2VCFImplementation()
        offset = 0
        
        vcfList, sampleNameList = msa2vcf._parseMSASampleList(offset, ref, [s1, s2, s3])
        for i in range(len(vcfList)):
            vcf = vcfList[i]
            self.logger.info("i:{} vcf:{}".format(i, vcf))
            e = str(expected[i])
            self.logger.info("      e:{}\n".format(e))            
#             cleanE = "".join(e.split()) # remove whitespace
            self.logger.info("AEDWIP:{}".format(str(vcf).split()))
            cleanV = " ".join(str(vcf).split()) 
            self.logger.info("cleanV:{}".format(cleanV))
            self.logger.info("cleanE:{}".format(e))
            self.assertEqual(cleanV, e)               
                
        self.logger.info("END\n")

    ################################################################################ 
    def testInsertAtBegining(self):
        self.logger.info("BEGIN")
        
        # TRUE              A
        # POS               1234
        ref = "ref1     ....ACTA"
        s1   =   "s1    GATT---G"
        
        offset = 0                
        msa2vcf = MSAVCF.MSA2VCFImplementation()        
        msaList = msa2vcf._parseMSASample(offset, ref, s1)
        for msa in msaList:
            self.logger.info("msa:\n{}".format(msa))     
        # TODO: check with Jonas is position correct?   
        e1 = "CHROM:ref1 sample:s1 POS:1 ID:. REF:A ALT:GATTA QUAL. FILTER:PASS INFO:. FORMAT:GT VariantType.insertion"
        self.assertEqual(str(msaList[0]), e1)       
        
        e2 = "CHROM:ref1 sample:s1 POS:4 ID:. REF:A ALT:G QUAL. FILTER:PASS INFO:. FORMAT:GT VariantType.polymorphism"
        self.assertEqual(str(msaList[1]), e2)       
 
        self.logger.info("END")
        
        
    ################################################################################ 
    def testInsertRunAtEND(self):
        self.logger.info("BEGIN")
        
        # TRUE          GATTA
        # POS           12345
        ref3 = "ref1    GATT.."
        s3   =   "s3    ----AC"
        
        offset = 0                
        msa2vcf = MSAVCF.MSA2VCFImplementation()        
        msaList = msa2vcf._parseMSASample(offset, ref3, s3)
        self.logger.info("msaList:\n{}".format(msaList))     
        # TODO: check with Jonas is position correct?   
        e3 = "CHROM:ref1 sample:s3 POS:4 ID:. REF:T ALT:TAC QUAL. FILTER:PASS INFO:. FORMAT:GT VariantType.insertion"
        self.assertEqual(str(msaList[0]), e3)        
        self.logger.info("END")      
          
    ################################################################################         
    def testInsertRun(self):
        self.logger.info("BEGIN")
        
        # TRUE          G   T
        # pos           1   2
        ref1 = "ref1    G..T"
        s1     = "s1    -AT-"

        msa2vcf = MSAVCF.MSA2VCFImplementation()
        offset = 0        
        
        msaList = msa2vcf._parseMSASample(offset, ref1, s1)
        self.logger.info("msaList:\n{}".format(msaList))
        e1 = "CHROM:ref1 sample:s1 POS:1 ID:. REF:G ALT:GAT QUAL. FILTER:PASS INFO:. FORMAT:GT VariantType.insertion"
        self.assertEqual(str(msaList[0]), e1)
        
        # TRUE          G
        # pos           1
        ref2 = "ref1    G..."
        s2   =   "s2    -ATT"
        msaList = msa2vcf._parseMSASample(offset, ref2, s2)
        self.logger.info("msaList:\n{}".format(msaList))
        e2 = "CHROM:ref1 sample:s2 POS:1 ID:. REF:G ALT:GATT QUAL. FILTER:PASS INFO:. FORMAT:GT VariantType.insertion"
        self.assertEqual(str(msaList[0]), e2)

        
        # TRUE          G   T
        # POS           1   2
        ref4 = "ref1    G..T"
        s4   =   "s4    -ATC"
        msaList = msa2vcf._parseMSASample(offset, ref4, s4)
        e4a = "CHROM:ref1 sample:s4 POS:1 ID:. REF:G ALT:GAT QUAL. FILTER:PASS INFO:. FORMAT:GT VariantType.insertion"
        e4b = "CHROM:ref1 sample:s4 POS:2 ID:. REF:T ALT:C QUAL. FILTER:PASS INFO:. FORMAT:GT VariantType.polymorphism"
        self.logger.info("msaList:\n{}".format(msaList))     
        self.assertEqual(str(msaList[0]), e4a)

        self.assertEqual(str(msaList[1]), e4b)

    ################################################################################ 
    def testTricky(self):
        self.logger.info("BEGIN")
        
        offset = 0
        msa2vcf = MSAVCF.MSA2VCFImplementation()   
        
        # true          GAATT
        # pos           12345
        ref1 = "ref1    GAATT"
        s1     = "s1    **--C"
        msaList = msa2vcf._parseMSASample(offset, ref1, s1)
        self.logger.info("msaList:\n{}".format(msaList))        
        e1 = "CHROM:ref1 sample:s1 POS:5 ID:. REF:T ALT:C QUAL. FILTER:PASS INFO:. FORMAT:GT VariantType.polymorphism"
        self.assertEqual(str(msaList[0]), e1)

        # true          GAATT GA
        # pos           12345 67
        ref2 = "ref2    GAATT.GA"
        s2     = "s2    **--C.-C"        

        msaList = msa2vcf._parseMSASample(offset, ref2, s2)
        self.logger.info("msaList:\n{}\n".format(msaList))
        e2 = "CHROM:ref2 sample:s2 POS:5 ID:. REF:T ALT:C QUAL. FILTER:PASS INFO:. FORMAT:GT VariantType.polymorphism" 
        self.assertEqual(str(msaList[0]), e2)
        
        e3 = "CHROM:ref2 sample:s2 POS:7 ID:. REF:A ALT:C QUAL. FILTER:PASS INFO:. FORMAT:GT VariantType.polymorphism"
        self.assertEqual(str(msaList[1]), e3)
        
        # true          GAATT GA G
        # pos           12345 67 8
        ref4 = "ref4    GAATT.GA*G"
        s4     = "s4    **--C.-CAT"  
        msaList = msa2vcf._parseMSASample(offset, ref4, s4)
        self.logger.info("msaList:\n{}\n".format(msaList))   
        e4 = "CHROM:ref4 sample:s4 POS:5 ID:. REF:T ALT:C QUAL. FILTER:PASS INFO:. FORMAT:GT VariantType.polymorphism" 
        self.assertEqual(str(msaList[0]), e4)
        
        e5 = "CHROM:ref4 sample:s4 POS:7 ID:. REF:A ALT:C QUAL. FILTER:PASS INFO:. FORMAT:GT VariantType.polymorphism" 
        self.assertEqual(str(msaList[1]), e5)
        
        e6 = "CHROM:ref4 sample:s4 POS:8 ID:. REF:G ALT:T QUAL. FILTER:PASS INFO:. FORMAT:GT VariantType.polymorphism"
        self.assertEqual(str(msaList[2]), e6)
                
        self.logger.info("END\n")

    ################################################################################ 
    def testTricky2(self):
        self.logger.info("BEGIN")
        # true        GA CT
        # ps          12 34
        ref1 = "r1    GA|CT"
        s1   = "s1    ---A-"
        
        offset = 0
        msa2vcf = MSAVCF.MSA2VCFImplementation()           
        msaList1 = msa2vcf._parseMSASample(offset, ref1, s1)
        self.logger.info("msalist1:\n{}".format(msaList1))
        e1 = "CHROM:r1 sample:s1 POS:3 ID:. REF:C ALT:A QUAL. FILTER:PASS INFO:. FORMAT:GT VariantType.polymorphism"
        self.assertEqual(str(msaList1[0]), e1)
        self.logger.info("END\n")
           

    ################################################################################ 
    def testVGConstructBug(self):
        """
        index file shortenBug2299_1,2180d.out.fasta.fai not found, generating...
error:[vg::Constructor] Variant/reference sequence mismatch: TC vs pos: 111: CC; do your VCF and FASTA coordinates match?
Variant: ref    111    .    TC    T    0    PASS    .
zero ind: 110 1-indexed: 111
        """
        self.logger.info("BEGIN")
        # true         A GTCCCC
        # pos          1 234567
        ref1 = "ref    A|GTCCCC"
        s1   = "s1     -|--.---"
        offset = 0
        msa2vcf = MSAVCF.MSA2VCFImplementation()           
        msaList1 = msa2vcf._parseMSASample(offset, ref1, s1)
        self.logger.info("msalist1:\n{}".format(msaList1))
        
        e1 = "CHROM:ref sample:s1 POS:3 ID:. REF:TC ALT:T QUAL. FILTER:PASS INFO:. FORMAT:GT VariantType.deletion"
        self.assertEqual(str(msaList1[0]), e1)
        
        self.logger.info("END\n")

    ################################################################################ 
    def testBug287(self):
        """
        E_gen.txt 
        error:[vg::Constructor] Variant/reference sequence mismatch: A vs pos: 287: C; do your VCF and FASTA coordinates match?
Variant: E*01:01:01:01    287    .    A    AA    0    PASS    .
zero ind: 286 1-indexed: 287
        """
        self.logger.info("BEGIN")
        
        ref1 = "ref    A.G"
        s1   =  "s1    ***"
        
        offset = 0
        msa2vcf = MSAVCF.MSA2VCFImplementation()           
        msaList1 = msa2vcf._parseMSASample(offset, ref1, s1)
        self.logger.info("msalist1:\n{}".format(msaList1))
                
        self.assertEqual(msaList1, [])
        
        # true          AGCA G
        #pos            1234 5
        ref2 = "ref2    AGCA.G"
        s2   =  "s2     **--C-"        
        msaList2 = msa2vcf._parseMSASample(offset, ref2, s2)
        self.logger.info("msalist2:\n{}".format(msaList2))        
        
        e2 = "CHROM:ref2 sample:s2 POS:4 ID:. REF:A ALT:AC QUAL. FILTER:PASS INFO:. FORMAT:GT VariantType.insertion"
        self.assertEqual(str(msaList2[0]), e2)

        self.logger.info("END\n")    
        
    ################################################################################ 
    def testTemplate3579(self):
        """
        + ../vg/bin/vg construct -r shortenBug3579_1,0d.out.fasta -v shortenBug3579_1,0d.out.vcf.gz
index file shortenBug3579_1,0d.out.fasta.fai not found, generating...
error:[vg::Constructor] Variant/reference sequence mismatch: C vs pos: 724: T; do your VCF and FASTA coordinates match?
Variant: ref    724    .    C    T    0    PASS    .
        """
        self.logger.info("BEGIN")
        
        # test cased based on file E*01:01:01:01.sequence.w100-E*01:03:02:07.sequence.w100 
        #true         AACC GG TT AC GTC
        #pos          1234 56 78 90 123
        ref1 = "r1    AACC.GG|TT|AC|GTC"
        s1   = "s1    **--.--|--|--|--T"
        
        offset = 0
        msa2vcf = MSAVCF.MSA2VCFImplementation()           
        msaList1 = msa2vcf._parseMSASample(offset, ref1, s1)  
        self.logger.info("msaList1:\n{}".format(msaList1))      
       
        e1 = "CHROM:r1 sample:s1 POS:13 ID:. REF:C ALT:T QUAL. FILTER:PASS INFO:. FORMAT:GT VariantType.polymorphism"
        self.assertEqual(str(msaList1[0]), e1)

        self.logger.info("END\n")        
      
    ################################################################################ 
    def testAGenBug288(self):
        """
        A_gen.txt. bug at position 288 in sample A*68:02:01:0
        """
        self.logger.info("BEGIN")
        
        # true        CCC AGACGCCGAGG
        # pos         123 45678901234
        ref1 = "r1    CCC.AGACGCCGAGG"
        s1   = "s1    --............-"        

        offset = 0
        msa2vcf = MSAVCF.MSA2VCFImplementation()           
        msaList1 = msa2vcf._parseMSASample(offset, ref1, s1)  
        self.logger.info("msaList1:\n{}".format(msaList1))  
                
        e1 = "CHROM:r1 sample:s1 POS:2 ID:. REF:CCAGACGCCGAG ALT:C QUAL. FILTER:PASS INFO:. FORMAT:GT VariantType.deletion"
        self.assertEqual(str(msaList1[0]), e1)
        
        self.logger.info("END\n")
      


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
