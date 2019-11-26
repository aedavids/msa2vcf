'''
Created on Oct 10, 2019

@author: Andrew Davidson (aedavids@ucsc.edu)
'''
import datetime
import logging
from msa2vcf.MSASample import MSASample, VariantType
from msa2vcf.VCFSample import VCFSample

class MSA2VCFImplementation(object):
    '''
    TODO:AEDWIP
    
    assume input file came from ftp://ftp.ebi.ac.uk/pub/databases/ipd/imgt/hla/alignments/
    
    ref: 
        https://samtools.github.io/hts-specs/VCFv4.3.pdf
        https://www.ebi.ac.uk/ipd/imgt/hla/nomenclature/alignments.html
        https://www.ebi.ac.uk/ipd/imgt/hla/nomenclature/index.html
    
    public functions:
        __init__
        parse(inputFilePath, outputFilePath)
    '''

    logger = logging.getLogger(__name__)

    ################################################################################            
    def parse(self, inputFilePath, outputFilePath):
        """
        assume input file came from ftp://ftp.ebi.ac.uk/pub/databases/ipd/imgt/hla/alignments/
        """        
        self.logger.info("BEGIN")   
        
        #FIXME: AEDWIP: make numCols a CLI argument
        numColens = 3
        referenceName, sequenceDict = self._load(inputFilePath, numColens)
        
        offset = 0 #-1 #1 #0
        referenceSequence = sequenceDict.pop(referenceName, None) # remove from dictionary
        
        # calculate the length of the reference sequence. Make sure we do not count
        # Multiple sequence aliment symbols for insertion and "|"
        cleanReferenceSequence = referenceSequence.replace(".", "")
        cleanReferenceSequence = cleanReferenceSequence.replace("|", "")
        referenceLength = len(cleanReferenceSequence)
        self.logger.debug("referenceSequence:\n{}".format(referenceSequence))
        self.logger.debug("len(referenceSequence):{}".format(len(referenceSequence)))
        self.logger.debug("len(cleanReferenceSequence):{}".format(referenceLength))                
        self.logger.debug("cleanReferenceSequence:\n{}".format(cleanReferenceSequence))

        # rebuild reference to match original file format
        reference = "{}\t{}".format(referenceName, referenceSequence)
        
        # rebuild sample list to match original file format
        sampleList = []
        for name, sequence in sequenceDict.items():
            sampleList.append("{}\t{}".format(name,sequence))
        vcfList, sampleNamesList = self._parseMSASampleList(offset, reference, sampleList)
        
        self._write2VCF(outputFilePath, vcfList, sampleNamesList, referenceName, referenceLength)

        self.logger.info("END\n")                

    ################################################################################            
    def _load(self, inputFilePath, numColensInSampleName=3):       
        """
        parse and pre-process
        turns out blocks in input file are not independent. 
        easy fix is to read in all samples and concatenate them into a single sample
        E.G. given a sample "s1" find all its segments in the various blocks and combine the
        
        arguments:
            inputFilePath
            
            numColensInSampleName: in
                ref: http://hla.alleles.org/nomenclature/naming.html
                example HLA-A*02:101:01:020
                3 ":" means full length sample includes exons and introns
                    
        return (referenceName, sequenceDict)
        """
        
        self.logger.info("BEGIN")
        with open(inputFilePath) as infp:            
            sequenceDict = dict() # key = sample name, value is sequence'
            referenceName = None
            for rawLine in infp:
                #strip any leading whitespace
                line = rawLine.lstrip()
                self.logger.info("line:XXX{}XXX".format(line))
            
                if line == None or line =="" or line.isspace() or line[0] == "#" \
                    or "gDNA" in line or line[0] == "|":
                    continue
                
                if "Please" in line:
                    # Please see http://hla.alleles.org/terms.html for terms of use.
                    self.logger.info("found please")
                    continue                
            
                if not referenceName:
                    # first sequence
                    referenceName, referenceSequence = self._cleanInputLine(line)
                    sequenceDict[referenceName] = referenceSequence
                else :
                    sequenceName, sequence = self._cleanInputLine(line)
                    if numColensInSampleName != sequenceName.count(":"):
                        continue
                    
                    if sequenceName not in sequenceDict:
                        sequenceDict[sequenceName] = sequence  
                    else:
                        sequenceDict[sequenceName] += sequence
                
        self.logger.info("END\n")                
        return (referenceName, sequenceDict)

        
    ################################################################################    
    def _cleanInputLine(self, line):
        """
        input:
            line
            example: "A*01:01:01:01     CAGGAGCAGA GGGGTCAGGG ATTGGGGAGT CCCAGCCTTG "

        return
            example: ("A*01:01:01:01", "CAGGAGCAGAGGGGTCAGGGATTGGGGAGTCCCAGCCTTG")
        """
        ll = line.strip().split()
        name = ll[0]
        seq = "".join(ll[1:])
        return(name, seq)
            
        
    ################################################################################            
    def __init__(self):
        '''
        Constructor
        '''
        pass
        
    ################################################################################                
    def _combine(self, numberOfSamples, sampleNameLookUp, variantType, msfsList):
        """
        A production rule. I.E. rewrites symbols from MSF to VCF format

        
        arguments:
            numberOfSamples: int. the length of the number of unique sample names. 
                            Used to create logic vector identifying which samples match the msfc object
                            
            sampleNameLookUp: dict(). key = sampleName, value = index into logic vector.
                            
            varientType: VarientType enumeration
                TODO: we probably do not need this. wait till we have everything working before we 
                clean up
            
            msfsList:  array of MSASample objects that have same ref, position
            
                example
                    [CHROM:r1 sample:s1 POS:2 ID:. REF:A ALT:T QUAL. FILTER:PASS INFO:. VariantType.polymorphism,
                     CHROM:r1 sample:s3 POS:2 ID:. REF:A ALT:T QUAL. FILTER:PASS INFO:. VariantType.polymorphism]
                                                   
                {'C': {'C': [CHROM:r1 sample:s2 POS:2 ID:. REF:A ALT:C QUAL. FILTER:PASS INFO:. VariantType.polymorphism]},
                       'T': [CHROM:r1 sample:s1 POS:2 ID:. REF:A ALT:T QUAL. FILTER:PASS INFO:. VariantType.polymorphism,
                             CHROM:r1 sample:s3 POS:2 ID:. REF:A ALT:T QUAL. FILTER:PASS INFO:. VariantType.polymorphism]}
                }

            
        returns: VCFSample object
        """
        lv = [0] * numberOfSamples
        for msfs in msfsList:
            idx = sampleNameLookUp[msfs.sampleName]
            lv[idx] = 1
            
        ret = VCFSample(msfs, lv)
        return ret

    ################################################################################            
    def _createAbstractSyntaxTree(self, msfsList):
        """
        creates an abstract syntax tree that makes it easy to create a topological ordering and
        productions for generating a vcf version of the original msf file
        
        arguments:
            msfsList: a list of MSASample objects
            
        returns (sampleNamesList, sampleNameLookUp, AST)
            sampleNamesList is a sorted list sample names in the msfsList
            sampleNameLookUp is a dictionary. key = sampleName value = index into  sampleNameList
            AST is the abstract syntax tree
            
            Keys can be easily sort so that the final VCF file is ordered by position
            
            inner most dictionary makes it easy to create VCF logical vector
            
            general structure of the AST:
                {key=ref+pos : value= {
                                        key=VarientType
                                        value= {
                                                key=ALT
                                                value=[MSASample]
                                      }
                                      
            example:
                input:
                ref    GAT
                 N1    -C-
                 N2    -T-
                 N3    -T-
            
            AST = {ref2 : {
                            polymorphism : {
                                            C : ["msfc for N1"],
                                            T : ["msfs for N2", "msfs for N3"]
                                          }
                          }
                   }
        """
        
        AST = dict() 
        sampleNamesSet = set()
        for msfs in msfsList:
            
            # fetch the variations dictionary for grouping key
            # POS: is not needed however it makes debugging easier
            # I.E. we can decompose the key
            groupingKey = msfs.refChromName + "POS:" + str(msfs.position)
            variationsDict = None
            if groupingKey in AST:
                variationsDict = AST[groupingKey]
            else:
                variationsDict = {
                                VariantType.unknown:{},
                                VariantType.polymorphism:{},
                                VariantType.insertion:{},
                                VariantType.deletion:{}
                }    
                AST[groupingKey] = variationsDict
               
            # fetch the dictionary for the msfs variation type
            varientDict = variationsDict[msfs.variantType]
               
            # add msfs to ALT dictionary
            altKey = msfs.alt
            altList = None # list of MSASample object that have same variant type, ref and position
            if altKey in varientDict :
                altList = varientDict[altKey]
            else:
                altList = [] 
                varientDict[altKey] = altList
                           
            altList.append(msfs) 
            sampleNamesSet.add((msfs.sampleName))
            
        
        # create an index from sample name to logic vector idx
        sampleNamesList = sorted( list(sampleNamesSet) )
        sampleNameLookUp = { sampleNamesList[i]:i for i in range(len(sampleNamesList)) }
        
        return (sampleNamesList, sampleNameLookUp, AST)
    
    ################################################################################ 
    def _parseDeletion(self, bases, i, numBases, position, refSeq, sampleSeq):
        # examples see testParseMSASample_DeleteSNP() and testDeleteRun()
         
        # startPosition is index into true reference sequence
        # i is the index into the MSA format reference sequence
        startPosition = position - 1
        if startPosition <= 0:
            startPosition = 1
         
        endOfDeletionPosition = i
        while endOfDeletionPosition + 1 < numBases and sampleSeq[endOfDeletionPosition + 1] == '.':
            endOfDeletionPosition += 1
            if refSeq[endOfDeletionPosition] in bases:
                position += 1                    
         
        if i - 1 >= 0: 
            msaRefChar = refSeq[i -1:endOfDeletionPosition + 1] 
            msaRefChar = msaRefChar. replace('.', '').replace('|', '') # see testAGenBug288()
            alt = msaRefChar[0]
        else :
            # deletion at beginning
            msaRefChar = refSeq[i:endOfDeletionPosition + 2]
            msaRefChar = msaRefChar. replace('.', '').replace('|', '') # see testAGenBug288()                    
            alt = msaRefChar[-1] #"." # unknown. i.e delete first nucleotide in ref
        
        # FIXME: I do not think we need if/else just use i += (endOfDeletionPosition - i) + 1
        if (endOfDeletionPosition - i) == 0 :
            # single nucleotide deletion
            i += 1
        else:
            # run of deletions
            i += (endOfDeletionPosition - i) + 1
            
        return (i,startPosition, msaRefChar, alt )            


    ################################################################################ 
    def _parseInsertion(self, bases, i, numBases, position, refSeq, sampleSeq):
        # startPosition is index into true reference sequence
        # i is the index into the MSA format reference sequence
        startPosition = position #FIXME: AEDWIP this was cut and past from delet by accident - 1
        if startPosition <= 0:
            startPosition = 1
            
        endInsertionPosition = i
        while endInsertionPosition + 1 < numBases and refSeq[endInsertionPosition + 1] == '.':                        
            endInsertionPosition += 1   
            if refSeq[endInsertionPosition] in bases:  # this should never happen
                position += 1                             
        
        if i - 1 >= 0: 
            msaRefChar = refSeq[i-1]                    
            alt = sampleSeq[i:endInsertionPosition + 1]
            alt = msaRefChar + alt
        else :
            # insertion at begining of reference
            msaRefChar = refSeq[endInsertionPosition + 1] #"." # unknown.                            
            alt = sampleSeq[i:endInsertionPosition + 1] + msaRefChar
#                     msaRefChar = "." # unknown.        
      
        i = endInsertionPosition + 1        
        
        return (i,startPosition, msaRefChar, alt )
    
    ################################################################################ 
    def _parseMSASample(self, offset, reference, sample): 
        '''
        FIXME: add doc
        
         arguments:
            offset: integer position offset. I.E. position of first base in reference
            
            reference:reference sequence line
                example: "A*01:01:01:01     CAGGAGCAGA GGGGTCAGGG ATTGGGGAGT CCCAGCCTTG "
                
            sample: sequence to compare
                example: "A*01:01:01:02    T--------- ---------- ---------- ---------- "
            
        returns a list of MSASample object
        '''
        ret = [] 
        
        # combine all the sequences into a single string with out any white space
        # make it easier to calculate the position
        refList = reference.strip().split()
        refSeq = "".join(refList[1:])
        sampleList = sample.strip().split()
        sampleSeq = "".join(sampleList[1:])
                
        # create place holders for 8 required columns in VCF
        refChromName = refList[0]
        sampleName = sampleList[0]
        if len(refSeq) != len(sampleSeq):
            self.logger.warning("skipping ref len(chrom:{})=={} != len(sampleName:{})=={}"
                                .format(refChromName, len(refSeq),  sampleName, len(sampleSeq)))
            return ret
        
        # clean up to prevent bugs
        refList = None 
        sampleList = None
        
        # init the required vcf values. 
        # The Multiple Sequence alignment format does not include this information
        ident = "." # unknown in VCF is '.'
        qual = "."  # unknown quality score
        filterArg = "PASS" # TODO: ???
        info = "." # empty
        formatArg = "GT"                 

        # position is the location in the a fasta a version reference sequence
        # E.G. nucleotide positions
        position = offset
        varientType = VariantType.unknown
            
        # begin parsing
        bases = set(["A", "C", "G", "T"])
        numBases = len(refSeq)
        i = 0
        
        while i < numBases:
            try:
                msaRefChar = refSeq[i]
            except IndexError as ie:
                self.logger.error("AEDWIP IndexError:\n{}".format(ie))
                
            if msaRefChar in bases or msaRefChar == "*":
                # * means unknown in multiple sequence alignment format
                # probably never occurs
                position += 1
                
            msfChar = sampleSeq[i]
            alt = "ERROR"
            varientType = VariantType.unknown
            
            if msaRefChar == "|" and msfChar == "-":
                i += 1
                            
            elif msfChar == "-":
                # sample matches reference
                i += 1
              
            #####################################################                                
            elif msfChar in bases and msaRefChar in bases:
                # good snp
                alt = msfChar
                varientType = VariantType.polymorphism
                i += 1
                msfs = MSASample(refChromName, sampleName, position, ident, msaRefChar, alt, qual, 
                                    filterArg, info, formatArg, varientType)
                ret.append(msfs)                              

                
            elif msfChar == "*" :
                # "*" means unknown   VCF use '.' i.e. no info  
                i += 1
            
            elif msaRefChar == "*":
                # '*' means unknown. probably never occurs
                i += 1
                position -= 1                                                   
            
            #####################################################                  
            elif msfChar == '.' and msaRefChar in bases:
                # FIXME: unit test fail if we use the refactored version of delete
#                 i, startPosition, msaRefChar, alt = self._parseDeletion(bases, i, numBases, position, refSeq, sampleSeq)                
#                 varientType = VariantType.deletion                
#                 msfs = MSASample(refChromName, sampleName, startPosition, ident, msaRefChar, alt, qual, 
#                                     filterArg, info, formatArg, varientType)
#                 ret.append(msfs)
                                
                varientType = VariantType.deletion
                # examples see testParseMSASample_DeleteSNP() and testDeleteRun()
                  
                # startPosition is index into true reference sequence
                # i is the index into the MSA format reference sequence
                startPosition = position - 1
                if startPosition <= 0:
                    startPosition = 1
                  
                endOfDeletionPosition = i
                while endOfDeletionPosition + 1 < numBases and sampleSeq[endOfDeletionPosition + 1] == '.':
                    endOfDeletionPosition += 1
                    if refSeq[endOfDeletionPosition] in bases:
                        position += 1                    
                  
                if i - 1 >= 0: 
                    msaRefChar = refSeq[i -1:endOfDeletionPosition + 1] 
                    msaRefChar = msaRefChar. replace('.', '').replace('|', '') # see testAGenBug288()
                    alt = msaRefChar[0]
                else :
                    # deletion at beginning
                    msaRefChar = refSeq[i:endOfDeletionPosition + 2]
                    msaRefChar = msaRefChar. replace('.', '').replace('|', '') # see testAGenBug288()                    
                    alt = msaRefChar[-1] #"." # unknown. i.e delete first nucleotide in ref
  
                # FIXME: I do not think we need if/else just use i += (endOfDeletionPosition - i) + 1
                if (endOfDeletionPosition - i) == 0 :
                    # single nucleotide deletion
                    i += 1
                else:
                    # run of deletions
                    i += (endOfDeletionPosition - i) + 1
                      
                msfs = MSASample(refChromName, sampleName, startPosition, ident, msaRefChar, alt, qual, 
                                    filterArg, info, formatArg, varientType)
                ret.append(msfs)                

            #####################################################                  
            elif msaRefChar == '.' and  msfChar in bases:
                i, startPosition, msaRefChar, alt = self._parseInsertion(bases, i, numBases, position, refSeq, sampleSeq)                
                varientType = VariantType.insertion                
                msfs = MSASample(refChromName, sampleName, startPosition, ident, msaRefChar, alt, qual, 
                                    filterArg, info, formatArg, varientType)
                ret.append(msfs)

            elif msfChar == "|" and msaRefChar == "|":
                # does not say how what this symbol means
                # https://www.ebi.ac.uk/ipd/imgt/hla/nomenclature/alignments.html
                # if we do not advance vg construct throws error
                # error:[vg::Constructor] Variant/reference sequence mismatch: G vs pos: 432: A; do your VCF and FASTA coordinates match?
                # Variant: ref    432    .    G    A    0    PASS    .
                # zero ind: 431 1-indexed: 432
                i += 1

            elif msaRefChar == '.' and  msfChar == ".":
                self.logger.info("AEDWIP msaRefChar == msfChar == '.' ignored sampleName:{} offset:{} position:{}"
                                    .format(sampleName, offset, position))
                i += 1    

            else:
                self.logger.warning("unexpected ref char:'{}' msf char:'{}' ignored sampleName:{} offset:{} position:{}"
                                 .format(msaRefChar, msfChar, sampleName, offset, position))
                i += 1
    
        return ret
    
    ################################################################################                    
    def _parseMSASampleList(self, offset, reference, sampleList):
        """
        FIXME: add doc string
        
        arguments
            offset: int. the sample position offset
            reference: string. reference sequence
            sampleList: list of string in https://www.ebi.ac.uk/ipd/imgt/hla/nomenclature/index.html format
        
        returns: (vcfList, sampleNamesList)
                    vcfList is a list of VCFSample objects
                    sampleNameList are the sample names
        """
        ret = []
        
        # parse/ tokenize the samples
        msfsList = []
        for i in range(len(sampleList)):
            msfs = self._parseMSASample(offset, reference, sampleList[i])
            # FIXME: do not append if empty
            msfsList= msfsList + msfs # append the msfs to our master list
         
        # create abstract syntax tree 
        sampleNamesList, sampleNameLookUp, AST = self._createAbstractSyntaxTree(msfsList)
        
        # create topologic ordering
        # vcf format requires sample to be ordered by ref chrom name and position
        # example of a key : "E*01:01:01:01POS:1062"
        orderedSampleKeys = sorted(list(AST.keys()), key=lambda k: int(k.split("POS:")[1]))

        # run productions
        numberOfSamples = len(sampleNamesList)
        for sampleKey in orderedSampleKeys:
            for vt in VariantType:
                msfsDict = AST[sampleKey][vt]
                for msfsList in msfsDict.values():
                    vcf = self._combine(numberOfSamples, sampleNameLookUp, vt, msfsList)
                    ret.append(vcf)
        
        return (ret, sampleNamesList)
        
            
    ################################################################################                
    """
    FIXME: add dco
    """        
    def _write2VCF(self, outputFilePath, vcfList, sampleNamesList, referenceName, referenceLength):
        # "a+" create if it does not exist else append
        with open(outputFilePath, "a+") as fp:
            fp.write("##fileformat=VCFv4.3\n")
            # 20090805
            now = datetime.datetime.now()
            today = "{}{}{}".format(now.year, now.month, now.day)
            fp.write("##fileDate={}\n".format(today))
            fp.write("##contig=<ID={},length={}>\n".format(referenceName, referenceLength))
            fp.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")
                
            #hdr = "\n" + VCFSample.header(sampleNamesList) + "\n"
            hdr = VCFSample.header(sampleNamesList) + "\n"
            fp.write(hdr)
            for vcf in vcfList:
                fp.write(str(vcf) + "\n")
