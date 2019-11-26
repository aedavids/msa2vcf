#!/bin/sh 

#
# make it easier to test from cli
# aedavids@ucsc.edu
#


#set -x  # turns shell debugging on
#set +x # turns shell debugging off

fullPathToSelf=$0
if [ $# -ne 2 ]; then
        scriptName=`basename $0`
    echo "error: usage $scriptName multipleSequenceAlignementsPath reference.fasta"
    echo "missing argument"
    exit 1
fi

which vg
foundVG=$?
if [ $foundVG -ne 0 ]; then
	echo "error vg not found in path"
	exit $foundVG
fi 

set -x  # turns shell debugging on
#set +x # turns shell debugging off

msfPath=$1
refFastaPath=$2

#
# use filePrefix to find artifacts
# example:
#	fileName="E_gen.txt"
#	filePrefix="E_gen"
msfFileName=`basename $msfPath`
msfFilePrefix=${msfFileName%.*}

vcfFile=$msfFilePrefix.vcf
dotFile=$msfFilePrefix.dot
pdfFile=$msfFilePrefix.vis.pdf
vgFile=$msfFilePrefix.vg
gzFile=$vcfFile.gz
tbiFile=$gzFile.tbi

#set -x  # turns shell debugging on
set +x # turns shell debugging off

#
# clean up any artifacts from previous runs
#
#set -x  # turns shell debugging on
set +x # turns shell debugging off

'rm' info.log errors.log
'rm' $vcfFile $gzFile $dotFile $pdfFile $vgFile $tbiFile errors.log info.log

# make sure vg construct creates a new index file
'rm' $refFastaPath.fai

set -x  # turns shell debugging on
#set +x # turns shell debugging off

#
# convert the multiple alignments to vcf 
#
echo --------- RUN CONVERTER --------------
msf2vcf $msfPath $vcfFile
echo

#echo ------------- FIX REF CHROM NAME -------------
##
## AEDWIP hack vg can not deal with refernce names like E*01:01:01:01
## 10/25/19 adam says we do not need to fix this, just make sure contig name is same as fasta header
##  hack the reference fasta so that the name is 'ref'
##(HLA_haplotype) $ grep contig E_gen.txt.vcf
####contig=<ID=E*01:01:01:01,length=3909>
##
##(HLA_haplotype) $ grep contig E_gen.txt.vcf | cut -d '=' -f 3
##E*01:01:01:01,length
##
##(HLA_haplotype) $ grep contig E_gen.txt.vcf | cut -d '=' -f 3 | cut -d , -f 1
##E*01:01:01:01
#origRefChromName=`grep contig $vcfFile | cut -d '=' -f 3 | cut -d , -f 1`
##echo 'E*01:01:01:01' | sed 's/\*/\\*/g'
#refChromName=`echo $origRefChromName | sed 's/\*/\\\*/g'`
#echo AEDWIP new refChromName = $refChromName
##cat $vcfFile | sed 's/E\*01:01:01:01/AEDWIP/g' >$vcfFile
#
#
#cat $vcfFile | sed "s/$refChromName/ref/g" > tmp
#mv tmp $vcfFile
#echo

## now undo the first replacement
## https://unix.stackexchange.com/a/188269
##cat $vcfFile | sed "s/ref/$refChromName/" > tmp
##awk 'NR==1,/ref/{sub(/ref/, "FOO")} 1' $vcfFile > tmp
#awk -v foo=$origRefChromName 'NR==1,/ref/{sub(/ref/, foo)} 1' $vcfFile > tmp
#mv tmp $vcfFile	


# 
# construct graph
#
echo ------- CONSTRUCT GRAPH ------------
# bgzip will create an output file that starts with name of inpput file and ending in '.gz'
bgzip $vcfFile
# tabix will create a file that starts with name of input file ends in '.tbi`
tabix -p vcf $gzFile

# vg construct redirect output into a file we control the name of the output file
# typical it ends in '.vg'
vg construct -r $refFastaPath -v $gzFile > $vgFile
exitStatus=$?
if [ $exitStatus -ne 0 ]; then
	echo ERROR vg contruct failed with exitCode $exitStatus
	# uncompress the file
	bgzip -d $gzFile	
	exit $exitStatus
fi
#../vg/bin/vg construct -r $refFastaPath -v $vcfFile > $vgFile

# uncompress the file
bgzip -d $gzFile
echo

#
# construct the dot file
#
echo ------------- VIEW GRAPH ---------------
vg view -d $vgFile > $dotFile
echo

#
# construct the pdf version of the graph
#
echo ------- CREATE PDF ---------
dot -Tpdf $dotFile > $pdfFile

open $pdfFile
 