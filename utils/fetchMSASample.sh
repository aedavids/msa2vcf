#!/bin/sh 



set -x  # turns shell debugging on
#set +x # turns shell debugging off

fullPathToSelf=$0
if [ $# -ne 2 ]; then
        scriptName=`basename $0`
    echo "error: usage $scriptName multipleSequenceAlignementsPath sampleName"
    echo "missing argument"
    exit 1
fi

msaPath=$1
sampleNameArg=$2
sampleName=`echo $sampleNameArg | sed 's/\*/\\\*/g'`

#
# grep finds sample lines
# awk remove the first col. i.e. the sample name
# tr remove the spaces with in the sequence
# awk combines the lines
#

# the pattern to grep must end with a space.
# there are sample name like 'E*01:10' we want to make sure we do not pick up 'E*01:01:01:10'
grep "${sampleName}[[:blank:]]" $msaPath | awk '{$1=""; print $0}' | tr -d ' ' | awk '{printf("%s",$0)}'


# to debug problem at position 432
# take output from ths script and use
# following to select sequence in positions 415 to 425
# cut -b 415-10

# (HLA_haplotype) $ grep 'E\*01:10[[:blank:]]' E_gen.txt
