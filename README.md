# msa2vcf
Converts a multiple sequence alignment file in the format described at https://www.ebi.ac.uk/ipd/imgt/hla/nomenclature/index.html to VCF.

## references:
* [https://www.ebi.ac.uk/ipd/imgt/hla/nomenclature/alignments.html](https://www.ebi.ac.uk/ipd/imgt/hla/nomenclature/alignments.html)
* [https://www.ebi.ac.uk/ipd/imgt/hla/nomenclature/index.html](https://www.ebi.ac.uk/ipd/imgt/hla/nomenclature/index.html)


## Code / Dependencies:
* The code is written in Python 3.
* Shared and distributed via Jupyter Notebooks.

## Data
* [ftp://ftp.ebi.ac.uk/pub/databases/ipd/imgt/hla/alignments/](ftp://ftp.ebi.ac.uk/pub/databases/ipd/imgt/hla/alignments/)

## Quick Start:



```
$ msa2vcf/msa2vcf 
usage: inputFile outputFile [options]

Converts a multiple sequence alignment file in the format described at
https://www.ebi.ac.uk/ipd/imgt/hla/nomenclature/index.html to VCF. See
logging.ini.json to change log level

positional arguments:
  inputFile   The file in MSF format.
  outputFile  The VCF version of the file

optional arguments:
  -h, --help  show this help message and exit
  --version   show program's version number and exit
  
  
$ cd msa2vcf/
$ export PYTHONPATH=`pwd`
$ msa2vcf/msa2vcf test/data/E_gen.txt ./E_gen.vcf
```

## Running unit tests
```
$ cd msa2vcf/
$ ls
LICENSE          README.md     msa2vcf/         setupLogging.py  test/
(base) $ 
$ export PYTHONPATH=`pwd`
$ cd test
$ python3 -m unittest discover -v
```

## How to Debug
see [howToDebugPositionBugs.md](howToDebugPositionBugs.md)
