# gene-puller
Reworked version of Torsten Seeman's gene-puller.pl script. Retrieves aligned gene sequences from assemblies.

##Author

Jason Kwong (@kwongjc)

##Dependencies
* Python 2.7.x
* BLAST+
* MUSCLE
* BioPython

##Usage

```
$ gene-puller.py -h
usage: 
  gene-puller.py --genes FASTA [OPTIONS] FASTA1 FASTA2 FASTA3 ... FASTAN

Searches FASTA files for specified query sequences and outputs ClustalW alignment

positional arguments:
  FASTA          FASTA file to search (required)

optional arguments:
  -h, --help     show this help message and exit
  --genes FASTA  File of query genes in FASTA format (required)
  --out FILE     Output file
  --version      show program's version number and exit
```

##Bugs

Please submit via the [GitHub issues page](https://github.com/kwongj/gene-puller/issues).  

##Software Licence

[GPLv3](https://github.com/kwongj/gene-puller/blob/master/LICENSE)
