# Pyphillin 
#### (pronounced pie-fillin)
## Functional profiles constructed from 16S amplicon sequence data.

`Pyphillin` is a Python-clone of the functional profiling approach published by SecondGenome in _Piphillin: Improved Prediction of Metagenomic Content by Direct Inference from Human Microbiomes_ published in PLOS ONE (2016) and _Piphillin predicts metagenomic composition and dynamics from DADA2-corrected 16S rDNA sequences._ published in BMC Genomics (2020). While the authors demonstrate the benefits of their method, the company has decided to remove access to the research community. Pyphillin is an open-source, freely available clone of Piphillin. 

## Getting Started

### Installation
`Pyphillin` will be installable trhough the python package index PyPI. This will download the library, create an executable called `pyphillin`, as well as download the necessary reference databases.
```python 
pip install pyphillin
```

### Running
The current beta version of `Pyphillin` is command-line based. Updates will be coming.

Note: the current version requires the prior installation of VSEARCH [https://github.com/torognes/vsearch](https://github.com/torognes/vsearch) for performing sequence alignment. 
The input files required are:
 - taxonomic abundance table in tsv format, with rows representing samples and columns representing ASVs/OTUs
 - representative sequences for the ASVs/OTUs in the taxonomic abundance table, fasta format
 
