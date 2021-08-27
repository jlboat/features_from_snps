# features_from_snps
A python script for designed to find any genomic features (genes, CDS, etc) near SNPs or arbitrary coordinates.

This program requires that an SQL database be generated from an assembly and annotation using https://github.com/daler/gffutils but is otherwise organism agnostic.

```bash
python feature_from_snps.py -h 

usage: feature_from_snps.py [-h] --input INPUT --input-type INPUT_TYPE --output OUTPUT [--gff GFF]
                            [--distance DISTANCE] [--fdr FDR] [--info INFO]
                            [--feature-type FEATURE_TYPE]

This script was originally designed to find any genes near SNPs found to be significant in a GWAS analysis performed in GAPIT but has since been generalize for any coordinates.

optional arguments:
  -h, --help            show this help message and exit
  --gff GFF             The GFF3 file for the genome or a gffutils database.
  --distance DISTANCE   The distance in kb from a SNP to search for genes (default: 10).
  --fdr FDR             The significance threshold for significant SNPs. Only applicable with
                        input-type gapit.
  --info INFO           The gene info file as obtained from Phytozome (TSV).
  --feature-type FEATURE_TYPE
                        The type of features desired from the annotation (i.e. gene,CDS).

required arguments:
  --input INPUT         The name of the input file (CSV).
  --input-type INPUT_TYPE
                        The type of input file: gapit, coordinate, range. Note: headers are
                        \*\*required\*\* - or first line lost. \*gapit - the unedited GAPIT output with
                        header\* \*coordinate - markerName,chromosome,position \*\*\* \*range -
                        markerName,chromosome,position,start,end\*
  --output OUTPUT       The output file to be created (TSV).
```

