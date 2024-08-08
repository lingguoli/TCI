# Tissue Contribution Index

### samples
A BED format file, each line at least consists of 3 tab separated columns: 
1. The chromosome name
2. The start coordinate
3. The end coordinate

An example file is supplied, *demo.bed*.
---

### Usage

```
Version:   v1.0.0
Usage:     TCI -i <BED file>  [options]

    -i     Input your BED file

    -o     Output dir.
           -Default: pwd

    -min   Minimum of DNA fragment lenth.
           -Default: 150.

    -max   Maximum of DNA fragment lenth.
           -Default: 210.

    -T     Specify your own TSS bed file.
           This file contains at least four columns,
           eg. "chr1  100000  102000  TSS_ID_00001".
           -Default: "./supplemental/hg38.tss.rename.UD1000.bed"

    -R     Specify your own tissue specific highly expressed gene table.
           This file contains at least two columns named
           "Specific high expression tissue" and "TSS ID".
           -Default: "./supplemental/TPM_atlas.csv"

    -S     If using single-end sequencing BED, please enable this option.
```

---
### Example
```
TCI -i demo.bed 
```

