 
# Conversion tsv encoding scheme

```
 0 is bad or not data (6)
 1 is non C 
 4 is CH non conv
 5 is CG non conv
 6 is CH conv
 7 is CG conv
 8 is CH peak
 9 is CG peak
 ```

if the file looks at "CH", then all G will be 1 

if the file is "GH" then all C will be 1

Cytosine can be in 3 contexts: CpG, CHG, or CHH (H is non-G nucleotides), "CH" means the CHG or CHH


# Threshold

low threshold: `/group/flchedingrp/mitochi/Work/Project/Ethan/240325/5_Analysis/2_deleted/67_72_PFC8TACINIT_240406/PBEH2_240409_L1_t45_w15_l50/.CALL`

- Conversion TSV File: `<file>.PEAK.out` or `<file>.NOPK.out`. 1st column is read name, 2nd column onwards are the conversion

- Their nucleotide can be found in their corresponding `<file>.PEAK` or `<file>.NOPK`, first line column 6 onwards.

Most important groupings from the name are strand, convtype, and peak/nopeak, e.g.:

`PBEH2_bismark_bt2.bam_genePBEH2_BCBC<BC>_PLASMID<PLASMID>_DESC<DESCRIPTION>_<strand>_<window>_<threshold>_<convtype>.<PEAK/NOPK>.out`

- strand: `(Pos|Neg)` -> this means the read is "positive" (same as provided fasta in fasta/amplicon) or "negative" strand (reverse of the fasta in fasta/amplicon). 
- convtype: `(CH|GH|CG|GC)` -> this looks at CH/GH (excluding C/G of CpG) or CG/GC (include C/G of CpG)
- peak or no peak `(PEAK|NOPK)`

For example: `PBEH2_bismark_bt2.bam_genePBEH2_BCBC43_PLASMIDT7INIT_VR20_DESCLINEARAPALI_TX_CORPA_DEPROTEINZED_Pos_15_0.45_CH.PEAK.out`
strand is `Pos`, convtype is `CH`, has PEAK

convtype example:

```
Sequence: ACCCAGCGCA
Pos*CH:   -CCC----C-
Neg*GH:   -----G----
Pos*CG:   -CCC--C-C-
Neg*GC:   -----G-G--

'-' means masked
```

# File labeling

There are 16 files of each unique `BC+gene+description` coz there are 2 strands (Pos|Neg), 4 convtype (CH|GH|CG|GC), and whether they have peak or not (PEAK|NOPK) => (2x4x2=16)
 
Note that the nucleotide found in the `<file>.PEAK` or `<file>.NOPK` file will be exactly the same for these 16 files. So for example these 4 will have the same nucleotides:
- `PBEH2_bismark_bt2.bam_genePBEH2_BCBC43_PLASMIDT7INIT_VR20_DESCLINEARAPALI_TX_CORPA_DEPROTEINZED_Pos_15_0.45_CH.PEAK`
- `PBEH2_bismark_bt2.bam_genePBEH2_BCBC43_PLASMIDT7INIT_VR20_DESCLINEARAPALI_TX_CORPA_DEPROTEINZED_Neg_15_0.45_GH.PEAK`
- `PBEH2_bismark_bt2.bam_genePBEH2_BCBC43_PLASMIDT7INIT_VR20_DESCLINEARAPALI_TX_CORPA_DEPROTEINZED_Pos_15_0.45_CG.PEAK`
- `PBEH2_bismark_bt2.bam_genePBEH2_BCBC43_PLASMIDT7INIT_VR20_DESCLINEARAPALI_TX_CORPA_DEPROTEINZED_Neg_15_0.45_GC.PEAK`