# Not run as part of the actual workflow but prior in order to create
# a speadsheet that defines samples and filepaths to use for the actual workflow
# run. 

import pandas as pd
import re
from pathlib import Path


HIGH = '/group/flchedingrp/mitochi/Work/Project/Ethan/240325/5_Analysis/2_deleted/67_72_PFC8TACINIT_240406/PBEH2_240409_L1_t55_w20_l100/.CALL'



TARGET_DIR = Path('/group/flchedingrp/mitochi/Work/Project/Ethan/240325/5_Analysis/2_deleted/67_72_PFC8TACINIT_240406/PBEH2_240409_L1_t45_w15_l50/.CALL')
CONVERT_FILE_RE = r'PBEH2_bismark_bt2.bam_genePBEH2_BCBC(\d+)_PLASMID(.+)_(Pos|Neg)_(\d.)_(\d.\d+)_(.+)\.(NOPK|PEAK)$'


# PBEH2_bismark_bt2.bam_genePBEH2_BCBC00_PLASMIDPFC9_NTBSPQI_13_DESCLINEARBSAI_NT_COSSB_Neg_15_0.45_CG.NOPK.out
match_labels = [
    'sample_id',
    'description',
    'strand',
    'window',
    'percent',
    'convert_type',
    'file_type'
]


match_list = []

for each_file in TARGET_DIR.iterdir():
    match = re.search(CONVERT_FILE_RE, str(each_file))
    if match:
        match_list.append(
            dict(zip(match_labels, match.groups()))
        )
        match_list[-1].update(
            {'filename': each_file.name}
        )

df = pd.DataFrame(match_list)
df = df.sort_values(['sample_id', 'file_type', 'window', 'percent'])

print(df.head())
df.to_csv('conversionFiles.tsv', sep='\t', index=False)