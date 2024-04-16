import pandas as pd
from pathlib import Path
import re

PREFIX = '/group/flchedingrp/mitochi/Work/Project/Ethan/240325/5_Analysis/2_deleted/67_72_PFC8TACINIT_240406/'


LOW_DIR = '/group/flchedingrp/mitochi/Work/Project/Ethan/240325/5_Analysis/2_deleted/67_72_PFC8TACINIT_240406/PBEH2_240409_L1_t45_w15_l50/PEAKS_GENOME'
HIGH_DIR = '/group/flchedingrp/mitochi/Work/Project/Ethan/240325/5_Analysis/2_deleted/67_72_PFC8TACINIT_240406/PBEH2_240409_L1_t55_w20_l100/PEAKS_GENOME'


MATCH_LABELS = [
    'sample_id',
    'description',
    'strand',
    'window',
    'percent',
    'convert_type',
    'file_type'
]

OUTPUT_DIR = '../resources/mergedFiles'
SAMPLE_DF_PATH = 'samples.tsv'

if not Path(OUTPUT_DIR).is_dir():
    Path(OUTPUT_DIR).mkdir()


def iterate_dir_for_peak_coord_files(dir):
    return [bed_file for bed_file in Path(dir).iterdir() if bed_file.suffix == '.bed']


def create_merge_df(bed_file):

    # Read the bed file and apply column names for easier merging later
    bed_df = pd.read_csv(str(bed_file), sep='\t', header=None)
    bed_df.columns = [
        'sample', 'start', 'end', 'PEAK', 'status', 'strand', 'convert_path'
        ]
    
    # extract the name of the conversion file (should be the same for all peaks)
    # stored in the bed file
    # NOTE: The filepaths are not absolute for all files 
    # need to change the strategy here. Going back to looking seperately
    convert_file = set(bed_df.convert_path).pop()

    # check for relative file path (some have absolute some are relative)
    convert_split = convert_file.split('/')[0]

    if convert_split == 'PBEH2_240409_L1_t55_w20_l100' or convert_split == 'PBEH2_240409_L1_t45_w15_l50':
        convert_file = str(Path(PREFIX).joinpath(convert_file))

    # Read the convert file into another dataframe
    convert_df = pd.read_csv(convert_file, sep='\t')
    convert_df = convert_df.reset_index()
    merge_df = bed_df.merge(convert_df, left_on='PEAK', right_on='index')

    return merge_df


def extract_sample_info_from_bed_name(bed_file):
    pull = r'PBEH2_bismark_bt2.bam_genePBEH2_BCBC(\d+)_PLASMID(.+)_(Pos|Neg)_(\d.)_(\d.\d+)_(.+)\.(NOPK|PEAK).genome.bed'
    match = re.search(pull, str(bed_file))
    if match:
        return dict(zip(MATCH_LABELS, match.groups()))


def generate_merge_df_filepath(output_dir, bed_file):

    return Path(output_dir).joinpath(
        Path(bed_file).with_suffix('.merge').name
    )



def main():

    call_types = (
        ('high', HIGH_DIR),
        ('low', LOW_DIR)
    )

    sample_list = []

    for each_type, each_dir in call_types:
        # get the bed files
        bed_files = iterate_dir_for_peak_coord_files(each_dir)

        for each_bed in bed_files:
            
            # Generate and write the merged peak / coordinate data to a new
            # file
            merge_df = create_merge_df(each_bed)
            merge_path = generate_merge_df_filepath(OUTPUT_DIR, each_bed)
            merge_df.to_csv(str(merge_path), sep='\t', index=False)

            sample_info = extract_sample_info_from_bed_name(each_bed)

            sample_dict = {
                'merge_path': str(merge_path),
                'call_type': each_type
            }

            sample_dict.update(sample_info)

            sample_list.append(sample_dict)
    
    sample_df = pd.DataFrame(sample_list)
    sample_df.sort_values(['sample_id', 'file_type', 'window', 'percent'])
    sample_df.to_csv(SAMPLE_DF_PATH, sep='\t', index=False)


if __name__ == '__main__':
    main()
    














