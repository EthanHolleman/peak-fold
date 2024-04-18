import pandas as pd



def main():

    comp_df = pd.read_csv(
        snakemake.input['comp_df'], sep='\t'
    )
    sample_df = pd.read_csv(
        snakemake.input['sample_df'], sep='\t'
    )

    comp_df['sample_id'] = comp_df.apply(
        lambda row: extract_sample_id(row),
        axis=1
    )

    merge = comp_df.merge(sample_df, left_on='sample_id', right_on='Sample_ID')

    merge.to_csv(snakemake.output['out'], sep='\t', index=False)



def extract_sample_id(composite_call_row):
    """Composite call tsv have sample description column but not sample ID alone
    This function pulls out that integer for easier merging with sample
    description table

    Args:
        composite_call_row (Series): Row from composite call df
    """

    return int(composite_call_row.split('_')[1].split('bcBC')[-1])



if __name__ == '__main__':
    main()

