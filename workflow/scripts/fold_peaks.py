import pandas as pd
from seqfold import fold, dg, dg_cache, dot_bracket


def main():

    peak_df = pd.read_csv(snakemake.input['peaks'], sep='\t')
    peak_df['fold'] = peak_df.apply(
        lambda row: add_rna_folding_column(row),
        axis=1
    )
    peak_df.to_csv(snakemake.output['out'], sep='\t', index=False)



def add_rna_folding_column(row):

    #print(row)

    def get_peak_seq():
        return row.sequence[row.start: row.end]

    def fold_peak(peak_seq):
        structs = fold(peak_seq)
        dot = dot_bracket(peak_seq, structs)
        return dot

    peak_seq = get_peak_seq()
    return fold_peak(peak_seq)



if __name__ == '__main__':
    main()