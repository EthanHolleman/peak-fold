import pandas as pd
import re

# STATS TO CALCULATE PER PEAK BASIS
# - percent paired vs unpaired
# - percent converted vs unconverted
# - unconverted bases in paired vs shuffle mean
# peak length

CONVERTED = set([8, 9])
UNCONVERTED = set([4, 5])


def main():

    peak_df = pd.read_csv(snakemake.input['comp_df'], sep='\t')
    peak_df['peak_call'] = peak_df.apply(
        lambda row: make_peak_call(row), axis=1
    )
    stat_df = pd.DataFrame(
        peak_df.apply(
            lambda row: calc_percent_stats(row), axis=1
        ).to_list()
    )

    peak_df['converted'] = stat_df['converted']
    peak_df['unconverted'] = stat_df['unconverted']
    peak_df['paired'] = stat_df['paired']

    peak_df['plasmid'] = peak_df.apply(
        lambda row: extract_plasmid_name(row),
        axis=1
    )

    peak_df.to_csv(snakemake.output['out'], sep='\t', index=False)




def make_peak_call(row):
    """Using the composite peak call values (which get stored as a string version
    of python list so need to convert back to list of ints) pull out the region
    that is part of a SMRF-seq peak.

    Args:
        row (Series): Row from peak dataframe
    """

    def list_string_to_list(list_string):
        return [int(s) for s in list_string.replace(' ', '')[1:-1].split(',')]

    call = list_string_to_list(row.comoposite_call)
    return call[row.start:row.end]


def calc_percent_stats(row):
    """Calculate basic bisulfite conversion and DNA folding stats and return
    them as a dictionary

    Args:
        row (Series): Row from peaks dataframe

    Returns:
        dict: Dictionary of statistics 
    """
    # Bisulfite converted vs unconverted basic stats

    total_convert = sum([
        1 for each_call in row.peak_call if each_call in CONVERTED
        ])
    total_unconverted = sum([
        1 for each_call in row.peak_call if each_call in UNCONVERTED
        ])
    
    total_convert_opps = total_convert + total_unconverted

    convert_percent = total_convert / total_convert_opps
    unconvert_percent = total_unconverted / total_convert_opps

    # RNA folding paired vs unpaired basic stats

    percent_paired = (row.fold.count('(') + row.fold.count(')')) / len(row.fold)

    return {
        'converted': convert_percent,
        'unconverted': unconvert_percent,
        'paired': percent_paired
    }


def extract_plasmid_name(row):
    """Get the name of the plasmid from the sample description

    Args:
        row (Series): Peak dataframe row

    Returns:
        str: Name of the plasmid peak is called on
    """

    plasmid = re.search(r'plasmid(.+)desc', row['sample']).groups()[0]
    return plasmid


if __name__ == '__main__':
    main()



