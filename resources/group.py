import pandas as pd

samples = pd.read_csv('/home/ethollem/workflows/rloop_fold/workflow/samples/samples.tsv', sep='\t')

grouped = samples.groupby(['call_type', 'sample_id', 'strand', 'file_type'])

dfs = [grouped.get_group(group) for group in grouped.groups]

for i, each_df in enumerate(dfs):
    each_df['group_id'] = i

with_groups = pd.concat(dfs)

with_groups.to_csv('samples.grouped.tsv', sep='\t', index=False)
