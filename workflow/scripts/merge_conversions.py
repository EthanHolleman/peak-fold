import pandas as pd

# All frames should have the same peaks called in them so use the first
# frame that we read to get the peak names


# We really just care about Cs that are called as peaks and those that
# are not. 7 or 8 means called in a peak so for each read get all of the
# different call types from the different files and make a "consensus" call
# which marks 7 or 8 if any of the calls have 7 or 8 in that position

def make_composite_call(peak_name, frames):

    peak = frames[0].loc[frames[0].PEAK_x == peak_name]

    def extract_read_seq():
        colnames = peak.iloc[:, 13:].columns.tolist()
        seq = [col.split('.')[0] for col in colnames]
        return ''.join(seq)


    def get_calls_from_frame(frame):

        # Sequence starts at column 14
        calls = peak.iloc[:, 13:]
        assert len(calls) == 1
        return calls.iloc[0].tolist()

    def get_id_info_from_frame(frame):

        info = peak.iloc[:, :13]
        return info
    
    # get all of the calls as a list
    calls = [get_calls_from_frame(each_frame) for each_frame in frames]
    info = get_id_info_from_frame(frames[0])
    seq = extract_read_seq()

    consensus_call = calls[0]
    calls = calls[1:]

    for each_call_list in calls:
        for i, each_call in enumerate(each_call_list):

            if each_call == 8 or each_call == 9:
                consensus_call[i] = each_call
    
    info['comoposite_call'] = [consensus_call]
    info['sequence'] = seq

    return info




def main():

    # files from the same sample, strand and call type
    sample_files = snakemake.input
    print(sample_files)
    frames = [pd.read_csv(sample, sep='\t') for sample in sample_files]

    print(frames)

    peaks = list(frames[0].PEAK_x)

    composite_calls = []

    for each_peak in peaks:
        try:
            c_call = make_composite_call(each_peak, frames)
            composite_calls.append(c_call)
        except AssertionError:
            continue
    
    c_calls_df = pd.concat(composite_calls)
    c_calls_df.to_csv(snakemake.output['out'], sep='\t', index=False)


if __name__ == '__main__':
    main()
    




