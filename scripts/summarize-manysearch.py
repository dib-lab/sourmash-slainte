#! /usr/bin/env python
"""
Combine multiple sourmash containment reports for isolates x many metagenomes.
"""
import sys
import csv
import argparse
import pandas as pd


def main():
    p = argparse.ArgumentParser()
    p.add_argument('manysearch_csv')
    p.add_argument('-o', '--output', required=True)
    args = p.parse_args()

    df = pd.read_csv(args.manysearch_csv)
    print(df.head())

    samples = list(sorted(set(df['match_name'])))
    print(f'samples: {samples}')

    queries = list(sorted(set(df['query_name'])))

    sample_dfs = {}
    for sample in samples:
        sample_df = df[df['match_name'] == sample]
        sample_df = sample_df[['query_name', 'containment']]
        sample_dfs[sample] = sample_df

    # go through and merge, retaining the column named `containment`;
    # rename the containment column to cont_{sample}.

    # do the first sample:
    sample = samples.pop(0)
    combined_df = sample_dfs[sample]
    combined_df.rename(columns={'containment': sample}, inplace=True)

    # and then... the rest!
    while samples:
        sample = samples.pop(0)
        sample_df = sample_dfs[sample]
        sample_df = sample_df[['query_name', 'containment']]
        sample_df.rename(columns={'containment': sample},
                         inplace=True)
        combined_df = combined_df.merge(sample_df, on='query_name',
                                        how='outer')

    combined_df.rename(columns={'query_name': 'genome'}, inplace=True)
    combined_df.fillna(value=0, inplace=True)

    print(combined_df.head())

    # get the list of data columns
    cols = combined_df.columns.tolist()
    cols.remove('genome')
    data_cols = list(cols)

    # add in a max column
    combined_df['max'] = combined_df[data_cols].max(axis=1)

    # sort!
    combined_df.sort_values(by='max', inplace=True, ascending=False)

    # save!
    combined_df.to_csv(args.output, index=False)


if __name__ == '__main__':
    sys.exit(main())
