#! /usr/bin/env python
"""
Do a clustered plot of isolates vs sweep metagenomes
"""
import argparse
import sys
import pandas as pd
import seaborn as sns


def main():
    p = argparse.ArgumentParser()
    p.add_argument('summary_csv', help="output of 'summarize-manysearch.py'")
    p.add_argument('-o', '--output', help="figure name", required=True)
    args = p.parse_args()

    df = pd.read_csv(args.summary_csv)
    df.drop(['max'], inplace=True, axis=1)

    df = df.set_index('genome')
    df = df.transpose()

    cluster_plot = sns.clustermap(df, xticklabels=True, yticklabels=True,
                                  figsize=(18, 8))
    cluster_plot.savefig(args.output)


if __name__ == '__main__':
    sys.exit(main())
