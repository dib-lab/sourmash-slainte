#! /usr/bin/env python
"""
Do a prefetch-style overlap analysis, but weight the overlap by metag hash
abund.
"""
import sys
import csv
import argparse
import sourmash
import numpy as np


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--genomes', nargs='+', required=True,
                   help="genome queries")
    p.add_argument('--metagenomes', nargs='+', required=True,
                   help="metagenomes to search")
    p.add_argument('-o', '--output', help="output CSV")
    p.add_argument('-k', '--ksize', default=21, type=int)
    args = p.parse_args()

    # load all the genomes!
    queries = []
    for filename in args.genomes:
        queries.extend(sourmash.load_file_as_signatures(filename,
                                                        ksize=args.ksize))

    assert all( (not q.minhash.track_abundance for q in queries) )

    columns = ['intersect_bp',
               'match_filename',
               'match_name',
               'match_md5',
               'query_filename',
               'query_name',
               'query_md5',
               'ksize',
               'moltype',
               'scaled',
               'f_unique_weighted',
               'n_unique_weighted_found',
               'sum_weighted_found',
               'total_weighted_hashes',
               'average_abund',
               'median_abund',
               'std_abund',
               ]
    out_fp = open(args.output, 'w', newline='')
    out_w = csv.DictWriter(out_fp, fieldnames=columns)
    out_w.writeheader()

    # go through metagenomes one by one
    for metag_filename in args.metagenomes:
        metag = sourmash.load_file_as_signatures(metag_filename,
                                                 ksize=args.ksize)
        metag = list(metag)
        assert len(metag) == 1
        metag = metag[0]

        # make sure metag has abundance!
        assert metag.minhash.track_abundance, filename

        # calculate total weighted hashes for use in denominator:
        total_sum_abunds = metag.minhash.sum_abundances

        # get a flattened copy for use in intersections...
        flat_metag = metag.minhash.flatten()

        # other info!

        results_template = dict(match_md5=metag.md5sum(),
                                match_name=metag.name,
                                match_filename=metag_filename,
                                ksize=args.ksize,
                                moltype=metag.minhash.moltype,
                                scaled=metag.minhash.scaled)

        # now, loop & get weighted containment for each genome
        for q in queries:
            intersect_mh = q.minhash.intersection(flat_metag)
            w_intersect_mh = intersect_mh.inflate(metag.minhash)

            abunds = list(w_intersect_mh.hashes.values())

            if abunds:
                mean = np.mean(abunds)
                median = np.median(abunds)
                std = np.std(abunds)
                overlap_sum_abunds = w_intersect_mh.sum_abundances
                f_sum_abunds = overlap_sum_abunds / total_sum_abunds
            else:
                mean = median = std = 0
                overlap_sum_abunds = 0
                f_sum_abunds = 0

            results_d = dict(intersect_bp=len(intersect_mh),
                             query_filename=q.filename,
                             query_name=q.name,
                             query_md5=q.md5sum(),
                             f_unique_weighted=f_sum_abunds,
                             n_unique_weighted_found=overlap_sum_abunds,
                             sum_weighted_found=overlap_sum_abunds,
                             total_weighted_hashes=total_sum_abunds,
                             average_abund=mean,
                             median_abund=median,
                             std_abund=std)
            results_d.update(results_template)

            out_w.writerow(results_d)

            print(len(intersect_mh), overlap_sum_abunds, total_sum_abunds,
                  f"{f_sum_abunds:.03f}")

    out_fp.close()


if __name__ == '__main__':
    sys.exit(main())
