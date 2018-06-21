import pandas as pd
import argparse as ap


def main(args):
    performance_file = args.performance_file
    perf_gain_file = args.perf_gain_file

    performance = pd.read_csv(performance_file, index_col=0, sep='\t')
    perf_gain = pd.read_csv(perf_gain_file, sep='\t')

    # List of disease gene sets (447)
    dsets = performance.index.values
    dsets_gain = perf_gain.index.values

    # Get average AUPRC score over all 447 disease sets
    performance.mean().sort_values(ascending=False).to_csv('avg_AUPRC_perf_score.csv', sep='\t')
    perf_gain.mean().sort_values(ascending=False).to_csv('avg_AUPRC_perf_gain_score.csv', sep='\t')

    # Sort after what indra recovered best and give indra's ranking for each disease set: 'disease set', 'indra's rank out of 22 sets' (resulting in 447 rows of ranking)
    indra_rank = []
    indra_rank_gain = []

    for dset in dsets:
        rnk = performance[performance.index == dset].T.sort_values(dset, ascending=False).T.columns.get_loc('INDRA Human-Grounded-Specific_only')
        indra_rank.append((dset, rnk + 1))

    for dset_gain in dsets_gain:
        rnkg = perf_gain[perf_gain.index == dset_gain].T.sort_values(dset_gain, ascending=False).T.columns.get_loc('INDRA Human-Grounded-Specific_only')
        indra_rank_gain.append((dset_gain, rnkg + 1))

    indra_rank_df = pd.DataFrame(indra_rank)
    indra_rank_df.rename(columns={0: 'Disease Set', 1:'INDRA ranking'}, inplace=True)
    indra_rank_df.sort_values('INDRA ranking').to_csv('indra_ranking_per_diseaseset.csv', sep='\t', index=False)

    indra_rank_gain_df = pd.DataFrame(indra_rank_gain)
    indra_rank_gain_df.rename(columns={0: 'Disease Set', 1: 'INDRA ranking'}, inplace=True)
    indra_rank_gain_df.sort_values('INDRA ranking').to_csv('indra_ranking_per_diseaseset_perfgain.csv', sep='\t', index=False)


if __name__ == '__main__':
    parser = ap.ArgumentParser()
    parser.add_argument('-p', '--performance_file', required=True, help='AUPRC perforamance/Z score stored as tab-csv file.')
    parser.add_argument('-g', '--perf_gain_file', required=True, help='AUPRC performance gain stored as tab separated csv file.')
    parser.add_argument('-v', '--verbose', help='More cowbell! (Default: True)', action='store_true')
    args = parser.parse_args()
    main(args)
