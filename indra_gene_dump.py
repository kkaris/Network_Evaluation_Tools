from indra.tools import assemble_corpus
import sys
import csv
import argparse as ap


def main(args):
    # This file takes about 32 GB to load
    if not args.infile:
        args.infile = './Data/indra_raw/bioexp_all_raw.pkl'
    if not args.outfile:
        args.outfile = './filtered_indra_network.sif'

    # Load statements from file
    stmts_raw = assemble_corpus.load_statements(args.infile)

    # Expand families, fix grounding errors and run run preassembly
    stmts_fixed = assemble_corpus.run_preassembly(
                    assemble_corpus.map_grounding(
                        assemble_corpus.expand_families(stmts_raw)))

    # Default filtering: specific (unique) genes that are grounded.
    stmts_filtered = assemble_corpus.filter_grounded_only(
                         assemble_corpus.filter_genes_only(stmts_fixed, specific_only=True))
    # Custom filters
    if args.human_only:
        stmts_filtered = assemble_corpus.filter_human_only(stmts_filtered)
    if args.filter_direct:
        stmts_filtered = assemble_corpus.filter_direct(stmts_filtered)

    binary_stmts = [s for s in stmts_filtered if len(s.agent_list()) == 2 and s.agent_list()[0] is not None]
    rows = []
    for s in binary_stmts:
        rows.append([ag.name for ag in s.agent_list()])

    # Write rows to .sif file
    with open(args.outfile, 'w', newline='') as csvfile:
        wrtr = csv.writer(csvfile, delimiter='\t')
        for row in rows:
            wrtr.writerow(row)


if __name__ == '__main__':
    try:
        from indra.tools import assemble_corpus
    except ImportError as e:
        print('Indra not installed! Maybe switch venv?')
        sys.exit(e)

    parser = ap.ArgumentParser()
    parser.add_argument('-i', '--infile', required=True, help='Input file with indra statements')
    parser.add_argument('-o', '--outfile', required=True, help='Output file with gene interactions')
    parser.add_argument('-ho', '--human-only', action='store_true')
    parser.add_argument('-fd', '--filter-direct', action='store_true')
    args = parser.parse_args()

    main(args)

