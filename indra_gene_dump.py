from indra.tools import assemble_corpus
from indra.assemblers import graph_assembler as ga
import csv

# This file takes about 32 GB to load
temp_path = '~/ASDF/data/bioexp_all_raw.pkl'
path = './data/bioexp_all_raw.pkl'

stmts_raw = assemble_corpus.load_statements(temp_path)

# We need unique, grounded, human genes only
stmts_filtered = assemble_corpus.filter_human_only(
                     assemble_corpus.filter_grounded_only(
                         assemble_corpus.filter_genes_only(stmts_raw, specific_only=True)))

binary_stmts = [s for s in stmts_filtered if len(s.agent_list()) == 2 and s.agent_list()[0] is not None]
rows = []
for s in binary_stmts:
    rows.append([ag.name for ag in s.agent_list()])

# Write rows to .sif file
with open('grounded_bioexp.sif', 'w', newline='') as csvfile:
    wrtr = csv.writer(csvfile, delimiter='\t')
    for row in rows:
        wrtr.writerow(row)

