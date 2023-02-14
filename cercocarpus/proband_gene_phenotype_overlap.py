#!/usr/bin/env python

import re
import sys
import json
import math
import argparse
import numpy as np
import pandas as pd
import seaborn as sns
import networkx as nx
import matplotlib.pyplot as plt
from matplotlib_venn import venn2
from scipy.stats import hypergeom
from networkx.drawing.nx_agraph import write_dot, graphviz_layout
from tqdm import tqdm
from pprint import pprint
from joblib import Parallel, delayed

description_text = (
    """
    Synopsis:
    
    pgpo_template.py --gene CHD7,CFTR --proband_terms hpt_terms_proband.tsv \
                     --phen2gene phenotype_to_genes.txt --json hp.json
    
    Description:
    
    Calculate the shared information content betweeen the HPO terms
    associated with a proband and a gene(s).
    
    Prints to STDOUT the details of the shared information content of the proband and
    gene(s) HPO phenotype overlap in tab-delimited format.
    """
)

def main():
    parser = argparse.ArgumentParser(
        description=description_text,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--gene', '-g',
                        help='A comma-separated list of candidate genes.')
    parser.add_argument('--proband_terms', '-t', dest='proband_terms_file',
                        help='A file containing the list of HPO IDs (first column) associated with the proband.')
    parser.add_argument('--phen2gene', '-p', dest='phen2gene_file',
                        help='The phenotype_to_genes.txt file from HPO')
    parser.add_argument('--json', '-j', dest='json_file',
                        help='The hp.json file from HPO')
    parser.add_argument('--jobs', '-n', type=int, default=2,
                        help='The number of jobs to run in parallel')
    # parser.add_argument('--argument', '-a',
    #                     help='argument is...')
    # parser.add_argument('--argument', '-a',
    #                     help='argument is...')
    args = parser.parse_args()

    sys.stderr.write('INFO : loading_data_file : ' + args.proband_terms_file)

    # Parse Text Files to Datafames
    gene_txt = args.gene
    gene_list = gene_txt.split(',')
    df_cnd = pd.DataFrame(gene_list, columns=['gene'])
    df_prb = pd.read_table(args.proband_terms_file)
    df_prb.drop_duplicates(subset='id', inplace=True)
    df_prb.set_index('id', inplace=True)
    prb_ids = set(df_prb.index.to_list())
    
    # Parse HPO Phenotype_to_Gene File to Dataframe
    df_p2g = pd.read_table(args.phen2gene_file,
                           skiprows=1,
                           header=None,
                           index_col=False,
                           names=['id', 'term', 'mim', 'gene', 'source', 'source_id'])
    df_p2g.drop_duplicates(subset=['id', 'gene'], inplace=True)
    df_p2g = df_p2g.set_index('id')
    
    # # Parse HPO Genes_to_Phenotype File to Dataframe
    # df_g2p = pd.read_table(g2p,
    #                        skiprows=1,
    #                        header=None,
    #                        index_col=False,
    #                        names=['gene_id', 'gene', 'id', 'term', 'freq_raw', 'freq_hpo', 'info_gd_source',
    #                               'gd_source', 'disease_id'])
    # df_g2p.drop_duplicates(subset=['id', 'gene'], inplace=True)
    # df_g2p = df_g2p.set_index('id')
    
    # # Parse HPO Phenotpye File to Dataframe
    # 
    # df_phen = pd.read_table(phen,
    #                        skiprows=5,
    #                        header=None,
    #                        index_col=False,
    #                        names=['disease_id', 'disease', 'qualifier', 'id', 'reference', 'evidence', 'onset',
    #                               'frequency', 'sex', 'modifier', 'aspect', 'biocuration'])
    # df_phen.drop_duplicates(subset=['id', 'disease_id'], inplace=True)
    # df_phen = df_phen.set_index('id')
    
    # Get term frequency and information content (ic)
    term_counts = pd.DataFrame(df_p2g.reset_index()['id'].value_counts())
    term_counts.rename(columns={'id': 'count'}, inplace=True)
    term_counts.index.rename('id', inplace=True)
    term_counts['freq'] = term_counts['count'] / term_counts['count'].sum()
    term_counts['ic'] = -1 * np.log10(term_counts['freq'])
    term_counts.loc['HP:0000118'] = [term_counts['count'].sum(), 1, 0]
    df_p2g = df_p2g.join(term_counts)
    #df_g2p = df_g2p.join(term_counts)
    
    # Parse HPO JSON
    with open(args.json_file) as json_fh:
        data = json.load(json_fh)
        
    # Parse Nodes from HPO JSON Data
    nodes = list()
    for node in data['graphs'][0]['nodes']:
        hp_id = node['id']
        hp_id = re.sub('.*\/', '', hp_id)
        hp_id = re.sub('_', ':', hp_id)
        if not re.match('^HP:', hp_id):
            continue
        hp_lbl = node['lbl']
        hp_def = 'No definition'
        if 'meta' in node.keys():
            if 'definition' in node['meta'].keys():
                hp_def = node['meta']['definition']['val']
                nodes.append((hp_id, {'term': hp_lbl, 'definition': hp_def}))
                
    # Parse Edges from HPO JSON Data
    edges = []
    for edge in data['graphs'][0]['edges']:
        sub = edge['sub']
        sub = re.sub('.*\/', '', sub)
        sub = re.sub('_', ':', sub)
        
        obj = edge['obj']
        obj = re.sub('.*\/', '', obj)
        obj = re.sub('_', ':', obj)
        
        edges.append((obj, sub))

    # Create NetworkX Graph of HPO
    G = nx.DiGraph()
    G.add_nodes_from(nodes)
    G.add_edges_from(edges)
    
    # Trim HPO to subgraph of 'Phenotypic abnormality'
    
    root = 'HP:0000118' # Phenotypic abnormality
    sub_ids = nx.descendants(G, 'HP:0000118')
    sub_ids.add('HP:0000118')
    pabG = G.subgraph(sub_ids)

    # Get subgraph/leaves of proband HPO ancestors
    prb_id_ancestors = set()
    for id in prb_ids:
        prb_id_ancestors.add(id)
        if id not in pabG.nodes:
            # WARN
            continue
        prb_id_ancestors.update(nx.ancestors(pabG, id))
        prb_id_ancestors.add('HP:0000118')
        prbG = pabG.subgraph(prb_id_ancestors)
        prb_leaves = [node for node in prbG.nodes() if prbG.out_degree(node)==0]

    #--------------------------------------------------------------------------------

    # Get subgraph/leaves of gene HPO ancestors
    gene = df_cnd.iloc[0]['gene']
    gene_ids = set(df_p2g.query('gene == @gene').index.to_list())
    gene_id_ancestors = set()
    for id in gene_ids:
        gene_id_ancestors.add(id)
        if id not in pabG.nodes:
            # WARN
            continue
        gene_id_ancestors.update(nx.ancestors(pabG, id))
    gene_id_ancestors.add('HP:0000118')
    genG = pabG.subgraph(gene_id_ancestors)
    gen_leaves = [node for node in genG.nodes() if genG.out_degree(node)==0]

    # Get subgraph of proband & gene IDs
    prb_gene_ids = prb_id_ancestors.union(gene_id_ancestors)
    pgG = pabG.subgraph(prb_gene_ids)
    pgG_ids = set(pgG.nodes)
    pgUG = nx.Graph(pgG)
    
    gene_chunks = np.array_split(gen_leaves, 8)
    all_lcas = []
    for prb_id in tqdm(prb_leaves):
        if prb_id not in pgG_ids:
            # WARN
            continue
        
        # Check for gene_ids in genG_ids
        gene_lcas = []
        lcas = Parallel(n_jobs=args.jobs)(delayed(get_lcas)(prb_id, genes, pgG, pgUG) for genes in gene_chunks)
        for lca in lcas:
            all_lcas.extend(lca)

    df_lcas = pd.DataFrame(all_lcas, columns=['prb_id', 'gene_id', 'lca_id', 'shpl'])
    df_lcas.set_index('lca_id', inplace=True)
    df_lcas = df_lcas.join(term_counts, how='left')
    df_lcas.sort_values(by='ic', ascending=False, inplace=True)
    df_lcas.drop_duplicates(subset='prb_id', keep='first', inplace=True)

    df_lcas = df_lcas.merge(df_p2g['term'], left_on='prb_id', right_index=True).drop_duplicates()
    df_lcas.rename(columns={'term': 'prb_term'}, inplace=True)

    df_lcas = df_lcas.merge(df_p2g['term'], left_on='gene_id', right_index=True).drop_duplicates()
    df_lcas.rename(columns={'term': 'gene_term'}, inplace=True)

    df_lcas = df_lcas.merge(df_p2g['term'], left_index=True, right_index=True).drop_duplicates()
    df_lcas.rename(columns={'term': 'lca_term'}, inplace=True)

    df_lcas.sort_values(by='ic', ascending=False, inplace=True)
    df_lcas.reset_index(inplace=True, drop=False)
    df_lcas.rename(columns={'index': 'lca_id'}, inplace=True)
    df_lcas = df_lcas.loc[:,['ic', 'prb_term', 'gene_term', 'lca_term', 'prb_id', 'gene_id', 'lca_id', 'shpl', 'count', 'freq']]
    
    print(df_lcas.to_csv(sep='\t'))
            
        # prb_gen_lcas = pd.DataFrame(columns=['prb_id', 'gen_id', 'lca_id'])
        # prb_gene_lcas.append(prb_gene_lca)
        #if prb_id == gene_id:
        #    ic = 0
        #    if prb_id in term_counts.index:
        #        ic = term_counts.loc[prb_id]['ic']
        #        prb_gen_lcas.append((prb_id, ic, (prb_id, gene_id)))
        #    return
        #
        #if gene_id not in pgG_ids:
        #    # WARN
        #    return
        #
        #if lca == 'HP:0000118':
        #    ic = 0
        #elif lca in term_counts.index:
        #    ic = term_counts.loc[lca]['ic']
        #    prb_gene_lcas = sorted(prb_gen_lcas, key=lambda x: x[1], reverse=True)
        #    gene_lcas.append(prb_gene_lcas[0])

    # sys.stderr.write('DEBUG : gene_lcas : ' + str(len(gene_lcas)))
    # gene_lcas = sorted(gene_lcas, key=lambda x: x[1], reverse=True)
    # for lca in gene_lcas:
    #     (lca_id, ic, pair) = lca
    #     prb_id, gen_id = pair
    #     lca_term = G.nodes[lca_id]['term']
    #     prb_term = G.nodes[prb_id]['term']
    #     gen_term = G.nodes[gen_id]['term']
    #     print('\t'.join([gene, prb_term, gen_term, lca_term, str(ic)]))

    # # ### Write out Graphvix Dot File of Proband Ancestors Graph
    # 
    # #write_dot(prbG,'proband.dot')
    # # dot -Tpng test.dot > test.png
    # 
    # # Trim proband terms Dataframe to leaf node terms
    # df_prb = df_prb[df_prb.index.isin(prb_trim)]
    # 
    # # Expand proband terms Dataframe to include ancestors 
    # df_ancestors = pd.DataFrame(prbG.nodes, columns=['id'])
    # df_ancestors.set_index('id', inplace=True)
    # df_ancestors['review'] = 1
    # df_ancestors['context'] = 'Ancestor'
    # df_ancestors = df_ancestors.join(df_p2g['term'].drop_duplicates(), how='left')
    # 
    # # Score All HPO Disease Genes Against Proband Terms with Hypergeometric Test
    # # Expect this step to take ~6 minutes
    # data = []
    # all_genes = set(df_p2g[df_p2g.index.isin(G.nodes())]['gene'].to_list())
    # for gene in tqdm(all_genes):    
    #     df_gene = df_p2g.query('gene == @gene')
    #     # df_inter = df_gene.join(df_prb, how='inner', lsuffix='_gene', rsuffix='_proband').reset_index()
    #     df_inter = df_gene.join(df_ancestors, how='inner', lsuffix='_gene', rsuffix='_proband').reset_index()
    #     df_inter.drop_duplicates(subset='id', inplace=True, ignore_index=True)
    # 
    #     x = len(df_inter.index.to_list())    # Intersection
    #     M = len(set(df_p2g.index.to_list())) # Total number of terms
    #     # n = len(df_prb.index.to_list())      # Proband terms
    #     n = len(df_ancestors.index.to_list())      # Proband ancestor terms
    #     N = len(df_gene.index.to_list())     # Gene terms
    # 
    # 
    #     pval = hypergeom.sf(x - 1, M, n, N)
    #     phred = -1 * math.log10(pval)
    #     
    #     data.append({'gene': gene, 'x': x, 'n': n, 'N': N, 'pval': pval, 'phred': phred})
    # 
    # df_data = pd.DataFrame(data)
    # 
    # # Summary statistics for hypergeometric test input/results
    # print(df_data.describe())
    # 
    # # Plot distribution of genes/proband 'term overlap' and -log10(pval) from hypergeometric test
    # sns.kdeplot(df_data['x']);
    # sns.kdeplot(df_data['phred']);
    # 
    # # Add Percentile Rank Based on -log10(p) of Hypergeometric Test
    # df_data['pct_rank'] = df_data['phred'].rank(pct = True)
    # 
    # # Basic Metrics on Dataset 
    # print("Candidate Genes:        " + str(len(df_cnd)))
    # print("Proband Terms:          " + str(len(df_prb)))
    # print("Proband Ancestor Terms: " + str(len(df_ancestors)))
    # print("Proband Ancestor Graph: " + str(len(prbG)))
    # 
    # # Table of hypergeometric details for candidate gene
    # df_cnd_mrg = df_data.merge(df_cnd, on='gene')
    # df_cnd_mrg.rename(columns={'x': 'overlap', 'n': 'prb_terms', 'N': 'gene_terms'}, inplace=True)
    # df_cnd_mrg.sort_values(by='pval', ascending=True, inplace=True)
    # 
    # # Display details for inidividual candidate genes
    # all_overlap = set()
    # cnd_gene_count = int(len(df_cnd))
    # for gene in df_cnd_mrg['gene'].to_list():
    #     df_gene = df_p2g.query('gene == @gene')
    #     df_inter = df_gene.join(df_ancestors, how='inner', lsuffix='_gene', rsuffix='_proband').reset_index()
    #     df_inter.drop_duplicates(subset='id', inplace=True, ignore_index=True)
    #     all_overlap |= set(df_inter['term_gene'])
    #     
    #     x = len(df_inter.index.to_list())     # Intersection
    #     M = len(set(df_p2g.index.to_list()))  # Total number of terms
    #     n = len(df_ancestors.index.to_list()) # Proband ancestor terms
    #     N = len(df_gene.index.to_list())      # Gene terms
    #     rank = round((df_data[df_data['gene'] == gene]['pct_rank'].iloc[0] * 100), 1)
    #     
    #     print('Gene           = ' + gene)
    #     print('Total terms    = ' + str(M))
    #     print('Proband terms  = ' + str(n))
    #     print('Gene terms     = ' + str(N))
    #     print('Shared terms   = ' + str(x))
    #     
    #     pval = hypergeom.sf(x - 1, M, n, N)
    #     pv_adj = pval * float(cnd_gene_count)
    #     phred = -1 * math.log10(pv_adj)
    #     print('pval           = ' + '%s' % float('%.2g' % pval))
    #     print('pv_adj         = ' + '%s' % float('%.2g' % pv_adj))
    #     print('')
    #     print('-log10(pv_adj) = ' + str(round(phred, 1)))
    #     print('Prcntl rank    = ' + str(rank))
    # 
    #     venn2([set(df_ancestors.index.to_list()),
    #            set(df_gene.index.to_list())],
    #           set_labels = ('Proband', gene))
    #     plt.show()
    # 
    #     print(df_inter[['id', 'term_gene']].sort_values(by='id'))
    # 
    #     print('\n\n--------------------------------------------------------------------\n\n')
    # 
    # print("ALL OVERLAPS\n------------\n")
    # all_overlap = list(all_overlap)
    # all_overlap.sort()
    # print('\n'.join(all_overlap))
    # 
    # # Show Hypergeometric Details for Top-10 Proband/Gene Overlaps
    # print(df_data.sort_values(by='pval'))
    # 
    # # Display details for individual genes in top-10 proband/gene overlaps
    # all_overlap = set()
    # top_genes = df_data.sort_values(by='pval').head(10)
    # top_genes = top_genes['gene'].to_list()
    # 
    # for gene in top_genes:
    #     df_gene = df_p2g.query('gene == @gene')
    #     df_inter = df_gene.join(df_ancestors, how='inner', lsuffix='_gene', rsuffix='_proband').reset_index()
    #     df_inter.drop_duplicates(subset='id', inplace=True, ignore_index=True)
    #     all_overlap |= set(df_inter['term_gene'])
    #     
    #     x = len(df_inter.index.to_list())     # Intersection
    #     M = len(set(df_p2g.index.to_list()))  # Total number of terms
    #     n = len(df_ancestors.index.to_list()) # Proband terms
    #     N = len(df_gene.index.to_list())      # Gene terms
    #     rank = round(df_data[df_data['gene'] == gene]['pct_rank'].iloc[0] * 100)
    #     
    #     print('Gene = ' + gene)
    #     print('Total terms   = ' + str(M))
    #     print('Proband terms = ' + str(n))
    #     print('Gene terms    = ' + str(N))
    #     print('Shared terms  = ' + str(x))
    #     
    #     pval = hypergeom.sf(x - 1, M, n, N)
    #     pv_adj = pval * cnd_gene_count
    #     phred = -1 * math.log10(pv_adj)
    #     print('pval           = ' + '%s' % float('%.2g' % pval))
    #     print('pv_adj         = ' + '%s' % float('%.2g' % pv_adj))
    #     print('')
    #     print('-log10(pv_adj) = ' + str(round(phred, 1)))
    #     print('Prcntl rank    = ' + str(rank))
    #     
    #     venn2([set(df_ancestors.index.to_list()),
    #            set(df_gene.index.to_list())],
    #           set_labels = ('Proband', gene))
    #     plt.show()
    #     
    #     print(df_inter[['id', 'term_gene']])
    #     
    #     print('\n\n--------------------------------------------------------------------\n\n')
    # 
    # print("ALL OVERLAPS\n------------\n")
    # all_overlap = list(all_overlap)
    # all_overlap.sort()
    # print('\n'.join(all_overlap))

def get_lcas(prb_id, genes, pgG, pgUG):

    lcas = []
    for gene_id in genes:
        lca = nx.lowest_common_ancestor(pgG, gene_id, prb_id)
        shpl = nx.shortest_path_length(pgUG, gene_id, prb_id)
        lcas.append((prb_id, gene_id, lca, shpl))
    return lcas
    
if __name__ == "__main__":
    main()

