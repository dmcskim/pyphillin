import pandas as pd
import argparse
import subprocess
import gzip
import pkg_resources

COPY_NUMBER_FILE = pkg_resources.resource_filename(__name__,\
                    'data/pyphillin_16S_copy_number.tsv.gz')
KEGG_PROFILES_FILE = pkg_resources.resource_filename(__name__,\
                    'data/pyphillin_KEGG_profiles.tsv.gz')
SEQUENCE_DB = pkg_resources.resource_filename(__name__,\
                    'data/pyphillin_16S_sequences.fasta.gz')

def align_sequences(args):
    # - align 16S sequences to reference
    cmd = ['vsearch', '--usearch_global', args.rep_fasta, '--db',\
           SEQUENCE_DB, '--id', str(args.pct_id),\
           '--top_hits_only', '--maxaccepts', '0', '--maxrejects', '0',\
           '--uc_allhits', '--blast6out', args.blast_out]
    #print(cmd)
    subprocess.call(cmd)
    return

def profile_functions(args):
    # needs bug abundances per sample, K and 16S abundances per bug
    #   replace sequences with functional profiles using formula
    tax_table = pd.read_csv(args.tax_table, sep='\t', index_col=0)
    total_read_cnt = tax_table.sum().sum()
    #print('tax table loaded')
    copy_num = pd.read_csv(COPY_NUMBER_FILE, index_col=0, sep='\t')
    #print('16S copy numbers loaded')
    kegg_profiles = pd.read_csv(KEGG_PROFILES_FILE, index_col=0, sep='\t')
    #print('KEGG profiles loaded')
    blast_res = pd.read_csv(args.blast_out, index_col=0, sep='\t', header=None)
    blast_res[1] = [x.split('|')[0] for x in blast_res[1]]
    #print('blast results loaded')

    # get unique mapping of ASV -> reference
    asv_map = {}
    asv_cnt = {}
    all_cnt = 0
    for tasv in blast_res.index.unique():
        ttemp = blast_res.loc[tasv,1]
        if type(ttemp) == str:
            asv_map[tasv] = [ttemp]
        else:
            asv_map[tasv] = blast_res.loc[tasv,1].unique().tolist()
        asv_cnt[tasv] = len(asv_map[tasv])
        all_cnt += len(asv_map[tasv])

    # only keep matched results
    keep = [x for x in tax_table.columns if x in asv_map.keys()]
    ndata = tax_table[keep]
    aligned_read_cnt = ndata.sum().sum()

    # divide ASV read count by number of top vsearch hits (asv_cnt)
    tcnt = [asv_cnt[x] for x in ndata.columns]
    ndata = ndata.div(tcnt)

    #print(aligned_read_cnt, total_read_cnt, aligned_read_cnt/total_read_cnt)
    print('Reads used: {0:0.2%}'.format(aligned_read_cnt/total_read_cnt))

    # rename to refseq id
    # ** Make empty table using all_cnt (# of columns needed)
    # ** For each ASV in ndata.columns, make duplicate of column and apply
    #    asv_map column labels
    #tcols = [asv_map[x] for x in ndata.columns]
    #ndata.columns = tcols
    tdata = pd.DataFrame()
    tcols = []
    ti = 0
    for asv,taxa in asv_map.items():
        for tax in taxa:
            tdata[ti] = ndata[asv]
            tcols.append(tax)
            ti += 1
    tdata.columns = tcols
    tdata.index = ndata.index

    # sum identical columns
    tdata = tdata.groupby(tdata.columns, axis=1).sum()
    # divide by 16S copy number
    cn = copy_num.loc[tdata.columns].values.T[0]
    tdata = tdata.div(cn)
    # multiply by KEGG profile table
    kegg_data = kegg_profiles.loc[tdata.columns]
    kegg_data.fillna(0, inplace=True)
    results = jdata.dot(kegg_data)
    # save results
    results.to_csv(args.out, index=True, sep='\t')

    # get full results
    fsamps, ftaxa, fabund, fkegg = [], [], [], []
    ttdata = tdata.copy()
    ttdata['sample'] = ttdata.index
    tmelt = pd.melt(ttdata, id_vars='sample', var_name='taxa',\
                    value_name='abundance')
    for i in tmelt.index:
        samp = tmelt.loc[i, 'sample']
        taxa = tmelt.loc[i, 'taxa']
        abund = tmelt.loc[i, 'abundance']
        kdata = kegg_data.loc[taxa]
        tkdat = list(abund*kdata.values)
        cnum = kdata.shape[0]
        fsamps += [samp]*cnum
        ftaxa += [taxa]*cnum
        fabund += tkdat
        fkegg += list(kdata.index)
    #store in dataframe and save
    full_results = pd.DataFrame({'samples':fsamps,\
                                 'features':ftaxa,\
                                 'KEGG_ID':fkegg,\
                                 'abundance':fabund})
    full_results.to_csv(args.full_out, index=True, sep='\t')
    return 

def main():
    parser = argparse.ArgumentParser(description="Create functional profile\
                                     from 16S data")
    #add arguments for abundance table, representative sequences,
    # additional parameters: % identity
    parser.add_argument("tax_table", help="Taxonomy abundance table.")
    parser.add_argument("rep_fasta", help="Representative sequences, fasta\
                        format.")
    parser.add_argument("--out", help="Outfile for functional profile.",\
                       default="pyphillin_results.tsv")
    parser.add_argument("--blast_in", help="Infile for blast results.")
    parser.add_argument("--blast_out", help="Outfile for blast results.",\
                       default="blast_results.tsv")
    parser.add_argument("--full_out", help="Outfile for full results.",\
                       default="full_results.tsv")
    parser.add_argument("-pid", "--pct_id", help="Percent identity cutoff.",\
                       default=0.97, type=float)
    args = parser.parse_args()

    #align sequences
    if args.blast_in is None:
        align_sequences(args)
    else:
        args.blast_out = args.blast_in

    #create functional profile
    profile_functions(args)
    return

if __name__ == '__main__':
    main()


