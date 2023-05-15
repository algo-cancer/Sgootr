'''
    post_filter.py <corrected_rc_cna.npz> <labels.csv> <[optional]whitelist.txt>
'''

import sys, numpy as np, pandas as pd
from datetime import datetime

if __name__ == "__main__":

    f = open(snakemake.log[0], 'w')
    sys.stderr = sys.stdout = f

    f.write('[{}] Sgootr is performing post-correction ' \
            'filtering of the input matrices\n'.format(datetime.now()))

    obj = np.load(snakemake.input[0], allow_pickle=True)
    _N, _M, _C, cells = obj['n'], obj['m'], obj['cna'], obj['rows']
    site_coverage_threshold = float(snakemake.params.site_coverage_threshold)
    cell_coverage_threshold = float(snakemake.params.cell_coverage_threshold)

    assert (site_coverage_threshold >= 0) and (site_coverage_threshold <= 1), \
           "site_coverage_threshold must be a decimal [0,1]."
    assert (cell_coverage_threshold >= 0) and (cell_coverage_threshold <= 1), \
           "cell_coverage_threshold must be a decimal [0,1]."

    if len(snakemake.input) == 3: # optional whitelist input
        with open(snakemake.input[2], 'r') as fi:
            whitelisted_cells = [x.strip() for x in fi.readlines()]
            whitelist = np.array([(x in whitelisted_cells) for x in cells])

    _df = pd.read_csv(snakemake.input[1], index_col=0)

    assert np.array_equal(np.isnan(_N), np.isnan(_M)), \
           "Methylated and unmethylated read matrices have different coverage."
    assert np.array_equal(_df['cell'].values, cells), \
           "Cells in labels file do not match those in read matrices."    

    '''
        given read error corrected input matrices, further remove (in 
        particular order):
        1. sites covering less than <site_coverage_threshold> of total number of
           cells
        2. cells covering less than <cell_coverage_threshold> of the remaining sites
    '''
    covered = ~np.isnan(_N)

    site_coverage = np.sum(covered, axis=0) 

    selected_sites = (site_coverage >= round(site_coverage_threshold * \
                                             _N.shape[0]))

    cell_coverage = np.sum(covered[:,selected_sites], axis=1)

    if len(snakemake.input) == 3: # if there are whitelisted cells
        selected_cells = (cell_coverage >= round(cell_coverage_threshold * \
                                                 np.sum(selected_sites))) | whitelist
    else:
        selected_cells = (cell_coverage >= round(cell_coverage_threshold * \
                                                 np.sum(selected_sites)))    

    f.write('[{}] Sgootr selected {} cells and {} sites\n'.format(datetime.now(), \
                                                                  np.sum(selected_cells), \
                                                                  np.sum(selected_sites)))
    np.savez(snakemake.output[0], n=_N[selected_cells,:][:,selected_sites], \
                                  m=_M[selected_cells,:][:,selected_sites], \
                                  cna=_C[selected_cells,:][:,selected_sites], \
                                  rows=obj['rows'][selected_cells], \
                                  cols=obj['cols'][selected_sites])

    f.write('[{}] Sgootr is subsetting label file\n'.format(datetime.now()))
    df = _df[_df['cell'].isin(cells[selected_cells])].reset_index(drop=True)
    df.to_csv(snakemake.output[1])

    f.write('[DONE]')
    f.close()
