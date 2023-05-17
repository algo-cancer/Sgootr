import argparse, yaml, os, networkx as nx, numpy as np, seaborn as sns
import pandas as pd, matplotlib.pyplot as plt
from os.path import exists

def parse_po(partial_order):
    return

def parse_configs(config, patient):

    with open(config, 'r') as f:
        configs = yaml.safe_load(f)
    directory = os.path.abspath(os.path.join(configs['OUTDIR'], patient))
    palettes = configs['PATIENTS'][patient]['palette']
    partial_order = configs['PATIENTS'][patient]['partial_order']   

    return directory, palettes, partial_order


def find_max_iter(directory):

    cur_max = 0
    while exists(os.path.join(directory,'t{}'.format(cur_max+1),'RF.txt')):
        cur_max += 1
    return cur_max


def find_star(max_iter, directory):

    with open(os.path.join(directory,'t{}'.format(max_iter),'RF.txt'), 'r') as f:
        RF_dists = [int(line.strip().split('\t')[-1]) for line in f.readlines()]

    df = pd.DataFrame({'iter':range(1,max_iter+1), 'RF($t_i$, $t_{i-1}$)': RF_dists})
    plt.figure(figsize=(6,3))
    ax = sns.lineplot(x=df['iter'], \
                      y=df['RF($t_i$, $t_{i-1}$)']).set(xlim=(0,max_iter+1), \
                                                        xticks=range(0,max_iter+1,5), \
                                                        xlabel='iteration')
    b = RF_dists[::-1]
    star = len(b) - np.argmin(b)
    plt.axvline(star, 0,225, color='r')
  
    plt.savefig(os.path.join(directory, 'RF.png'), bbox_inches='tight')
 
    return star
    


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Get Sgootr output.')
    parser.add_argument('-c', required=True, type=str, \
                        help='path to the configuration file')
    parser.add_argument('-p', required=True, type=str, \
                        help='patient name as appeared in the configuration file')

    args = parser.parse_args()
    directory, palettes, partial_order = parse_configs(args.c, args.p)

    #########################################################
    # First, we find iteration * given results from available
    # iterations
    #########################################################

    # find maximum iteration with available information, accounting
    # for potiential early stopping due to lack of shared sites between
    # some pair of cells due to iterative pruning
    max_iter = find_max_iter(directory)  
    star = find_star(max_iter, directory)
    

    #########################################################
    # Write to results.pdf
    #########################################################
    html_path = os.path.join(directory, 'results.html')
    with open(html_path, 'w') as f:
        f.write(("<h1 id='Header'>Results from <code>Sgootr</code></h1>\n"
                 "<h2><code>Sgootr<code> has determined t{} as the output "
                 "iteration t*</h2>\n"
                 "<p><img src='{}' alt='Tree distances across iterations' "
                 "style='max-height: 250px'/></p>\n"
                 "<h2>Here is the tree <code>Sgootr</code> obtained with leaves colored according "
                 "to: </h2>\n").format(star, 'RF.png'))
        for label in palettes.keys():
            f.write(("<h2>{}</h2>\n"
                     "<p><img src='{}' alt='Tree with leaves colored according to {}' "
                     "style='max-width: 800px'/>"
                     "</p>\n").format(label, \
                                      os.path.abspath(os.path.join(directory, \
                                                                   't{}'.format(star), \
                                                                   '{}.png'.format(label))), \
                                      label))
        if partial_order != None:
            f.write('')# TODO
            
    pdf_path = os.path.join(directory, 'results.pdf')
    os.system('wkhtmltopdf {} {} 2&> /dev/null'.format(html_path, pdf_path))
    os.system('rm {}'.format(html_path))
    print('t* = t{}, see the following summary file for detail:\n{}'.format(star, pdf_path))
