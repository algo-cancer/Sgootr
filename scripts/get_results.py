import argparse, yaml, os, skbio, networkx as nx, numpy as np, seaborn as sns
import pandas as pd, matplotlib.pyplot as plt
from os.path import exists


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
    plt.close()
 
    return star
    

def sort_labels(directory, partial_order, palettes):

    df = pd.read_csv(partial_order)
    labels, adj = list(df.columns), df.values
    palette = [palettes[k] for k in palettes if set(palettes[k].keys()) == set(labels)][0]
    
    G = nx.from_numpy_matrix(adj, create_using=nx.DiGraph)
    G = nx.relabel_nodes(G, {i:labels[i] for i in range(len(labels))})
    tsorted = list(nx.topological_sort(G))

    df = pd.read_csv(os.path.join(directory, 'labels.csv'))
    cell_labels = [df[c].values for c in df.columns if set(df[c].values) == set(tsorted)][0]

    return tsorted, palette, cell_labels


def bottom_up(status_dict, node, leaf_labels):

    if node.is_tip(): # if node is a leaf, corresponding to a single cell
        status_dict[node.name] = {leaf_labels[int(node.name)]}, 0
        return status_dict[node.name]
    
    u,v = node.children
    
    Ru, Cu = bottom_up(status_dict, u, leaf_labels)
    Rv, Cv = bottom_up(status_dict, v, leaf_labels)
    
    if len(Ru.intersection(Rv)) == 0:
        status_dict[node.name] = Ru.union(Rv), Cu+Cv+1
        return status_dict[node.name]
    else:
        status_dict[node.name] = Ru.intersection(Rv), Cu+Cv
        return status_dict[node.name]


def fitch_hartigan(tree, leaf_labels):

    status_dict = {}
    bottom_up(status_dict, tree.root(), leaf_labels)

    return status_dict 


def top_down(node, parent_status, status_dict, order, migration_graph):

    s, _ = status_dict[node.name] 
    ss = np.array(list(s))
    
    if node.is_root():
        status = ss[np.argmin([order.index(s_) for s_ in ss])]
    else:
        if parent_status in ss:
            status = parent_status
        else:
            status = ss[np.argmin([order.index(s_) for s_ in ss])]
            migration_graph[order.index(parent_status)][order.index(status)] += 1
    
    if not node.is_tip():
        top_down(node.children[0], status, status_dict, order, migration_graph)
        top_down(node.children[1], status, status_dict, order, migration_graph)


def infer_migration(tree, status_dict, order):

    migration_graph = np.zeros((len(order), len(order)), dtype=int)
    top_down(tree.root(), None, status_dict, order, migration_graph)
    
    return migration_graph


def plot_migration(directory, order, migration_graph, palette):

    # save migration graph in a file
    m_mat = os.path.join(directory, 'migration_graph.csv')
    np.savetxt(m_mat, migration_graph, delimiter=',', header=','.join(order), fmt='%d')

    # create visualization for the migration graph
    G = nx.from_numpy_matrix(migration_graph, parallel_edges=True, create_using=nx.MultiDiGraph())
    G = nx.relabel_nodes(G, {i:order[i] for i in range(len(order))})
    pos = nx.kamada_kawai_layout(G)

    plt.figure()

    nx.draw_networkx_nodes(G, \
                           pos, \
                           node_color = [palette[loc] for loc in order], \
                           node_size = 2500, \
                           alpha = 1)
    nx.draw_networkx_labels(G, \
                            pos, \
                            font_size=12)
    ax = plt.gca()
    ax.set_xlim([1.3*x for x in ax.get_xlim()])
    ax.set_ylim([1.3*y for y in ax.get_ylim()])
    ax.axis('off')

    for e in G.edges:
        ax.annotate("", \
                    xy=pos[e[0]], \
                    xycoords='data', \
                    xytext=pos[e[1]], \
                    textcoords='data', \
                    arrowprops=dict(arrowstyle='<-', \
                                    connectionstyle='arc3, rad={}'.format(str(.3*e[2]+0.1)), \
                                    shrinkA=25, shrinkB=25, patchA=None, patchB=None))
    m_png = os.path.join(directory, 'migration_graph.png')
    plt.savefig(m_png, bbox_inches='tight')
    
    plt.close()

    return m_mat, m_png


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Get Sgootr output.')
    parser.add_argument('-c', required=True, type=str, \
                        help='path to the configuration file')
    parser.add_argument('-p', required=True, type=str, \
                        help='patient name as appeared in the configuration file')

    args = parser.parse_args()
    directory, palettes, partial_order = parse_configs(args.c, args.p)

    ##########################################################
    ## First, we find iteration * given results from available
    ## iterations
    ##########################################################

    # find maximum iteration with available information, accounting
    # for potiential early stopping due to lack of shared sites between
    # some pair of cells due to iterative pruning
    max_iter = find_max_iter(directory)  
    star = find_star(max_iter, directory)
   

    ##########################################################
    ## Next, if a partial order is given in the config.yaml
    ## file, run the partial-order-informed Fitch-Hartigan
    ## algorithm to infer a migration history
    ##########################################################

    if partial_order != None:
        assert exists(partial_order), "Partial order file does not exist"
        t = skbio.TreeNode.read(os.path.join(directory, \
                                             't{}'.format(star), \
                                             'tree.nwk'))
        order, palette, cell_labels = sort_labels(directory, partial_order, palettes)
        status_dict = fitch_hartigan(t, cell_labels) 
        migration_graph = infer_migration(t, status_dict, order)

        m_mat, m_png = plot_migration(directory, order, migration_graph, palette)


    ##########################################################
    ## Write to results.pdf
    ##########################################################
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
            f.write(("<h2>Visualization for the migration history inferred " \
                     " by <code>Sgootr</code> accounting for partial order</h2>\n" \
                     "<h3>See {} for the adjacency matrix representation.</h3>\n" \
                     "<p><img src='{}' alt='Migration history inferred' " \
                     "style='max-height: 250px'/></p>\n").format(m_mat, m_png))

            
    pdf_path = os.path.join(directory, 'results.pdf')
    os.system("wkhtmltopdf -q {} {}".format(html_path, pdf_path))
    os.system('rm {}'.format(html_path))
    print('t* = t{}, see the following summary file for detail:\n{}'.format(star, pdf_path))
