This is the repository for `Sgootr`, a tool that jointly infers a tumor lineage tree and selects lineage-informative features using single-cell methylation sequencing data.

![Schema Figure for Sgootr](/assets/sysarch.png)

# Table of Contents

  1. [Getting Started](#start) 
     * [Setting Up](#setup): how to download the required tools and programs
     * [Example](#example): a step-by-step guide to reproduce our main result on metastatic colorectal cancer patient CRC01 [^1]
  2. [How to Use `Sgootr`](#manual)
     * [Configurations](#config): descriptions of program parameters in [`config.yaml`](https://github.com/algo-cancer/Sgootr/blob/main/config.yaml)
     * [Files](#files) 
       * [Input](#input): content and format of input files to `Sgootr`
       * [Ouput](#output): content and format of output files to `Sgootr`
       * [Intermediary](#intermediary): content and format of select intermediary file to `Sgootr`
     * [The Optional Biclustering Module](#biclustering): how to perform the optional biclustering procedure
  3. [Support and Contact](#support)

<a name="start"></a>
# Getting Started

To help you get started with using `Sgootr`, we will first describe how to set up the required tools and programs, then lead you through an example that reproduces our main result on metastatic colorectal cancer patient CRC01 [^1].

<a name="setup"></a>
## Setting Up

Follow instructions to [install `conda`](https://conda.io/projects/conda/en/latest/user-guide/install/).

Follow instructions to [install `Snakemake`](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html).

Follow instructions to [install `R`](https://www.r-project.org/).

Additionally, if you would like to use the biclustering module of `Sgootr`, follow instructions to [install `Gurobi`](https://support.gurobi.com/hc/en-us/articles/360044290292-How-do-I-install-Gurobi-for-Python-), retrieve a `Gurobi` license (which is [free for academics](https://www.gurobi.com/academia/academic-program-and-licenses/)), and [set up the license](https://www.gurobi.com/documentation/9.5/quickstart_mac/retrieving_and_setting_up_.html) on your computing environment.

Then:

```console
 $ git clone https://github.com/algo-cancer/Sgootr.git
 $ cd Sgootr
```

You should be all set!


<a name="example"></a>
## Example 

Once you have finished [setting up](#setup), we can apply `Sgootr` to the main dataset of interest from our paper, multiregionally-sampled metastatic colorectal patient CRC01 made available by Bian *et al.* [^1]. First, download the preprocessed data `data.tar.gz` into `Sgootr/` from [here](https://umd.box.com/v/sgootr-crc01), then we run `Sgootr` with the [default configurations](https://github.com/algo-cancer/Sgootr/blob/main/config.yaml):

```console
 $ tar -xf data.tar.gz
 $ snakemake --cores <number of cores> --use-conda
```

Once the run has finished, to reproduce panels in Figure 2 in our paper, in the same working directory `Sgootr/`:

```console
 $ conda env create -f envs/sgootr.yml
 $ conda activate sgootr
(sgootr) $ python scripts/get_results.py -c config.yaml -p CRC01
(sgootr) $ conda deactivate
```

This should generate the report `CRC01/results.pdf` that summarizes the results.

<a name="manual"></a>
# How to Use `Sgootr`

If you have successfully run [the example](#example) above, let's now see how you can run `Sgootr` with your own data. If you had trouble running [the example](#example) above, please [contact us](#support). 

We will describe the configurations and input files used by `Sgootr`, and the final and select intermediary files of `Sgootr`.

<a name="config"></a>
## Configurations

`Snakemake` executes the `Sgootr` rules according to the program parameters in [`config.yaml`](https://github.com/algo-cancer/Sgootr/blob/main/config.yaml). We describe them below.

 **Parameter**                  | **Default Value** | **Symbol Used in Paper** | **Description** 
--------------------------------|-------------------|--------------------------|----------------
 `OUTDIR`                       | `./`              | | absolute path to directory in which `Sgootr` will create a directory containing the program output, ending with '/'; default is the current working directory
 `DEFAULT_CNA`                  | `2`               | | default copy number (e.g. 2 in a diploid human genome), a natural number 
 `E`                            | `.01`             | $\epsilon$ | per-base sequencing error rate, a realistically small decimal \[0,1)
 `SITE_COVERAGE_THRESHOLD`      | `.66`             | | used in post-correction filtering, removing CpG sites covering fewer than this fraction of cells, a decimal \[0,1] 
 `CELL_COVERAGE_THRESHOLD`      | `.66`             | | used in post-correction filtering, removing cells covering fewer than this fraction of remaining sites, a decimal \[0,1]
 `P00`                          | `.33`             | P($0^c$) [^\*] | prior probability for the homozygous unmethylated status, a decimal (0,1)
 `P10`                          | `.33`             | P(mixed) [^\*]| prior probability for the heterozygous methylated status, a decimal (0,1)
 `P11`                          | `.33`             | P($1^c$) [^\*]| prior probability for the homozygous methylated status, a decimal (0,1)
 `P`                            | `.5`              | $p$ [^\*]| probability of sampling the methylated allele from a pair of heterozygously methylated alleles, a decimal (0,1) 
 `D0010`                        | `.5`              | dist(00,10), dist(10,00) | distance between homozygous unmethylated status and heterozygous methylated status, a positive number
 `D0011`                        | `1`               | dist(00,11), dist(11,00) | distance between homozygous unmethylated status and homozygous methylated status, a positive number
 `D1011`                        | `.5`              | dist(10,11), dist(11,10) | distance between heterozygous methylated status and homozygous methylated status, a positive number
 `ALGO`                         | `N`               | | distance-based algorithm to construct the tree at each iteration, **(F)** astME [^2] or **(N)** eighbor-joining [^3]
 `MAX_ITER`                     | `40`              | | the maximum number of iterations of pruning `Sgootr` will perform
 `KAPPA`                        | `.1`              | $\kappa$ [^a]| the fraction of sites `Sgootr` will prune each iteration, a decimal (0,1) 
 `PARTITION_VALIDITY_THRESHOLD` | `.5`              | $\omega$ | fraction of cells with heuristically called methylation status within a subtree for a partition to be valid, a decimal (0,1)
 `MINIMUM_SUBTREE_SIZE`         | `.05`             | $\delta$ | fraction of total cells a subtree should have for a partition to be valid, a decimal (0,1)


<a name="files"></a>
## Files

Here we will describe the content and format for various input, output, and intermediary files or information.

<a name="input"></a>
### Input

The input files (or information) for each patient will be specified under `PATIENTS:` in `config.yaml`, one entry per patient. Users may specify multiple patients at once by enumerating multiple entries under `PATIENTS:`, and when the program is run, `Sgootr` will be applied to all specified patients simultaneously [^6]. For each patient entry, one needs to specify the following:

**Input**         | **Required?** | **Description**
------------------|---------------|----------------
`methylated_rc`   | ✔️ | `.npz` file with fields `m` (cell-by-site `numpy` matrix where each entry contains the number of reads supporting a the _methylated_ status of a particular CpG site in a cell [^b]), `rows` (`numpy` array of cell names), and `cols` (`numpy` array of CpG site names)
`unmethylated_rc` | ✔️ | `.npz` file with fields `m` (cell-by-site `numpy` matrix where each entry contains the number of reads supporting a the _unmethylated_ status of a particular CpG site in a cell [^b]), `rows` (`numpy` array of cell names), and `cols` (`numpy` array of CpG site names)
`cna`             |   | `.npz` file with fields `m` (cell-by-site `numpy` matrix where each entry contains the called copy number of a particular CpG site in a cell [^c]), `rows` (`numpy` array of cell names), and `cols` (`numpy` array of CpG site names); leave blank if copy number calls are not available
`whitelist`       |   | `.txt` file containing cells users wish to keep regardless of filtering rules, one cell name per line; leave blank if not applicable
`root`            | ✔️ | name of the cell users wish to use as root for the constructed single cell methylation lineage tree
`labels`          | ✔️ | `.csv` file with header `,cell,<label_type_1>,<label_type_2>,...`, where each subsequent row contains the index (starting at `0`), cell name, and value for `label_type_1`, `label_type_2`, etc. for the cell; cell names should appeear in the same order as they do in `methylated_rc`, `unmethylated_rc`, and `cna` (if available) `rows` attribute
`partial_order`   |   | `.csv` file with header corresponding to the distinct locations of origin of single cells given by one of the label types in `labels`, and with values being the binary adjacency matrix representation of the directed acyclic graph of the given partial order; leave blank if one does not wish to infer the migration history with `Sgootr` or if location labels are not available. 
`palette`         | ✔️ | color palette according to which `Sgootr` will visualize the constructed tree; one entry per provided label type column in `labels`, and each entry is a mapping from each distinct value in the label type to a hexadecimal color code for cells of that label value

Please see the sample input files from the [example](#example) above and the provided `config.yaml` file for a concrete example input to `Sgootr`.

To obtain summary results, we run command `python scripts/get_results.py -c <config_file> -p <patient>`. We describe the inputs below:

**Input**      | **Description**
---------------|-----------------------------------------------
`config_file`  | typically, it would be `config.yaml` in the `/Sgootr` directory; however, if alternative configuration files are used, please specify the absolute path
`patient`      | the key for an entry under `PATIENTS:` in `config_file`; result files will be located in the output directory for that patient


<a name="output"></a>
### Output

Here we describe the output files `Sgootr` produces.

**File**              | **Description**
----------------------|----------------
`results.pdf`         | `.pdf` report summarizing the `Sgootr` output
`RF.png`              | visualization of Robinson-Foulds distances [^4] between trees `Sgootr` constructed in adjacent iterations, and the red vertical line marks iteration `t*`
`{label_type}.png`    | visualization of `t*`, where the single cells at the leaves are colored according to the color palette specified for `label_type`
`migration_graph.csv` | (available if a `partial_order` file is specified in `config.yaml`) `.csv` file with header corresponding to the distinct locations of origin of single cells given in the `partial_order` file, and with values being the adjacency matrix representation of the migration graph
`migration_graph.png` | (available if a `partial_order` file is specified in `config.yaml`) visualization of the migration history inferred by `Sgootr`

<a name="intermediary"></a>
### Intermediary

Here we describe select intermediary files in the output directory that may be of interest for the users.

**File**                        | **Description**
--------------------------------|----------------
`input.npz`                     | post-correction, post-filtering input to the iterative procedure of `Sgootr`, with fields `n` (cell-by-site _methylated_ read count `numpy` matrix), `m` (cell-by-site _unmethylated_ read count `numpy` matrix), `cna` (cell-by-site copy number `numpy` matrix), `rows` (list of cells), and `cols` (list of CpG sites)
`status_likelihoods.npz`        | `.npz` file with fields `p00`, `p10`, and `p11`, each being a cell-by-site `numpy` matrix where each entry containing the computed likelihoods of the observed reads for a CpG site in a cell being drawn from underlying homozygous unmethylated (referred to as P($0^c$ \|reads) in our paper), heterozygous (referred to as P(mixed\|reads) in our paper), or homozyous methylated alleles (referred to as P($1^c$\|reads) in our paper), respectively
`t{iter}/tree.nwk`              | the single cell methylation tumor lineage tree `Sgootr` constructs at iteration `iter`, where the cells are named by their index in the `rows` field of `input.npz` 
`t{iter}/{label_type}.png`      | visualization of `t{iter}`, where the single cells at the leaves are colored according to the color palette specified for `label_type`
`t{iter}/RF.txt`                | the lower triangular matrix (excluding the diagonal) of an {`iter`+1}-by-{`iter`+1} matrix, where an entry in row `i` and column `j` denotes the Robinson-Foulds distance [^4] [^5] between the single cell methylation lineage trees `Sgootr` constructed at iteration `i` and iteration `j`
`t{iter}/site_mask.npz`         | `.npz` file with field `mask`, which is a `numpy` array whose dimension is the number of total CpG sites considered in the iterative procedure of `Sgootr`. Entry `mask[j]` will be `np.inf` if CpG site `j` is used in constructing the single cell methylation tumor lineage tree at iteration `iter`; otherwise, entry `mask[j]` will be the partucular iteration site `j` was pruned out
`t{iter}/persistence_scores.npz`| `.npz` file with fields `scores` (`numpy` array of dimension equal to the number of CpG sites in the `cols` field of `input.npz`, where entry `j` corresponds to the persistence score of site `j` measured on the tree `t{iter-1}`) and `nodes` (`numpy` array of dimension equal to the number of CpG sites in the `cols` field of `input.npz`, where entry `j` corresponds to the name of internal node whose induced bipartitions result in the measured persistence score for site `j`)

<a name="biclustering"></a>
## The Optional Biclustering Module

Aside from the main procedure of `Sgootr`, we introduce in our paper an optional biclustering module to perform coordinated cell and CpG site filtering. See the Methods section of our paper for a description of the formulation, and see Supplemental Section S1 of our paper for further analysis and benchmarking of the formulation. 

To run biclustering, first ensure that [Gurobi has been properly setup](#setup) and that [the `conda` environment for `Sgootr` has been activated](#example), then from the `Sgootr/` directory:

```console
(sgootr) $ python scripts/bicluster.py -i <input> -o <output> -a <alpha> -b <beta> -c <threads> -t <runtime>
```

The required inputs are described below:

**Input**         | **Symbol Used in Paper** | **Description**
------------------|--------------------------|----------------
`input`           | $M$                      | `.npz` file with fields `m` (binary cell-by-site `numpy` matrix, where `m[i][j] = 1` if site $j$ is covered (by at least one read) in cell $i$, and `m[i][j] = 0` otherwise)
`output`          |                          | desired path to `.npz` file the biclustering command will generate with fields `rows` (`numpy` array of indices of selected cells) and `cols` (`numpy` array of indices of selected sites), with which users can obtain the submatrices
`alpha`           | $\alpha$ [^&]| biclustering parameter denoting the fraction of cells the users would like to retain in the resulting submatrix; a decimal (0,1]
`beta`            | $\beta$  [^&]| biclustering parameter denoting the fraction of sites the users would like to retain in the resulting submatrix; a decimal (0,1]
`threads`         | | number of threads the user would like to use for biclustering; a natural number
`runtime`         | | number of seconds expended by `Gurobi` before optimization returns witt a `TIME_LIMIT` status

<a name="support"></a>
# Support and Contact

Thank you for using `Sgootr`. The software and corresponding documentations are maintained by the research group of [Dr. S. Cenk Sahinalp](https://algo-cancer.github.io/). If you have encountered issues with the tool, please report on the [issue forum](https://github.com/algo-cancer/Sgootr/issues) or contact Yuelin Liu [[email]](yuelin@cs.umd.edu). Please always refer to the [official GitHub repository of `Sgootr`](https://github.com/algo-cancer/Sgootr) for the most updated version of the software and relevant documentation.

<!-- References -->

[^1]: Bian, S., Hou, Y., Zhou, X., Li, X., Yong, J., Wang, Y., Wang, W., Yan, J., Hu, B., Guo, H., Wang, J.,
Gao, S., Mao, Y., Dong, J., Zhu, P., Xiu, D., Yan, L., Wen, L., Qiao, J., Tang, F., Fu, W.: Single-cell multiomics sequencing and analyses of human colorectal cancer. Science **362**(6418), 1060-1063 (Nov 2018). [https://doi.org/10.1126/science.aao3791](https://doi.org/10.1126/science.aao3791)

[^2]: Lefort, V., Desper, R., Gascuel, O.: FastME 2.0: A Comprehensive, Accurate, and Fast Distance-Based Phylogeny Inference Program. Molecular Biology and Evolution **32**(10), 2798-2800 (Oct 2015). [https://doi.org/10.1093/molbev/msv150](https://doi.org/10.1093/molbev/msv150)

[^3]: Saitou, N., Nei, M.: The neighbor-joining method: a new method for reconstructing phylogenetic trees. Molecular Biology and Evolution **4**(4), 406-425 (Jul 1987). [https://doi.org/10.1093/oxfordjournals.molbev.a040454](https://doi.org/10.1093/oxfordjournals.molbev.a040454)

[^\*]: See Supplemental Section S8 of our paper for tips on how to choose these parameters.

[^a]: See Supplemental Section S9 for our paper for tips on how to choose these parameters.

[^b]: We distinguish a value of `0` and a value of `np.nan` in the cell-by-site read count matrix. If the `methylated_rc['m']` matrix has value of `0` in cell `i` and site `j`, it means site `j` in cell `i` is covered by some reads but all of them are unmethylated; however, if the `methylated_rc['m']` matrix has the value `np.nan` in cell `i` and site `j`, it means site `j` is not covered at all in cell `i`.

[^c]: Similar to the read count matrices, we distinguish a value of `0` and a value of `np.nan` in the cell-by-site copy number matrix, where `0` will denote a complete loss, and `np.nan` denotes the absence of a copy number call (and the copy number for that site in the the cell will be assumed to be the default).

[^4]: Robinson, D.F., Foulds, L.R.: Comparison of phylogenetic trees. Mathematical biosciences **53**(1-2), 131-147 (February 1981). [https://doi.org/10.1016/0025-5564(81)90043-2](https://doi.org/10.1016/0025-5564(81)90043-2)

[^5]: Sul, S.J., Williams, T.L.: An experimental analysis of robinson-foulds distance matrix algorithm. European Symposium on Algorithms, 793-804 (September 2008).

[^6]: One outstanding issue of `Sgootr` is that, if tree construction fails before `MAX_ITER` is reached, the `snakemake` command fails ungracefully. While users may still run the  `python scripts/get_results.py -c <config_file> -p <patient>` command to get results for that patient, if multiple patients are specified in the `config.yaml` file and tree construction fails for one of the patients, the command fails for all patients, potentially terminating the pruning procedures for those patients prematurely. In this case, one should comment out the specifications for the patient on which tree construction fails in `config.yaml`, then run the pipeline with `snakemake --cores <number of cores> --use-conda` again. This will continue the procedure for the remaining patients. Alternatively, one may run `Sgootr` for one patient at a time.

[^&]: See Supplemental Section S1 of our paper for tips on how to choose these parameters.
