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
-$ git clone https://github.com/algo-cancer/Sgootr.git
-$ cd Sgootr
```

You should be all set!


<a name="example"></a>
## Example 

Once you have finished [setting up](#setup), we can apply `Sgootr` to the main dataset of interest from our paper, multiregionally-sampled metastatic colorectal patient CRC01 made available by Bian *et al.* [^1]. First, download the preprocessed data `data.tar.gz` into `Sgootr/` from [here](https://umd.box.com/v/sgootr-crc01), then we run `Sgootr` with the [default configurations](https://github.com/algo-cancer/Sgootr/blob/main/config.yaml):

```console
-$ tar -xf data.tar.gz
-$ snakemake --cores <number of cores> --use-conda
```

Once the run has finished, to reproduce panels in Figure 2 in our paper, in the same working directory `Sgootr/`:

```console
-$ conda env create -f envs/sgootr.yml
-$ conda activate sgootr
(sgootr)-$ python -m ipykernel install --user --name sgootr
(sgootr)-$ conda deactivate
```

Open `notebooks/CRC01-Analysis.ipynb`, then follow the interactive notebook.

<a name="manual"></a>
# How to Use `Sgootr`

If you have successfully run [the example](#example) above, let's now see how you can run `Sgootr` with your own data. If you had trouble running [the example](#example) above, please [contact us](#support). 

We will describe the configurations and input files used by `Sgootr`, and the final and select intermediary files of `Sgootr`.

<a name="config"></a>
## Configurations

`Snakemake` executes the `Sgootr` rules according to the program parameters in [`config.yaml`](https://github.com/algo-cancer/Sgootr/blob/main/config.yaml). We describe them below.

 **Parameter**                  | **Default Value** | **Symbol Used in Paper** | **Description** 
--------------------------------|-------------------|--------------------------|----------------
 `OUTDIR`                       | `./`              | | path to directory in which `Sgootr` will create a directory containing the program output, ending with '/' 
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

<a name="input"></a>
### Input

<a name="output"></a>
### Output

<a name="intermediary"></a>
### Intermediary

<a name="support"></a>
# Support and Contact

Thank you for using `Sgootr`. The software and corresponding documentations are maintained by the research group of [Dr. S. Cenk Sahinalp](https://algo-cancer.github.io/). If you have encountered issues with the tool, please report on the [issue forum](https://github.com/algo-cancer/Sgootr/issues) or contact Yuelin Liu [[email]](yuelin@cs.umd.edu).

<!-- References -->

[^1]: Bian, S., Hou, Y., Zhou, X., Li, X., Yong, J., Wang, Y., Wang, W., Yan, J., Hu, B., Guo, H., Wang, J.,
Gao, S., Mao, Y., Dong, J., Zhu, P., Xiu, D., Yan, L., Wen, L., Qiao, J., Tang, F., Fu, W.: Single-cell multiomics sequencing and analyses of human colorectal cancer. Science **362**(6418), 1060-1063 (Nov 2018). [https://doi.org/10.1126/science.aao3791](https://doi.org/10.1126/science.aao3791)

[^2]: Lefort, V., Desper, R., Gascuel, O.: FastME 2.0: A Comprehensive, Accurate, and Fast Distance-Based Phylogeny Inference Program. Molecular Biology and Evolution **32**(10), 2798-2800 (Oct 2015). [https://doi.org/10.1093/molbev/msv150](https://doi.org/10.1093/molbev/msv150)

[^3]: Saitou, N., Nei, M.: The neighbor-joining method: a new method for reconstructing phylogenetic trees. Molecular Biology and Evolution **4**(4), 406-425 (Jul 1987). [https://doi.org/10.1093/oxfordjournals.molbev.a040454](https://doi.org/10.1093/oxfordjournals.molbev.a040454)

[^\*]: See Supplemental Section S8 of our paper for tips on how to choose these parameters.

[^a]: See Supplemental Section S9 for our paper for tips on how to choose these parameters.
