### Dataset list in [paper](https://www.nature.com/articles/s41587-022-01353-8#data-availability)

- Data generated as part of this study are available as [supplementary tables](https://static-content.springer.com/esm/art%3A10.1038%2Fs41587-022-01353-8/MediaObjects/41587_2022_1353_MOESM3_ESM.xlsx) or from [link](http://dig-cancer.csail.mit.edu/)

```
wget -c -r --no-parent --progress --no-check-certificate https://cb.csail.mit.edu/cb/DIG/downloads/ 

(base) [j_wang@n02 downloads]$ tree -L 2
.
|-- [ 541]  annotions
|   |-- [   6]  coding
|   |-- [ 490]  index.html
|   |-- [   6]  noncoding
|   |-- [   8]  splicing
|   `-- [  26]  style.css
|-- [3.6G]  dig_data_files
|   |-- [144M]  element_data.h5
|   |-- [180M]  gene_data.h5
|   |-- [152M]  genome_counts.h5
|   |-- [3.0G]  hg19.fasta
|   |-- [3.5K]  hg19.fasta.fai
|   |-- [ 912]  index.html
|   |-- [144M]  sites_data.h5
|   `-- [  26]  style.css
|-- [ 473]  examples
|   |-- [  13]  PCAWG_ICGC_subset
|   |-- [  12]  PCAWG_all
|   |-- [ 418]  index.html
|   `-- [  26]  style.css
|-- [ 555]  mutation_files
|   |-- [  31]  Dietlein_2019
|   |-- [   5]  PCAWG
|   |-- [ 502]  index.html
|   `-- [  13]  megacohorts
`-- [2.4G]  mutation_maps
    |-- [ 67M]  Adenocarcinoma_tumors_SNV_MNV_INDEL_msi_low.Pretrained.h5
    |-- [ 67M]  Biliary-AdenoCA_SNV_MNV_INDEL.Pretrained.h5
    |-- [ 77M]  Bladder-TCC_SNV_MNV_INDEL.Pretrained.h5
    |-- ....

```

- Browsable mutation maps for 37 cancer types are provided at [link](https://resgen.io/maxsh/Cancer_Mutation_Maps/views)

- PCAWG data are available from  [link](https://dcc.icgc.org/releases/PCAWG/)

```
(base) [j_wang@n02 downloads]$ wget -c -r --no-parent --progress --no-check-certificate https://cb.csail.mit.edu/cb/DIG/downloads/mutation_files/PCAWG/ICGC_only/

```
- Hartwig Medical Foundation data are available from  [link](https://database.hartwigmedicalfoundation.nl/)

- Whole-exome sequencing data compiled by Dietlein et al. are available from  [link](http://www.cancer-genes.org/)

- Targeted sequencing data are available from  [link](https://www.cbioportal.org/)

- The list of genes in the Cancer Gene Census is available at  [link](https://cancer.sanger.ac.uk/cosmic/download)

- Regions and model predictions for each fold of each cancer for which a deep learning model was trained. Regions that were filtered are also listed along with their respective model predictions as [Supplementary Data](https://static-content.springer.com/esm/art%3A10.1038%2Fs41587-022-01353-8/MediaObjects/41587_2022_1353_MOESM4_ESM.zip): <br>

- 36mer Mappability: [link](https://genome.ucsc.edu/cgi-bin/hgFileUi?db=hg19&g=wgEncodeMapability) and [download page](https://hgdownload.soe.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeMapability/)

- Replication Timing from ten cell lines from ENCODE : in the [supplementary tables - page T3](https://static-content.springer.com/esm/art%3A10.1038%2Fs41587-022-01353-8/MediaObjects/41587_2022_1353_MOESM3_ESM.xlsx) 

- 723 chromatin marks in 111 tissues from Roadmap Epigenomics : in the [supplementary tables - page T3](https://static-content.springer.com/esm/art%3A10.1038%2Fs41587-022-01353-8/MediaObjects/41587_2022_1353_MOESM3_ESM.xlsx) 

About more details, please refer to the [download script](https://github.com/jinxin-wang/RepDigDrive/blob/main/scripts/download_data_resources.sh)

### Extended data in paper

1. [Detailed overview of the Dig model](https://www.nature.com/articles/s41587-022-01353-8/figures/5) :

![Dig model](https://media.springernature.com/full/springer-static/esm/art%3A10.1038%2Fs41587-022-01353-8/MediaObjects/41587_2022_1353_Fig5_ESM.jpg)

2. [Epigenetic input features](https://www.nature.com/articles/s41587-022-01353-8/figures/6)

<img src="https://media.springernature.com/full/springer-static/esm/art%3A10.1038%2Fs41587-022-01353-8/MediaObjects/41587_2022_1353_Fig6_ESM.jpg" width="500" >

3. [Cryptic splice SNV enrichment](https://www.nature.com/articles/s41587-022-01353-8/figures/7)

<img src="https://media.springernature.com/full/springer-static/esm/art%3A10.1038%2Fs41587-022-01353-8/MediaObjects/41587_2022_1353_Fig7_ESM.jpg" width="500" >

4. [distribution of activating mutations in gene-tumor pairs](https://www.nature.com/articles/s41587-022-01353-8/figures/8)

<img src="https://media.springernature.com/full/springer-static/esm/art%3A10.1038%2Fs41587-022-01353-8/MediaObjects/41587_2022_1353_Fig8_ESM.jpg" width="500" >


