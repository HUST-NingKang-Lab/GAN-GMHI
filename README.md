# GAN-GMHI
GAN-GMHI framework consists of three stages, constructing a dataset containing phenotype and batch information for all samples, and then GAN guiding the batch effect correction of raw data, the corrected datasets are output as the training data set for GMHI prediction (see Figure). It is worth noting that the datasets to be batch-corrected by GAN must be classified based on the phenotype first, and the sub-data sets of each phenotype are regrouped according to the batch. To ensure that the unwanted technical variations among different datasets are eliminated, but the biological differences between different phenotypes are not diminished.
<img src="figures/1-1.jpg">
## Data
We have performed a comprehensive analysis on 2,636 healthy and 1,711 non-healthy (including 12 disease phenotypes) individuals’ stool metagenomes from 34 published studies (Gupta, et al.). The taxonomy abundance table is available at the "/data" directory.
## Scripts
The source codes used to reproduce all the results of this study is available at the "/scripts" directory.
