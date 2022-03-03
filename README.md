# NetBCE
NetBCE Enables Accurate Prediction of Linear B-Cell Epitopes with Interpretable Deep Neural Network

Activated B-lymphocytes (B cells) produce antibodies that bind with specific antigens, and are a key component in vertebrate immune responses. Thus, identification of B-cell epitopes (BCEs) plays an essential role in the development of peptide vaccines, immuno-diagnostic reagents and antibody production. Here, we obtained over 1.3 million B cell assays with experimentally identified BCE regions from IEDB database. Through quality control procedures, an experimentally well-characterized dataset was compiled, containing more than 126,000 experimentally epitope-containing regions from 3567 protein clusters. Numerous widely used sequence and structural features was encoded and benchmark tested by six conventional machine-learning algorithms. The result shown that different types of features displayed various accuracies for B cell epitope prediction and sequence features had superior performance compared to structural features. To learn a more efficient and interpretive representation of the epitope sequence hierarchically, a ten-layer deep learning framework, named NetBCE, was implemented to predict B cell epitopes. NetBCE achieved high accuracy and robust performance with the average AUC values of 0.8455 by 5-fold cross validation through automatically learning informative classification features. In comparison, NetBCE outperformed conventional machine learning methods and other existing tools based on the curated independent dataset, and achieved a≥8.84% improvement of AUC value for the B cell epitope prediction compared to other tools. To elucidate the capability of hierarchical representation by NetBCE, we visualized the epitopes and non-epitopes using UMAP based on the feature representation at varied network layers. We found the feature representation came to be more discriminative along the network layer hierarchy.

# Installation
Download NetBCE by
```
git clone https://github.com/BioDataStudy/NetBCE
```
Installation has been tested in Linux with Python 3.7.
Since the package is written in python 3x, python3x with the pip tool must be installed.
NetBCE uses the following dependencies: numpy, scipy, pandas, h5py, keras version=2.3.1, tensorflow=1.15. You can install these packages by the following commands:
```
conda create -n NetBCE python=3.7
pip install pandas
pip install numpy
pip install scipy
pip install h5py
pip install -v keras==2.3.1
pip install -v tensorflow==1.15
```

# Data collection and processing
To establish a reliable model, an experimentally well-characterized dataset was compiled as follows (Figure 1). From Immune Epitope Database (available at https://www.iedb.org/), the most comprehensive database holding the largest number of experimentally identified epitopes and non-epitopes, we obtained over 1.3 million B cell assays. Each entry contained an antigen protein sequence with a marked region, in the following called “epitope-containing region”, that was an experimentally verified epitope or non-epitope. Protein sequences were retrieved from National Center for Biotechnology Information (NCBI) and UniProt database based on the antigen protein IDs provided in the epitope entry. We preprocessed and filtered the dataset in several criteria (Figure 1). First, identical protein sequences were integrated and all related information about epitope-containing region were aggregated. Second, sequence redundancy for those proteins of non-identical but highly similar was cleared. Using CD-HIT program, all proteins were clustered with a threshold of 90% sequence similarity. For each cluster, only the protein containing the largest number of epitope-containing region was retained. To ensure high confidence of the dataset, each epitope assay was carefully explored and regarded as positive hit only when it has been tested as positive in two or more different B-cell assays, whereas region tested in at least two assays and not tested positive in any assay were regarded as non-epitopes. In addition, all candidate epitope were excluded that had lower than 5 or greater than 25 amino acid residues from the dataset since the number of these epitopes only account for about 1%, whereas including them may result in outliers during model development. Overall, the non-redundant dataset for training and testing contained 28,714 positive and 98,065 negative epitope-containing region from 3567 protein sequence clusters, respectively.

![image](https://github.com/BioDataStudy/NetBCE/blob/main/data/github_1.jpg)

# Performance
NetBCE outperformed conventional machine learning methods by increasing the area under curve (AUC) value by a range of 8.77-21.58% in same training dataset. Moreover, we compared NetBCE with other existing tools based on the curated independent dataset. NetBCE had high performance with the AUC values of 0.8438 on the independent dataset, and achieved a ≥ 8.84% improvement of AUC value for the B cell epitope prediction compared to other tools.

![image](https://github.com/BioDataStudy/5mC-Finder/blob/2d195b681b89259e738c0ba3bcce5dee25c2c08e/prediction/performance.png)

# Interpretability
To decipher the capability of the hierarchical representation and learning, we visualized the m5U and non-m5U sites using UMAP (Uniform Manifold Approximation and Projection) method based on the feature representations uncovered at different network layers. the learned features turn into more and more discriminative along the layer hierarchy, with m5U and non-m5U sites mixed at the input layer without clear boundary, culminating with a clear separation in the output layer.

![image](https://github.com/BioDataStudy/5mC-Finder/blob/99a4038ca69585ac5e23dae074a9f296d66850d7/umap/Uamp_testing.png)

# Motifs
The kernels in the first convolutional layer distinguish important weight matrices over the input sequences to recognize significant patterns. Therefore, we decoded all filters in the convolutional layer of 5mU-Finder and converted them into motifs. As a result, a total of 135 informative motifs were characterized, and several novel RNA mo-tifs with the consensus sequence were discovered, such as CxGGGAxC and GGGxUCG. The clustering analysis of the informative motifs re-vealed their enrichments and higher distributions in positive instances .

![image](https://github.com/BioDataStudy/5mC-Finder/blob/main/motif/Figure%202.jpg)

# Usage
Please cd to the 5mC-Finder/prediction/ folder which contains predict.py.
Example: 
```
python predict.py -f ../testdata/test.fasta -o ../testdata/test_result.txt
```
For details of other parameters, run:
```
python predict.py --help
```

# Citation
Please cite the following paper for using: Deciphering the location and consensus patterns of RNA 5-methyluridine sites by deep learning. Submission 2021.
