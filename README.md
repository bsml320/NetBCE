# NetBCE
NetBCE Enables Accurate Prediction of Linear B-Cell Epitopes with Interpretable Deep Neural Network

Activated B-lymphocytes (B cells) produce antibodies that bind with specific antigens, and are a key component in vertebrate immune responses. Thus, identification of B-cell epitopes (BCEs) plays an essential role in the development of peptide vaccines, immuno-diagnostic reagents and antibody production. Here, we obtained over 1.3 million B cell assays with experimentally identified BCE regions from IEDB database. Through quality control procedures, an experimentally well-characterized dataset was compiled, containing more than 126,000 experimentally epitope-containing regions from 3567 protein clusters. Numerous widely used sequence and structural features was encoded and benchmark tested by six conventional machine-learning algorithms. The result shown that different types of features displayed various accuracies for B cell epitope prediction and sequence features had superior performance compared to structural features. To learn a more efficient and interpretive representation of the epitope sequence hierarchically, a ten-layer deep learning framework, named NetBCE, was implemented to predict B cell epitopes. NetBCE achieved high accuracy and robust performance with the average AUC values of 0.8455 by 5-fold cross validation through automatically learning informative classification features. In comparison, NetBCE outperformed conventional machine learning methods and other existing tools based on the curated independent dataset, and achieved a≥8.84% improvement of AUC value for the B cell epitope prediction compared to other tools. To elucidate the capability of hierarchical representation by NetBCE, we visualized the epitopes and non-epitopes using Uniform Manifold Approximation and Projection (UMAP) based on the feature representation at varied network layers. We found the feature representation came to be more discriminative along the network layer hierarchy.

# Installation
Download NetBCE by
```
git clone https://github.com/BioDataStudy/NetBCE
```
Installation has been tested in Linux with Python 3.7.
Since the package is written in python 3x, python3x with the pip tool must be installed.
5mC-Finder uses the following dependencies: numpy, scipy, pandas, h5py, keras version=2.3.1, tensorflow=1.15. You can install these packages by the following commands:
```
conda create -n NetBCE python=3.7
pip install pandas
pip install numpy
pip install scipy
pip install h5py
pip install -v keras==2.3.1
pip install -v tensorflow==1.15
```

# Performance
NetBCE outperforms 84 conventional machine-learning predictors, with the area under curve (AUC) value greater than 0.97 in 10-fold cross-validations and independent test.

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
