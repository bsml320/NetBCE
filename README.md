# NetBCE
NetBCE Enables Accurate Prediction of Linear B-Cell Epitopes with Interpretable Deep Neural Network

Activated B-lymphocytes (B cells) produce antibodies that bind with specific antigens, and are a key component in vertebrate immune responses. Thus, identification of B-cell epitopes (BCEs) plays an essential role in the development of peptide vaccines, immuno-diagnostic reagents and antibody production. Here, we obtained over 1.3 million B cell assays with experimentally identified BCE regions from IEDB database. Through quality control procedures, an experimentally well-characterized dataset was compiled, containing more than 126,000 experimentally epitope-containing regions from 3567 protein clusters. Numerous widely used sequence and structural features was encoded and benchmark tested by six conventional machine-learning algorithms. The result shown that different types of features displayed various accuracies for B cell epitope prediction and sequence features had superior performance compared to structural features. To learn a more efficient and interpretive representation of the epitope sequence hierarchically, a ten-layer deep learning framework, named NetBCE, was implemented to predict B cell epitopes. NetBCE achieved high accuracy and robust performance with the average AUC values of 0.8455 by 5-fold cross validation through automatically learning informative classification features. In comparison, NetBCE outperformed conventional machine learning methods and other existing tools based on the curated independent dataset, and achieved a≥8.84% improvement of AUC value for the B cell epitope prediction compared to other tools. To elucidate the capability of hierarchical representation by NetBCE, we visualized the epitopes and non-epitopes using UMAP based on the feature representation at varied network layers. We found the feature representation came to be more discriminative along the network layer hierarchy.

# Installation
Download NetBCE by
```
git clone https://github.com/bsml320/NetBCE
```
Installation has been tested in Linux server (CentOS Linux release 7.8.2003 (Core)) with Python 3.7.
Since the package is written in python 3x, python3x with the pip tool must be installed.
NetBCE uses the following dependencies: numpy, scipy, pandas, h5py, keras version=2.3.1, tensorflow=1.15. You can install these packages by the following commands:
```
conda create -n NetBCE python=3.7
pip install pandas
pip install numpy
pip install scipy
pip install h5py
pip install plotly
pip install dominate
pip install -v keras==2.3.1
pip install -v tensorflow==1.15
```

# Data collection and processing
To establish a reliable model, an experimentally well-characterized dataset was compiled as follows. From Immune Epitope Database (available at https://www.iedb.org/), the most comprehensive database holding the largest number of experimentally identified epitopes and non-epitopes, we obtained over 1.3 million B cell assays. Each entry contained an antigen protein sequence with a marked region, in the following called “epitope-containing region”, that was an experimentally verified epitope or non-epitope. Protein sequences were retrieved from National Center for Biotechnology Information (NCBI) and UniProt database based on the antigen protein IDs provided in the epitope entry. We preprocessed and filtered the dataset in several criteria. First, identical protein sequences were integrated and all related information about epitope-containing region were aggregated. Second, sequence redundancy for those proteins of non-identical but highly similar was cleared. Using CD-HIT program, all proteins were clustered with a threshold of 90% sequence similarity. For each cluster, only the protein containing the largest number of epitope-containing region was retained. To ensure high confidence of the dataset, each epitope assay was carefully explored and regarded as positive hit only when it has been tested as positive in two or more different B-cell assays, whereas region tested in at least two assays and not tested positive in any assay were regarded as non-epitopes. In addition, all candidate epitope were excluded that had lower than 5 or greater than 25 amino acid residues from the dataset since the number of these epitopes only account for about 1%, whereas including them may result in outliers during model development. Overall, the non-redundant dataset for training and testing contained 28,714 positive and 98,065 negative epitope-containing region from 3567 protein sequence clusters, respectively.

![image](https://github.com/BioDataStudy/NetBCE/blob/main/data/github_1.jpg)

# NetBCE construction
With such a sufficient training dataset of B cell assays, the deep neural network can automatically learn informative classification features, making it very appropriate for B cell epitope prediction. Thus, a ten-layer deep learning framework, named NetBCE, was implemented to predict B cell epitopes. The epitope sequences were encoded and taken as input. Then, convolution-pooling module was followed to make feature extraction and representation. LSTM layer was added for retaining features from a long duration facilitate the model to catch the combinations or dependencies among residues at different positions. Lastly, an attention layer was joined to link the LSTM layer and the output layer. 

![image](https://github.com/BioDataStudy/NetBCE/blob/main/models/github_2.jpg)

# Performance
To evaluate the prediction performance of NetBCE, the 5-fold CV was performed on the training dataset. The ROC curves were drawn and the corresponding AUC values were calculated. We found that NetBCE had high performance with the average AUC values of 0.8455 by 5-fold CV, with a range from 0.8379 to 0.8528. Since the number of epitopes and non-epitopes were not balanced in the training dataset, we also performed PR analysis and calculated the corresponding AUC values. The PR curve indicates the trade-off between the amount of false positive predictions compared to the amount of false negative predictions. NetBCE achieved average PR AUC values of 0.6165, suggesting our model had great potential in predicting functional epitopes the high precision. We compared the performance of NetBCE with other 6 ML-based methods regarding AUC value (AB, DT, GB, KNN, LR and RF). We observed that average AUC of NetBCE was 8.77-21.58% higher than those of the other six ML-based methods. In addition, we compared NetBCE with other existing tools based on the curated independent dataset. NetBCE had high performance with the AUC values of 0.8400 on the independent dataset, and achieved a ≥ 22.06% improvement of AUC value for the B cell epitope prediction compared to other tools

![image](https://github.com/BioDataStudy/NetBCE/blob/main/models/github_3.jpg)

# Interpretability
To elucidate the capability of hierarchical representation by NetBCE, we visualized the epitopes and non-epitopes using UMAP (Uniform Manifold Approximation and Projection) method based on the feature representation at varied network layers. We found the feature representation came to be more discriminative along the network layer hierarchy. More specifically, the feature representations for epitopes and non-epitopes sites were mixed at the input layer. As the model continues to train, epitopes and non-epitopes tend to occur in very distinct regions with efficient feature representation. 

![image](https://github.com/BioDataStudy/NetBCE/blob/main/Interpretability/github_4.jpg)

# Usage
Please cd to the NetBCE/prediction/ folder which contains predict.py.
Example: 
```
python NetBCE_prediction.py -f ../testdata/test.fasta -o ../result/test_result
```
For details of other parameters, run:
```
python NetBCE_prediction.py --help
```

# NetBCE analysis report
Based on to the model constructed in this study, we developed a software to provide function for linear B-cell epitope prediction. The software of NetBCE is available at https://github.com/bsml320/NetBCE. NetBCE provided and visualized the prediction results in an interactive html file using the Python, PHP, JavaScript and Bootstrap package with an easily readable and interpretable manner. Users can input the candidate proteins in a FASTA format. In addition, user needs to select one or more peptide lengths so that the software can construct a library of candidate epitope peptides. For an example output page, our software provides a probability score for each candidate peptide, and its value ranges from 0 to 1. All prediction results can be copied, printed and downloaded in 3 formats, including “CVS”, “Excel” and “PDF”. Our software additionally provided two an interactive html plot showing the distribution of lengths and scores for all candidate peptides. 

![image](https://github.com/BioDataStudy/NetBCE/blob/main/prediction/github_5.jpg)

# Citation
Please cite the following paper for using: NetBCE Enables Accurate Prediction of Linear B-Cell Epitopes with Interpretable Deep Neural Network. Submission 2022.
