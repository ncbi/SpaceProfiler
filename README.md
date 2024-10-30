# SpaceProfiler
SpaceProfiler: A computational method for translating information between different (but related) types of data
#### Pre-required installations before running SpaceProfiler
Python libraries pandas, numpy, seaborn, scipy, matplotlib and sklearn are prerequired to be installed before running SpaceProfiler
#### Input Data to SpaceProfiler
Quantity vector of a cancer-related phenotype across patents, corresponding bulk gene expression data, and spatially resolved transcriptomic (SRT) data as cvs files
```
BulkExpression = pd.read_csv('Cancer-Bulk-ExprData.csv', index_col=0)    # Load the bulk gene expression matrix (columns are genes and rows are patients)
PhenotypeVector = pd.read_csv('Cancer-Phenotype.csv', index_col=0)       # Load the vector of phenotype quantity matched with the rows of bulk gene expression matrix (a vector across patients)
DataSRT = pd.read_csv('SRTsample-SRT-ExprData.csv', index_col=0) # Load the SRT gene expression matrix (columns are spots and rows are genes)  
MetaData = pd.read_csv('SRTsample-SRT-spot-locations.csv', index_col=0)  # Load the spot locations corresponding to SRT gene expression matrix
```
#### To Run SpaceProfiler
SpaceProfiler takes a phenotype quantity vector and a bulk gene expression matrix as bulk data matched across patients, and SRT gene expression matrix and spot locations matched ascross spots  
```
 [EigenGene, EigenPatient, Result]= SpaceProfiler(PhenotypeVector, BulkExpression, DataSRT, MetaData)  # Predict the phenotype quantity on spots of tissue slice  
```
#### Output of SpaceProfiler
As result, returns the Eigen-Gene vector, Eigen-Patient vector, Cosine similarity vector (Predicted phenotype quantity on spots) and Scatter plot for the predicted phenotype quantity on the spot locations  
```
EigenGene.to_csv('Eigen-Gene.csv')                  # Save the Eigen-Gene vector as a cvs file
EigenPatient.to_csv('Eigen-Patient.csv')            # Save the Eigen-Patient vector as a cvs file
Result['Cosine-Similarity'].to_csv('Predicted-Phentype-On-Spots.csv')  # Save the predicted phenotype quantity vector as a cvs file
plt.savefig('Prediction-Plot.png', dpi=300,  bbox_inches = 'tight')    # Save the prediction plot for phenotype quantity
```
#### Python [Package](code)
* The SpaceProfiler method is implemented in python and the codes are available as [Python Code](code/SpaceProfiler.py) and [Jupyter Notebook](code/SpaceProfiler.ipynb) modules.

#### Bulk data sets [Data Bulk](Bulk-data)
For each of the four cancer data sets (BRCA, COAD LUAD and LUSC), the phenotype quantity of patients (Tumor-Purity, Hazard, Stemness and Proliferation) and bulk gene expression matrices whose columns represent the genes and rows represent patients
#### SRT data sets [Data SRT](SRT-data)
For each of the SRT samples, the spatially resolved transcriptomics (SRT) gene expression data matrix whose columns represent spots and rows represent genes, and spot locations denoting the spacial location of each spot in the tissue slice

#### Results of the bulk data analysis [Result](result)
The following files proves the results for the analysis on bulk data of cancer patients.
* Lists of the Eigen-Genes for cancers and phenotypes with the predicted values over patients [Eigen-Genes](result/Eigen-Genes.xlsx).
* Lists of the Eigen-Patients for cancers and phenotypes with the predicted values over patients [Eigen-Patients](result/Eigen-Patients.xlsx).
* GO terms enriched for the sorted list of top predictor genes [GO terms for predictors](result/GO-terms.xlsx). The GO terms enriched for top negative and positive predictors in Eigen-Patient.
#### Results of the SRT data analysis [Result](result)
The following files proves the results for the analysis on spacial transcriptomics (SRT) data
* Predicted phenotype quantity on spots [Result](result/Cosine-Similaries.xlsx). Prediction based on similarity between Eigen-Patient and SRT gene expression in each spot.
* Plots for showing the prediction results over spatial location [Result](result)
