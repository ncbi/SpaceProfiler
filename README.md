# SpaceProfiler
SpaceProfiler: A computational method for translating translating information between different (but related) types of data
#### Pre-required installations before running SpaceProfiler
Python libraries pandas, numpy, seaborn, scipy, matplotlib and sklearn are prerequired to be installed before running SpaceProfiler
#### Input Data to TranNet
Quantity vector of a cancer-related phenotype across patents, corresponding bulk gene expression data, and spatially resolved transcriptomic data as cvs files
```
BulkExpression = pd.read_csv('Cancer-Bulk-ExprData.csv', index_col=0) # Load the bulk gene expression matrix (columns are genes and rows are patients)
PhenotypeVector = pd.read_csv('Cancer-Phenotype.csv', index_col=0)     # Load the vector of phenotype quantity matched with the rows of bulk gene expression matrix (a vector across patients)
DataSRT = pd.read_csv('SRTsample-SRT-ExprData.csv', index_col=0) # Load the SRT gene expression matrix (columns are spots and rows are genes)  
MetaData = pd.read_csv('SRTsample-SRT-spot-locations.csv', index_col=0) # Load the spot locations corresponding to SRT gene expression matrix
```
#### To Run SpaceProfiler
SpaceProfiler takes two matrices described above and return the transition weight matrix as output
```
SpaceProfiler(PhenotypeVector, BulkExpression, DataSRT, MetaData)  # Predict the phenotype quantity on spots of tissue slice  
```
#### Output of SpaceProfiler
As result, returns the Eigen-Gene vector, Eigen-Patient vector, Cosine similarity vector (Predicted phenotype quantity for spots) and plot for the predicted phenotype quantity on the spot locations  
```
EigenGene.to_csv('Eigen-Gene.csv')                  # Save the Eigen-Gene vector as a cvs file
EigenPatient.to_csv('Eigen-Patient.csv')            # Save the Eigen-Patient vector as a cvs file
Result['Cosine-Similarity'].to_csv('Predicted-Phentype-On-Spots.csv')  # Save the predicted phenotype quantity vector as a cvs file
plt.savefig('Prediction-Plot.png', dpi=300,  bbox_inches = 'tight')    # Save the prediction plot for phenotype quantity
```
#### Python [Package](code)
* The SpaceProfiler method is implemented in python and the codes are available as [Python Code](code/SpaceProfiler.py) and [Jupyter Notebook](code/SpaceProfiler.ipynb) modules.

#### Bulk data sets [Data Bulk](Bulk-data)
For each of the four cancer data sets (BRCA, COAD LUAD, LUSC), the phenotype quantity of patents (Tumor-Purity, Hazard, Stemness and Proliferation) and bulk gene expression matrices whose columns represent the genes and rows represent patients.
#### SRT data sets [Data SRT](SRT-data)
For each of the SRT 11 samples, the spatially resolved transcriptomics (SRT) gene expression data matrix whose columns represent spots and rows represent genes. The location of each spot in tissue slice corresponding to each SRT sample.

#### Results of the analysis on four cancer data sets [Result](result)
The following files proves the results for the analysis on BRCA, COAD, LUAD, LUSC cancer data.
* Lists of the Eigen-Genes for cancers and phenotypes with the predicted values over patients [Eigen-Genes](result/Eigen-Genes.xlsx).
* Lists of the Eigen-Patients for cancers and phenotypes with the predicted values over patients [Eigen-Patients](result/Eigen-Patients.xlsx).
* GO terms enriched for the sorted list of top predictor genes [GO terms for predictors](result/GO-terms.xlsx). The GO terms enriched for top negative and positive predictors in Eigen-Patient.
* Predicted phenotype quantity on spots [Result](result/Predicted Phenotype Phenotype Quantity.xlsx). Prediction based on similarity between Eigen-Patient and SRT gene expression in each spot.
* Plots showing the prediction results [Result](result)
