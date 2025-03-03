# Bioinfo Project 2025

This repository contains materials for the bioinformatics project (2025 L3 INFO, Paris-Saclay). The objective of this project is to identify biomarkers associated with ALS (Amyotrophic Lateral Sclerosis) disease. To achieve this, you will have access to RNA-Seq sequencing data from post-mortem brain cortex biopsies of individuals diagnosed with ALS and those without the disease.

The data for this project originate from the study titled "Postmortem Cortex Samples Identify Distinct Molecular Subtypes of ALS: Retrotransposon Activation, Oxidative Stress, and Activated Glia" by Tam et al. The complete study is available [here](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6866666/).

## Introduction

This README will serve as your guide throughout the project and your analyses. Please follow all the steps listed below. Unless specified otherwise within this document (e.g., mandatory instructions), you are granted complete freedom in choosing how to accomplish each task, particularly in regards to coding techniques.

You will work in pairs from the beginning. While discussions with other pairs are encouraged, direct code sharing (such as copy/paste) is strictly prohibited. Throughout the project, meticulously record both your personal contribution and that of your partner. This information will be required during the final session, which will include a roughly 15-minute oral presentation (details to follow).

In the initial sessions, you will receive lessons or introductions on concepts that may be new to you. Once these lessons begin, please pause all other activities and attend the sessions, even if you are already familiar with the concepts being discussed.

At the conclusion of your project, you are required to submit the following components via email to philippe.rinaudo@universite-paris-saclay.fr:
- All your code,
- A comprehensive report detailing the results of your analyses,
- A text file (either .csv, .txt, or similar) containing a ranked list of your top 100 genes (with the most significant gene listed first).

The subject line of the email should read: "Projet_bioinfo_2025_lastname1_firstname1_lastname2_firstname2". Additionally, please include an archive named lastname1_firstname1_lastname2_firstname2" containing all the required elements mentioned above.

The submission deadline for your project is 9 APRIL 2025.

Note that it is permissible to integrate your code and report within a single notebook, as well as to combine the notebook with standalone code files. However, the report must be consolidated into a single document for submission. The evaluation criteria will prioritize the following, in order:
- Adherence to the provided instructions,
- The clarity of your code,
- The clarity of your report,
- The relevance and efficiency of your code,
- The significance and accuracy of your findings.

## Last Session:

The final session, scheduled for 2nd April 2025, will be dedicated to presentations:
- You will have 15 minutes to present,
- Discuss the current progress of your work,
- Outline your next steps, unless your project is already complete,
- Have your code and a draft of your report ready for discussion,
- There is no need to prepare formal slides or other presentation materials,

Aim to speak for approximately 5 minutes, leaving about 10 minutes for questions and answers.

## Mandatory General Instructions:
- Programming Language: Your project must be coded in Python exclusively,
- Dependencies: Only Python modules are allowed for use, including Jupyter notebooks,
- Archive Content: Extracting your archive should create a folder named "nom_prenom" (your folder),
- Code Execution: Ensure your code runs correctly when your folder is placed within this repository by using relative paths only,
- File Management: Your code must not create any files outside of your designated folder,
- Code Comments: All code should be thoroughly commented,
- Coding Standards: Adhere to PEP 8 guidelines as closely as possible for coding style,
- Testing: Include unit tests within your code,
- Design Paradigm: Your code should be designed in an object-oriented manner as much as possible,

Note: The primary goal is to conduct an analysis. You are not required to develop a bioinformatics pipeline that is robust to data format changes. However, it's crucial to ensure the accuracy of all analyses by any necessary means.

## Final Note:
Within this document, you will find a series of steps intended to guide your analysis. Some steps provide clear instructions on what to do, while others are more open-ended, requiring investigation and thoughtful consideration. Before implementing any analyses you conceive, present them for approval to ensure they align with project goals. Do not rush through the available steps; follow the session pace instead. For more challenging steps, I will guide the discussion and thought process.

# Step 1 - Data Preprocessing

You have access to the raw data for the study in the "Data" folder. This folder contains two types of information:
- The RNA counts for each sample,
- An annotation file providing details on the experiment and, more importantly, the samples.
Begin by downloading the data and initiating the preprocessing.

## Gather RNA Counts

To analyze the samples, you will need to consolidate them into a single Python object. A common approach is to create a dataframe (or any "table-like" structure) where each row represents a "sample" and each column represents a "gene". It's crucial to thoroughly test your dataset to ensure that any future modifications can be easily identified and corrected.

Remember, your code must follow an object-oriented design.

Your code can be specifically tailored to the current state of the data, meaning it should be designed to work with this particular dataset rather than being universally applicable. Below is an example of how you can load the data.

```python
import pandas as pd
import glob
import re

path = "./Data" # the path of the data

pdList = [] # variable to temporary store all dataframes (one for each txt file)
# For all txt file
for fname in glob.glob(path+"/*.txt"):
    df = pd.read_table(fname) # put the file in a dataframe
    sample_name = re.search("GSM\d+", fname).group() # search the name of the sample in the file name
    df.rename(index= df["gene/TE"], inplace=True) # rename the index (=rows) using the column containing the gene name
    df.drop(columns=df.columns[0], axis=1, inplace=True) # drop the first column containing the gene name, no more need
    df.rename(columns={ df.columns[0]: sample_name }, inplace = True) # rename the column (there is only one column at this step) using the sample name
    pdList.append(df) # add the current dataframe in the list
data_matrix = pd.concat(pdList, 1) # concat all dataframe in 1 dataframe
data_matrix = data_matrix.transpose() # transpose the dataframe to get a more standard shape (samples x variables)
```
This code should be integrated into an object-oriented design. Additionally, it assumes that each file is correct and free of errors, which is an aspect that must be verified to ensure thoroughness.

## Gather Sample Annotations

The sample annotations are consolidated into a single "xml" file. Initially, open this file with any text editor to familiarize yourself with its structure. Then, determine which pieces of information are pertinent to your analysis.

Subsequently, construct a dataframe (or any "table-like" structure) where each row represents a sample and each column an annotation. Ensure to rigorously test your dataset (combining gene counts and annotations) to promptly identify and address any discrepancies in subsequent steps.

To parse an XML file, you can employ the "xml.etree.ElementTree" library. It will require manually exploring the file (using a text editor) to understand the structure and identify all relevant sections. The samples are encapsulated within blocks named "Sample," with additional information located in other blocks that you will need to discern.

Below is an example on how to create a dataframe with a single column corresponding to the "Cns_subregion".

```python
data_annotation = pd.DataFrame(columns = ['Sample_id', 'Cns_subregion']) # initialisation of the dataframe
xtree = et.parse('./Data/GSE124439_family.xml') # create a variable containing the xml in a tree shape
xroot = xtree.getroot() # get the root of the tree to start the exploration of the tree/xml
# for each element named "sample" that can be found from the root
for child in xroot.iter("{http://www.ncbi.nlm.nih.gov/geo/info/MINiML}Sample"):
    temp_sample_id = child.attrib['iid'] # the attribut of this node contains the sample id ()
    # for each element named "Characteristics" that can be found from the current sample
    for child2 in child.iter("{http://www.ncbi.nlm.nih.gov/geo/info/MINiML}Characteristics"):
        if(child2.attrib["tag"] == "cns subregion"):
            temp_cns_subregion = child2.text.replace('\n', '')
    temp_df = pd.DataFrame({'Sample_id': [temp_sample_id], 'Cns_subregion': [temp_cns_subregion]})
    data_annotation = pd.concat([data_annotation, temp_df])
```

## Create Initial Preprocessing Functions
At this stage, your code should already include at least one class, accompanied by associated getters and setters. By getters and setters, I refer to methods that facilitate access to an instance's attributes. It is considered best practice to restrict direct access to the instance attributes, allowing you to manage how they can be updated. In Python, this can be achieved by using the prefix "__" before the attribute name, which makes the attribute private and accessible only within its class through defined methods.
```python
class ALS_RNAseq:
    def __init__(self):
        self.__data_matrix = []
```

By adopting this approach, the attribute "__data_matrix" becomes inaccessible for modification or even reading from outside the class. To enable access to this attribute, it is necessary to develop a getter (or setter) method specifically for this purpose. These methods provide a controlled way of accessing and updating the private attribute, ensuring that any changes to "__data_matrix" are handled appropriately within the class's defined interface.

```python
class ALS_RNAseq:
    def __init__(self):
        self.__data_matrix = []

    def get_data_matrix(self):
        return self.__data_matrix
```

In this instance, the getter method might not seem particularly useful, but adopting this approach is considered good practice. It establishes a foundation for more advanced functionality and maintains the integrity and security of the data within your class.

Think about other preprocessing functions that could be beneficial for subsequent steps. For instance, functions that verify the completeness of required annotations or that allow for the subsetting of your data based on specific annotation criteria (e.g., extracting a sub-dataframe for "control" samples only). Additionally, consider functions that can modify one attribute (such as the data matrix) and update related attributes (like the annotations) accordingly.

Lastly, it might be useful to have a convenient method to "check" or examine your objects. Therefore, defining a "print" function can be particularly helpful. This can be achieved by implementing the "str" method within your class, which alters how the "print" function behaves when applied to your class objects. This method allows you to define a custom string representation of your objects, making it easier to understand their state or contents at a glance.

```python
class ALS_RNAseq:
    def __init__(self):
        self.__data_matrix = []

    def get_data_matrix(self):
        return self.__data_matrix
    
    def __str__(self):
        return "put here what you want to print"
```

# Step 2 - Descriptive Analysis

Descriptive analysis encompasses all kinds of direct data descriptions, such as calculating means, standard deviations, and generating histograms and box plots.

## Sample Description:

For each sample, calculate the mean (across all genes), the median, and the standard deviation. Devise an efficient method to report these data comprehensively. Utilize the annotation data to characterize your entire dataset. Detail the number of "disease groups," the variety of sample "sources," and the count of samples per individual, among other aspects. This overview should assist you in planning the next steps, enabling you to group your data correctly when comparing subsets of samples or to identify and mitigate potential biases.

Your goal is to concisely present all relevant information about the samples, providing a comprehensive yet summarized view of your dataset. To achieve this, consider the use of the following visualizations:

- Bar Charts for "Disease Groups" and Sample "Sources": Create bar charts to show the distribution of samples across different disease groups and sources. This will help you quickly visualize the diversity and spread within your dataset.
- Histograms for Measures of Central Tendency and Dispersion: Use histograms to represent the distribution of means, medians, and standard deviations calculated for each sample. This will provide an overview of the variability in gene expression across your study.
- Summary Tables: Prepare tables to summarize key information, such as the number of samples per disease group, source, and individual. Also include a summary table of descriptive statistics (mean, median, standard deviation) for each sample.
- Box Plots: For a more detailed analysis of gene expression distribution, use box plots for each disease group or sample source. This allows for the visualization of medians, quartiles, and outliers.
- Heatmaps: Although a detailed statistical analysis is not required at this stage, a heatmap showing gene expression across different samples or groups can provide a preliminary visualization of expression patterns.

## RNA Counts Description:

For each gene, calculate the mean (across all samples), the median, and the standard deviation. Develop an effective strategy to report these data and offer your initial interpretations. Note that you may need to employ graphical representations and undertake data transformation or manipulation to achieve a clear understanding of the data.

Given that samples represent different individuals, including ALS patients and control subjects, contemplate the subsets of data that could be analyzed and the rationale for selecting these subsets. Conduct the descriptive analysis for these chosen subsets.

To facilitate this analysis, consider the following steps and visualizations:

- Gene Expression Overview: Start with basic statistics (mean, median, standard deviation) for each gene across all samples. This foundational analysis helps identify genes with high variability, potentially indicating differential expression related to disease status or other factors.
- Data Transformation: Given the likely skewed distribution of gene expression data, log transformation (i.e., just log the counts) or other normalization methods might be necessary to stabilize variance and improve the interpretability of statistical analyses.
- Graphical Representations:
-- Box Plots: Display the distribution of expression levels for each gene across ALS patients and control subjects. Box plots are excellent for visualizing the central tendency, dispersion, and outliers within each group.
-- Histograms: Use histograms to illustrate the distribution of expression levels for selected genes, helping to understand the skewness or bimodality of the data.
-- Heatmaps: Create heatmaps to compare gene expression patterns across samples. This is particularly useful for identifying groups of genes that behave similarly across ALS patients versus controls.

## Begin Your Report:

At this juncture, you should have already started drafting your report. Before proceeding further, take this opportunity to refine your code and enhance your report, ensuring clarity and coherence in your presentation of the analyses conducted so far.

# Step 3 - PCA
As you may have observed, the number of genes is far too high to compare all samples using all genes with simple analyses. 
The PCA is a classical first step analysis in those cases, and offers (among other things) a good way to visualize your data.
To understand what a PCA is, let's check at my favorite youtube channel:
[StatQuest: PCA Step-by-Step.](https://www.youtube.com/watch?v=FgakZw6K1QQ) 
We will review the video together, wait for me please.

To implement a PCA in python, a simple way is to use the PCA function in the sklearn.decomposition package. 
Scikit-learn is a wonderfull Python library, and contains a lot of "must-have" features needed for a data-scientist. 
Take some time to visite the official [website.](https://scikit-learn.org/stable/)
For a pratical python PCA tutorial, let's check again a [Josh Starmer's video](https://www.youtube.com/watch?v=Lsue2gEM9D0&ab_channel=StatQuestwithJoshStarmer).

Before doing a PCA on your data, it is mandatory to "center" your data (so that each gene has a mean equal to zero) and it is advise to "scale" your data (so that each gene has a standard deviation of 1).
This way, you ensure that the PCA is accuratly done (with the centering) and all genes are considered base on their relative variability (with the scaling) and not their absolute value. It is not advise to consider multiple variables with different scale all together.

You can use the following code before doing a PCA (X will be used in the PCA after this code)

```python
from sklearn.preprocessing import StandardScaler
scaler = StandardScaler()
X = scaler.fit_transform(my_data) # my_data being your dataframe containing your genes in column and samples in row
```

Now, perform a PCA and plot the samples using their coordinates in the first PCs. 
TIPs: to select the good number of PCs, compute the percenatage of variance their capture.
Use the annotations to color your plots, and see if you can already observe some kind of signal.

(Bonus) PCA is also good way to find outliers. 
Outliers are samples that are greatly different from the other samples. 
The difference should be "huge", so that only experimental errors could explain it.
Using the PCA and visualization, look at possible outliers.

# Step 4 - tSNE and UMAP (optional)
Another (more recent) good vizualization tool for high dimensional data is the [t-SNE](https://www.youtube.com/watch?v=NEaUSP4YerM&ab_channel=StatQuestwithJoshStarmer), and its little brother, the [UMAP](https://www.youtube.com/watch?v=eN0wFzBA4Sc&t=482s&ab_channel=StatQuestwithJoshStarmer). 
The advantage of this two methods is that they can reduce the dimension of your data using a desired number of components (2 most of the time), not (too much) leaving alway a part of your data variability (in theory). 
On the other hand, they do not preserve large distance interpretation, so that only "local similarities" must be interpreted (e.g., outliers are much more difficult to spot). 
UMAP tends to preserve much better large distances, but still not reach the PCA in this topic.

Try to implement a t-SNE and/or a UMAP. 
UMAP can be implemented using the "umap" module, whereas t-SNE has a scikit-learn implementation. 

Compare this visualition vs the PCA one.

# Step 5 - Univariate Analysis
The next step in our project involves identifying genes whose expression significantly varies between groups of patients with ALS and healthy individuals. To accomplish this, we will utilize the DESeq2 package, a powerful tool designed for differential count data analysis from high-throughput sequencing experiments, such as RNAseq. DESeq2 facilitates this analysis by comparing gene expression levels across samples, accounting for biological variability. Importantly, before conducting the differential analysis, DESeq2 normalizes the data to ensure that comparisons are made on a consistent scale, adjusting for differences in library sizes or sequencing depth across samples.

To begin, it is necessary to install DESeq2 in your Python environment. Although DESeq2 is primarily an R library, Python versions are available, which we will use [PyDESeq2 Github](https://github.com/owkin/PyDESeq2). Ensure that these tools are installed beforehand. Next, please refer to the examples provided on the package's GitHub page to conduct your comparison(s).

### Multiple testing
The p-values give the probability to obtain the observed the difference in mean of the RNA count if there is actually no difference in reality. 
So low p-values can be interpreted as: "unexpected observations". 
However, the more observations we make, the more likely we are to have surprising observations.
So we should account for the "multiple testing" we do.
In general, when exploring data (and not looking for a confirmation), controling the false discovery rate (fdr) is the way to go. 
The fdr "limits" the number of false positive discoveries. 
For example, a fdr of 5% tries to ensure that only 5% of our results will be false positives. 
Note that FDR is automatically calculated by DESeq2.

### Graphical representation
Following the differential expression analysis and identification of significant genes using DESeq2, it's highly beneficial to visualize the results, enhancing the interpretability and communication of our findings. A particularly effective way to do this is through the creation of a volcano plot. A [volcano plot](https://en.wikipedia.org/wiki/Volcano_plot_(statistics)) plot displays the statistical significance of the expression changes (usually as -log10 of the p-value) against the magnitude of change (log2 fold change), allowing for an immediate visual assessment of the data. Genes that are significantly upregulated or downregulated are easily identifiable as they appear to the right and left of the plot, respectively, and further from the bottom, indicating higher levels of significance. This visualization not only helps in quickly spotting genes of potential interest but also in presenting a compelling graphical summary of the differential expression analysis. Creating a volcano plot can be done using various plotting libraries available in Python, such as Matplotlib or Seaborn.

## STEP 6 - Multivariate Analysis - Elastic-Net (Option 1)
DESeq2 is an excellent tool because it offers a way to estimate the relevance of each gene (to distinguish ALS vs control samples for example). 
With this analysis, we are able to rank the genes from the most relevant to the less usefull for our question.
However, this analysis looks at the genes one by one, and doesn't try to take advantage of gene combinations. 

We can extend univariate analysis, using multiple genes simultaneously, with the help of logistic regression. 
The idea behind a logistic regression is to use the RNA counts of a sample as variables (for all genes) and try to predict the sample group (ALS or control for example).
To do so, we try to find a linear combination of the variables that compute a probability for each possible classes (ALS or control).

Regularized logistic regressions are extentions of standard logistic regressions that penalized variables with no clear interest.
In another words, they try to attribute a coefficient of 0 to variables with low interest (no useful signal).
Therefore this methods are very good candidates for biomarker selections. 

Among these regularized regressions, [elastic-net](https://www.youtube.com/watch?v=1dKRdX9bfIo&ab_channel=StatQuestwithJoshStarmer) is one of the most versatile. 
It can handle "wide" dataset (more variables than samples), correlation (do not have to choose between very similar variables) and are easily configurable (parameter tuning).

Warning: like the PCA, the elastic-net algorithm takes multiple variables into account. 
Therefore, variables with high values and high standard deviation can introduce bias (being artifically more important). 
So, like PCA, I advise you to normalize the data before using elastic-net.

(Optional) Also, you can use the output of the PCA and perform the elastic-net on this transformed data.
Interpretations will be more difficult (you will have to look at loadings of your components at the end of the analysis) but potentially more powerfull.

### Elastic-Net implementation
You can use the method "LogisticRegression" from scikitlearn (sklearn.linear_model) to implement an elastic-net analysis:
```python
from sklearn.linear_model import LogisticRegression
elasticNet = LogisticRegression(penalty='elasticnet', solver='saga', l1_ratio=0.5, C=0.5).fit(x, y)
```
Where x contains the data and y the groups (2 groups only).
Two parameters can be fine tune: 
- l1_ratio: the ration between l1 and l2 regularizations
- C: inverse of regularization strength

### Fine tuning of parameters
To fine tune parameters while avoiding overfitting, we can perform a [cross validation](https://www.youtube.com/watch?v=fSytzGwwBVw&ab_channel=StatQuestwithJoshStarmer): 
```python
from sklearn.linear_model import LogisticRegression
elasticNet = LogisticRegressionCV(penalty='elasticnet', cv= 3, solver='saga', l1_ratios=[0.25,0.5,0.75], Cs=[0.1,0.5], scoring= 'accuracy').fit(x, y)
```
You have to specify the number of folds (here cv=3). 
Choose a number of fold adapted to the data. 
By default, the folds are stratified (they split the dataset taking the class balance into account, something that we absolutely want here).

You also have to specify the scoring function used to evaluate and compare the models (accuracy in this example). 
Check the [documentation](https://scikit-learn.org/stable/modules/classes.html#module-sklearn.metrics) to get a list of scoring functions and choose one adapted to the our data.
Have a special look to [balanced accuracy](https://scikit-learn.org/stable/modules/generated/sklearn.metrics.balanced_accuracy_score.html#sklearn.metrics.balanced_accuracy_score), [f1 score](https://scikit-learn.org/stable/modules/generated/sklearn.metrics.f1_score.html#sklearn.metrics.f1_score), [ROC AUC](https://scikit-learn.org/stable/modules/generated/sklearn.metrics.roc_auc_score.html#sklearn.metrics.roc_auc_score), and [the matthews correlation coefficient](https://scikit-learn.org/stable/modules/generated/sklearn.metrics.matthews_corrcoef.html#sklearn.metrics.matthews_corrcoef).

Once you have your elastic-net model (with the best parameters), you can look at its overall performance.
```python
model_prediction = elasticNet.predict(x)
accuracy_train_dataset = sklearn.metrics.accuracy_score(y, model_prediction)
```
Here I used the accuracy metric, but remember to use (at least) the metric you used during the tuning.

### Testing your model
We built and evaluated our model on the same dataset. 
Even if we used a cross-validation to choose the parameters, the final model used all the avalaible data. 
Therefore we have no idea how the model performs on new data. 
For this purpose, we need a "test set". 

Download the "data set", process the data (warning: there might be some differences in the xml file) and use this new dataset to evaluate the model.
I advise you to use different scores (begining with the one used to choose the parameters).

### Get the best genes
The elastic-net model contains (by definition) all the coefficients of the regression, i.e., the importance of each gene to differentiate the ALS samples from the control samples.
```python
elasticNet.coef_
```
You can rank these coefficients (using their absolute values) to obtain a good list of candidates.

## STEP 7 - XGBOOST (Option 2)
There are many other methods to analyze multivariate data, and some of them are more or less equivalent. 
A different class of methods compare to the regularization regression is the class of decision tree based methods. 
Among this class, the "ensemble" classes, i.e., using a set of decision trees, are very used nowadays because they are very powerful.
In particular, the xgboost method is extremely popular and used in pretty much all data analysis competitions.
For a detailled explanation refer to this link [xgboost explanation](https://www.youtube.com/watch?v=OtD8wVaFm6E&ab_channel=StatQuestwithJoshStarmer).

The xgboost implementation for biomarker selection is pretty straightforward using the package xgboost: 
```python
from xgboost import XGBClassifier
model = XGBClassifier()
model.fit(X_train, y) # xtrain is your data matrix, and y is your group labels (ctrl or ALS)
```
Then you can extract the feature (=gene) importance of your model:
```python
print(model.feature_importances_)
```
You can also used the built-in funtion to plot the feature importance: 
```python
from xgboost import plot_importance
plot_importance(model)
```
As for the elastic-net, remember to:
- Standardized your data
- Test your model on an independant test data set (there is no fine tuning in the previous code, so cross validation is not mandatory)

## STEP 8 (final)- Report your results
At this step, you should have different lists of candidate biomarkers. 
Propose a final ranked list, making use of all your results, containing the TOP "100" candidates.
This list will be compared to the state of art of ALS biomarkers but also other related diseases.
The state of art lists will be avalaible in this repository after the 4 may. 

Good luck !
