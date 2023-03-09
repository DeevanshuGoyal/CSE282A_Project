# CSE282A_Project - Diet Optimization for Dysbiotic Microbiome

Diet has been shown to affect the composition of microbiome by changing the number and activity
of different microbial communities by providing or depleting the necessary nutrients. This impor-
tant role suggests that a modified and personalized diet could be utilized to help re-establish a
healthy microbiome from a diseased state. Given the heterogeneity of these communities, and their
differences across individuals and time, this approach presents two main challenges:

* Defining the healthy vs. disease state microbiome profile and; 
* Finding an optimal diet that can most efficiently restore the healthy gut microbiome. 

In this project, we aim to address the second challenge.

# Algorithmic approach

In this project, we implement three algorithmic approaches to find the set of nutrients that can most efficiently restore the gut microbiome diversity (as per the reward function defined in the report). These approaches are Naive randomised selection, Gibbs' sampling, and Randomised divide and conquer, which are implemented as laid out in [CSE282A_project.py](https://github.com/DeevanshuGoyal/CSE282A_Project/blob/main/CSE282A_project.py).

# Repository guide

* [data](https://github.com/DeevanshuGoyal/CSE282A_Project/tree/main/data): This directory contains the raw data for the human gut microbiome study
  * [Data description_400.docx](https://github.com/DeevanshuGoyal/CSE282A_Project/blob/main/data/Data%20description_400.docx): Contains details on what data do the different files represent alongside additional instructions on how to process it
  * [taxonomy_400.csv](https://github.com/DeevanshuGoyal/CSE282A_Project/blob/main/data/taxonomy_400.csv): Provides the relative abudance(RA) % data for all 2,827 samples (columns) and the total of 400 taxons (or ASVs) (rows)
  * [metadata.csv](https://github.com/DeevanshuGoyal/CSE282A_Project/blob/main/data/metadata.csv): Contains information on the classification of a particular sample into 'NORMAL', 'DEVIANT', and 'OTHERS' categories. 'NORMAL' samples constitute the reference dataset for a healthy gut microbiome population. Each of the 'DEVIANT' samples constitutes an unhealthy and imbalanced gut microbiome population.
  * [nim-aminoacidsD_400.csv](https://github.com/DeevanshuGoyal/CSE282A_Project/blob/main/data/nim-aminoacidsD_400.csv): Provides the nutitional impact of the nutrients in this table (amino acids as a source of carbon, nitrogen and energy) for each of the 400 taxons (or ASVs)
  * [nim-aminoacids_400.csv](https://github.com/DeevanshuGoyal/CSE282A_Project/blob/main/data/nim-aminoacids_400.csv): Provides the nutitional impact of the nutrients in this table (amino acids as protein building blocks) for each of the 400 taxons (or ASVs)
  * [nim-sugars_400.csv](https://github.com/DeevanshuGoyal/CSE282A_Project/blob/main/data/nim-sugars_400.csv): Provides the nutitional impact of the nutrients in this table (carbohydrates as a source of carbon and energy) for each of the 400 taxons (or ASVs)
  * [nim-vitamins_400.csv](https://github.com/DeevanshuGoyal/CSE282A_Project/blob/main/data/nim-vitamins_400.csv): Provides the nutitional impact of the nutrients in this table (vitamins as precursors of essential metabolic cofactors) for each of the 400 taxons (or ASVs)
  * [Phenotypes_legend.xlsx](https://github.com/DeevanshuGoyal/CSE282A_Project/blob/main/data/Phenotypes_legend.xlsx): contains explanation of binary phenotypes and abbreviations for all nutrients considered in this study
