

# Multilayer-Network-Centrality
This repository contains scripts to calculate centrality of genes in a multi-tissue network.

<object data="https://github.com/tarunkumariitm/Multilayer-Network-Centrality/blob/main/MultiCens.pdf" type="application/pdf" width="700px" height="700px">
    <embed src="https://github.com/tarunkumariitm/Multilayer-Network-Centrality/blob/main/MultiCens.pdf">
        <p>This browser does not support PDFs. Please download the PDF to view it: <a href="http://yoursite.com/the.pdf">Download PDF</a>.</p>
    </embed>
</object>


### Step 1: Generate coexpression based supra-adjacency matrix from multitissue data
Run process_GTEx_v8_eQTL.R to generate a supra-adjacency matrices for tissue-tissue pairs. In the R script, edit the hormones list and keep the hormones for which the supra-adjacency matrix is required.
The script uses GTEx v8 dataset to build co-expression network


### Step 2: Calculate gene centrality

Before executing the code, make sure you have all the packages from requirements.txt

Open the following files in a jupyter-notebook and execute all the cells sequentially. 

1. insulin_source_genes.ipynb: Finds centrality of all genes in the source tissue. For example, in case of Insulin, this script finds centrality of all the genes in Pancreas tissue. The script uses supra-adjacency matrix from step 1 and requires insulin-responding genes in the Muscle tissue.
2. insulin_target_genes.ipynb: Finds centrality of all genes in the target tissue. For example, in case of Insulin, this script finds centrality of all the genes in the Muscle tissue. The script uses supra-adjacency matrix from step 1 and requires insulin-producing genes in the Pancreas tissue.
