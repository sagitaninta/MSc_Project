### Appendix of code
Link to MSc project github repository: (https://github.com/DebGreer/MSc_Project)  
Links to code for all sections of the pipeline and producing graphs and statistics are included below: 

#### Creating reference files
[Identifying ancestral alleles script](https://github.com/DebGreer/MSc_Project/blob/master/anc_seq_v1.py)

[Creating a conservation score reference file explainer/command summary](https://github.com/DebGreer/MSc_Project/blob/master/Creating%20a%20conservation%20score%20reference%20file.md)  
[R liftover script](https://github.com/DebGreer/MSc_Project/blob/master/liftOver_apocrita.R)

[Parallel script for running sift score script for each chromosome](https://github.com/DebGreer/MSc_Project/blob/master/parallel_sift_all.sh)
[Creating SIFT score reference file script](https://github.com/DebGreer/MSc_Project/blob/master/sift_ref_file_v3.sh)

#### Per sample pipeline
Obtaining posterior probabilities:
[Getting posterior probabilities from gentoype likelihoods script](https://github.com/DebGreer/MSc_Project/blob/master/all_genotype_likelihoods_v3.py)  

Calculating load scores from BAM file and reference score file:  
[Conservation load score script](https://github.com/DebGreer/MSc_Project/blob/master/mut_load_calculator_v3.py)

[SIFT load score script](https://github.com/DebGreer/MSc_Project/blob/master/sift_calculator_v2.py)

Overall pipeline:
[Example array script](https://github.com/DebGreer/MSc_Project/blob/master/array_pipeline_combined_modern_set4_1.sh)

#### Statistical analyses, figures and metadata

[R analysis script - plots, stats, sample map](https://github.com/DebGreer/MSc_Project/blob/master/Analysis_MSc_v6.R) 

[Coverage calculation](https://github.com/DebGreer/MSc_Project/blob/master/coverage_calculator.sh)

[Heterozygosity graph](https://github.com/DebGreer/MSc_Project/blob/master/Heterozygosity_graph_v3.R)

[Prior comparison](https://github.com/DebGreer/MSc_Project/blob/master/Prior_comparison_v2.R)
