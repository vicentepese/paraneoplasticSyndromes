# Genetic and HLA study of paraneoplastic syndromes 

## Use permit and author declaration 
This project was created under Stanford School of Medicine's Sleep Unit with Dr. Emmanuel Mignot as Principal Investigator. Any data belonging to the project, including but not limited to genetic markers, HLA imputations, or final diagnoses have been removed from the repository to protect patients and comply with data usage regulation. No Patient Health Information protected under HIPPA regulation was included in the project. The author discourages the use of the code for other purposes but the one intended in this project but discloses the methods utilized to serve as an inspiration for other researchers.  

## Introduction 

Paraneoplastic syndromes (PS) are a group of rare disorders triggered by an abnormal immune response to a tumor or "neoplasm". PSs are hypothesized to occur when cancer-fighting cells antibodies or white blood cells (T-cells) attack healthy cells in the neural system in a process called cross-reactivity. This process takes place when cell surface proteins that are naturally found in other parts of the body (the neural system) are synthesized by the tumor, activating an immune response.

In this project, we studied the genetic profile of anti-LGI1 encephalitis, a subtype of limbic encephalitis (LE) characterized by the presence of LGI1 antibodies in the Cerebral Spinal Fluid (CSF). It is characterized by a myriad of clinical presentations and symptom specificity including seizures, memory problems, irritability, depression, confusion, and dementia. Typically, LE was considered a PS but new antibodies unassociated with tumors in LE patients suggest the possibility of an autoimmune explanation to this condition.

The principal aim of this study was to:
1. Identify Single Nucleotide Polymorphisms (SNP) associated with LGI1 LE.
2. Identify Human Leukocyte Antigen (HLA) associated with LGI1 LE.
3. Identify specific Amino Acids (AA) associated with LGIE LE.

# Genome-Wide Association Study

Genome-Wide Association Studies (GWAS) is an observational study of a genome-wide set of genetic variants in a control versus case paradigm to identify variants associated with a trait (i.e. a condition in cases). If one of the variants (one allele) is significantly more present in cases it is said to be *associated* with the disease. The associated SNPs mark a region in the genome that may influence the risk of the disease. 

We studied the associated genetic variants of LGI1 LE by performing a GWAS. The code belonging to this part of the study is located in the directory *GWAS*. For more information, please refer to the directory. 

# Human Leukocyte Antigen analysis

The Human Leukocyte Antigen System (HLA) is a group of proteins encoded by the Major Histocompatibility Complex (MHC) gene complex in humans located in the 6<sup>th</sup> chromosome.

<p align="center">
  <img src="https://upload.wikimedia.org/wikipedia/commons/7/77/HLA.svg">
</p>


These cell-surface proteins are responsible for the regulation of the immune system and for triggering immune responses when a foreign organism is encountered. There are three different classes of MHC encoded by multiple regions of HLA with different functions:
1. MHC Class I: encoded in HLA-A, B, and C. These proteins are responsible for presenting peptides from inside the cell. For instance, if the cell is infected by a foreign organism (i.e. a virus) the HLA system will bring fragments of the virus to the surface of the cell so that it is identified by the immune system as an infected cell and is destroyed by killer T-cells.
2. MHC Class II: encoded in HLA-DP, DM, DO, DQ and DR. These proteins present antigens from outside the cell to T-cells, and stimulate the multiplication of T-helper cells which in turn stimulate antibody-producing B-cells to produce antibodies to that specific antigen.
3. MHC Class III: encodes components of the complement system. 

The HLA system is highly polymorphic - that is, there are many different types of alleles. For this reason, we identified which alleles are associated with the disease. In addition, we narrowed down to Amino Acid association with LGI1 LE and studied other aspects related to hetero- and homozygosity.he code belonging to this part of the study is located in the directory *HLA_analysis*. For more information, please refer to the directory. 
