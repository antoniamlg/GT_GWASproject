# GWAS Pipeline Documentation

- [GWAS Pipeline Documentation](#gwas-pipeline-documentation)
  - [Steps Overview](#steps-overview)
  - [Datanumbers](#datanumbers)
- [0. Filter by phenotype \[X\]](#0-filter-by-phenotype-x)
- [1. Initial per-chip QC (split) \[X\]](#1-initial-per-chip-qc-split-x)
      - [As we do not have that many individuals (now ~`1071`), I will start with the SNP-level QC first, to preserve as many individuals as I can.](#as-we-do-not-have-that-many-individuals-now-1071-i-will-start-with-the-snp-level-qc-first-to-preserve-as-many-individuals-as-i-can)
  - [1.1 SNP-level QC \[X\]](#11-snp-level-qc-x)
    - [1.1.1 Call rate/Missingness \[ \]](#111-call-ratemissingness--)
      - [Results](#results)
    - [1.1.2 Filter SNPs that are missing in 100% of individuals + with 0.05 threshold \[X\]](#112-filter-snps-that-are-missing-in-100-of-individuals--with-005-threshold-x)
    - [1.1.3 Replot missingness to see if i actually did something. \[X\]](#113-replot-missingness-to-see-if-i-actually-did-something-x)
    - [1.1.4 Merge chips \[X\]](#114-merge-chips-x)
    - [1.1.5 Sex-imputation (TRY) \[X\]](#115-sex-imputation-try-x)
    - [1.1.5 MAF \& HWE \[X\]](#115-maf--hwe-x)
  - [1.2 Sample-level QC \[ \]](#12-sample-level-qc--)
    - [1.2.1 Sex checks \[ \]](#121-sex-checks--)
    - [1.2.2 Missingness \[ \]](#122-missingness--)
    - [1.2.3 Heterozygosity outliers \[ \]](#123-heterozygosity-outliers--)
      - [Decide on thresholds](#decide-on-thresholds)
    - [1.2.4 Relatedness / Duplicates \[ \]](#124-relatedness--duplicates--)
  - [1.3 Population Structure/ Stratification/ Batch Effects](#13-population-structure-stratification-batch-effects)
    - [1.3.1 PCA](#131-pca)
          - [After QC, merge chips again](#after-qc-merge-chips-again)
  - [1. Impute Sex](#1-impute-sex)
  - [2. Check Sex](#2-check-sex)
  - [How do do .md files](#how-do-do-md-files)
    - [code blocks](#code-blocks)
    - [bullet lists](#bullet-lists)
    - [tables](#tables)
    - [links](#links)
    - [images](#images)
    - [text formatting](#text-formatting)

## Steps Overview
| Step         | Command                              | Description                  |
|--------------|--------------------------------------|------------------------------|
| Split IID by chip   | splitbychip.ipynb (Python)    | extract IIDs specific to one chip |
|    |                      |          |

Also check out this link [QC README](https://github.com/kaspermunch/PopulationGenomicsCourse/blob/master/Exercises/GWAS_QC/step_by_step_tutorial.md).

## Datanumbers
| Chip Name             | # Individuals | # Variants          | #Ind. after phenotype filtering  | #Var. after SNP QC (hwe + maf) | # Ind. after Sample QC
|-----------------------|---------------|---------------------|----------------------------------|------|-----|
| whole dataset         | 2009          | 1376653             | 1071                             | 
| HTS_iSelect_HD        | 587           |                     | 275                              | 513670 | 264 |
| Illumina_GSAs         | 248           |                     | 128                              | 410086 | 121 |
| OmniExpress           | 291           |                     | 170 -> 133 after sex-check       | 569238 | 125 |
| OmniExpress_plus      | 484           |                     | 266                              | 84250 | 258 |
| unknown_chip          | 399           |                     | 232                              | 17 | - (no heterozygosity check) |

# 0. Filter by phenotype [X]
Our height.txt file only contains the phenotype of `1106` individuals. We cannot use the remaining ones without a phenotype for anything so before doing any QC, splitting etc. filter them out.

```bash
# create keep file from individuals
# duplicate column one to get the following format: IID FID HEIGHT
awk '{print $1, $1}' height.txt > /faststorage/project/populationgenomics/students/amlg/project/data/shit_phenotypes_kickedout/phenokeep.txt

# then use plink to filter the b-files
plink --bfile shit_phenotypes_kickedout/pheno_filtered --keep splitbychip/chipX.keep --make-bed --out splitbychip/chipX
```

* --keep: 1071 people remaining
*  At least 2 duplicate IDs in --keep file
*  1376653 variants and 1071 people pass filters and QC

# 1. Initial per-chip QC (split) [X]
- why: different chips = different SNPs, genotyping errors etc.

Python script [find script](https://github.com/antoniamlg/GT_GWASproject/blob/main/Scripts/splitbychip.ipynb) to:
1. Read the metadata and split IIDs by chip.
2. Create `.keep` files for each chip group.

Use `.keep` files to filter bed file into chip specific files:
```bash
plink --bfile gwas_data --keep /faststorage/project/populationgenomics/students/amlg/project/data/splitbychip/chipX.keep --make-bed --out /faststorage/project/populationgenomics/students/amlg/project/data/splitbychip/chipX_data
```

This results in 6 different b-file sets. Each of them only containing data specific to the chip. I decided, to not use the all the datapoints with the "chip" label as they do not contain any useful data.

#### As we do not have that many individuals (now ~`1071`), I will start with the SNP-level QC first, to preserve as many individuals as I can. ####

## 1.1 SNP-level QC [X]

### 1.1.1 Call rate/Missingness [ ]
Run the --missing command to generate the GWA-data.lmiss with the missing data rate for each SNP. <br>
Use R to make a histogram of the missing data rates (F_MISS).
Run the test-missing command and make a list of all the names of all SNPs where the differential missingness p-value is less than 1e-5. Save the list as fail-diffmiss-qc.txt.

```bash
plink --bfile splitbychip/unknown_chip_data --missing --out ../SNPmiss/unknown_chip_SNPmiss
ls -lht
```

See visualization script [here](https://github.com/antoniamlg/GT_GWASproject/blob/main/Scripts/SNPmiss_QC.ipynb).
[X] Before continuing to the next step, do a plot to see where to set the threshold. Apparently, 0.05 removes almost everything which is annoying.

A lot of care must be taken when filtering variants, since we can lose potential variants associated with the phenotypic trait.

#### Results

* use `0.9` as a filtering threshold to first get rid of SNPs which are missing in all individuals (data cleaning)
* use `0.05` -> If a SNP has a high F_MISS (e.g., > 0.05), it means that a large proportion of individuals have missing data for that SNP.
* the `unknown_chip` does not have any data below 0.2% missing data rate. So the idea is, to just exclude the whole chip. Does that make any sense though? <br>
=> sooo next step: [X] do the actual filtering

### 1.1.2 Filter SNPs that are missing in 100% of individuals + with 0.05 threshold [X]

```bash
plink --bfile chipX --geno 0.05 --make-bed --out chipX...
```

=> to not do it on every single dataset [find script here](https://github.com/antoniamlg/GT_GWASproject/blob/main/Scripts/QC_geno_filter.sh)

### 1.1.3 Replot missingness to see if i actually did something. [X]
I plotted the datasets after the filtering to check if everything worked (it worked)

### 1.1.4 Merge chips [X]
We need to run the sex-imputation on all the datasets, so we at least have a chance to impute sex for the `omics` chip which does not have any information about the sex. <br>
Sooooooooo, apparently you cannot really merge bfiles back together, if they don't have exactly the same SNPs

### 1.1.5 Sex-imputation (TRY) [X]
If it works - amazing.
If it doesn't - just continue with QC and don't care too much (-*Bjarke*).

I tried the 
```bash
plink --bfile OmniExpress_geno005 --check-sex --out OmniExpress_sexcheck
```
and it outputs
```bash
(popgen) [amlg@s21n34 SNPmiss005]$ plink --bfile OmniExpress_geno005 --check-sex --out OmniExpress_sexcheck
PLINK v1.90b6.21 64-bit (19 Oct 2020)          www.cog-genomics.org/plink/1.9/
(C) 2005-2020 Shaun Purcell, Christopher Chang   GNU General Public License v3
Logging to OmniExpress_sexcheck.log.
Options in effect:
  --bfile OmniExpress_geno005
  --check-sex
  --out OmniExpress_sexcheck

515438 MB RAM detected; reserving 257719 MB for main workspace.
625388 variants loaded from .bim file.
170 people (62 males, 71 females, 37 ambiguous) loaded from .fam.
Ambiguous sex IDs written to OmniExpress_sexcheck.nosex .
Using 1 thread (no multithreaded calculations invoked).
Before main variant filters, 170 founders and 0 nonfounders present.
Calculating allele frequencies... done.
Total genotyping rate is 0.971161.
625388 variants and 170 people pass filters and QC.
Note: No phenotypes present.
Error: --check-sex/--impute-sex requires at least one polymorphic X chromosome
locus.
```

This means
* none of the SNPs on the X chromosome are present in your dataset and variable (polymorphic)
* you cannot perform the `plink sex-check`
  
Check for errors in the Omni_Express dataset:
1. Check if you have X chromosome SNPs: 
```bash
   awk '$1 == 23' OmniExpress_geno005.bim | wc -l
```
   <br> returns 0, so the dataset does not contain any X-chromosome SNPs - which means **sex imputation is not possible from this dataset**

2. Keep sex as-is and exclude ambigous:
   Dataset has 
* 62 males
* 71 females
* 37 ambiguous -> fam file has sex info for most of individuals
   So we just remove them.

```bash
   plink --bfile OmniExpress_geno005 --remove OmniExpress_sexcheck.nosex --make-bed --out OmniExpress_clean
```


### 1.1.5 MAF & HWE [X]

After doing all this, remove all low-quality SNPs
```bash
plink --bfile OmniExpress_clean \
      --maf 0.01 \
      --hwe 1e-6 \
      --make-bed \
      --out OmniExpress_qc
```

In addition to removing SNPs identified with differential call rates between cases and controls, this command removes SNPs with call rate less than 95% with --geno option and deviation from HWE (p<1e-6) with the --hwe option. It also removes all SNPs with minor allele frequency less than a specified threshold using the --maf option.

## 1.2 Sample-level QC [ ]

### 1.2.1 Sex checks [ ]

Don't use for *Omics*-chip - it does not have any info on sex

[?] Is the number of individuals on the Omics chip == ambiguous individuals? Or do we have more individuals somewhere without any sex? -> also some on the other chips

```bash
plink --bfile unknown_chip_snpqc --check-sex --out ../SampleQC_sexcheck/unknown_chip_samplesex
```
Only the unknown-chip data produced a nosex-output file. So we will filter these one out but then we are also done with the sex-check :)

Remove proplematic sex with
```bash
plink --bfile unknown_chip_snpqc --remove ../SampleQC_sexcheck/unknown_chip_samplesex.nosex --make-bed --out unknown_chip_snpqc_cleanedsex
```

### 1.2.2 Missingness [ ]

Look at section 1.1.1, we already did that there.
Calculated it again anyways to be sure.
```bash
plink --bfile GWA-QC --missing --out GWA-QC
```

### 1.2.3 Heterozygosity outliers [ ]

Heterozygosity helps you:
* Detect individuals with unusual genotypes (e.g. contamination, sample swaps, inbreeding)
* Combine with missingness for sample QC
* Identify outliers that could skew your results

```bash
plink --bfile GWA-QC --het --out GWA-QC 
```
-> failed for the unknown-chip data: Error: --het requires at least one polymorphic autosomal marker
` awk '{print $1}' unknown_chip_snpqc_cleanedsex.bim | sort | uniq -
```bash
# check if we even have autosomes
 awk '{print $1}' unknown_chip_snpqc_cleanedsex.bim | sort | uniq -c
# ~out: 17 24
```
This means, the unknown chip data contains 17x a SNP mapping to chromosome 24. Chromosome 24 corresponds to the Y chromosome. So we have no polymorphic autosomal SNPs in this dataset.

calculate the observed heterozygosity rate per individual using the formula:

Het = (N(NM) − O(Hom))/N(NM) -> look at the exercises!

* then, do a plot of *F-values* and find outliers manually in Python/R
* find plot [here](https://github.com/antoniamlg/GT_GWASproject/blob/main/Scripts/heterozygosityVSmissingness.ipynb) 

Each point is one sample. We have 802 samples left.
* X-axis: heterozygosity rate
* Y-axis: missing genotype rate
* colors: different genotyping chips

Normal samples should cluster together with **low missing rate** aka y~0. Reasonable heterozygosity lies around ~0.30 ± 0.05. <br>
Outlier samples to check:
* High missingness: points with high y-value
* High heterozygosity: may indicate sample contamination or mix-ups (x > 0.35-0.4) 
* low heterozygosity: could suggest inbreeding or poor genotyping (x < 0.25) -> shit, that's the whole Illumina chip :(
* far from cluster: anything not clustering with the rest of the group

So the points above 0.2 missingness rate are problematic and also the one green point at (1,1), which has 100% missingness and very high heterozygosity.

#### Decide on thresholds ####
* missingness > 0.05 (already did that but anyways)
* remove samples more than 3 standard deviations from the mean

Make a file with the FID and IID of all individuals that have a genotype missing rate >=0.03 or a heterozygosity rate that is more than 3 s.d. from the mean. Then use plink to remove these individuals from the data set.
(Exercise suggestion) <br>

The code for how these .remove files were created lies [here](https://github.com/antoniamlg/GT_GWASproject/blob/main/Scripts/heterozygosityVSmissingness.ipynb).
I will now remove these outliers from the files.

```bash
plink --bfile GWA-QC --remove wrong_het_missing_values.txt --make-bed --out GWA-QC
```

### 1.2.4 Relatedness / Duplicates [ ]

* what is this command? flags individuals with pi_hat > 0.185
```bash
plink --bfile chipX_qc4 --genome --min 0.185 --out chipX_related

```

1. prune the data
```bash
plink --bfile GWA-QC --indep-pairwise 500kb 5 0.2 --out GWA-QC
```
2. Calculate IBD between each pair of individuals
```bash
plink --bfile GWA-QC --extract GWA-QC.prune.in --genome --min 0.185 --out GWA-QC
```
3. Remove a member from each of the pairs that are too closely related from the data set. To keep it simple you can just always remove the individual mentioned first. *

4. Remove individuals with wrong IBD.
```bash
plink --bfile  GWA-QC --remove wrong_ibd.txt --make-bed --out GWA-QC
```
## 1.3 Population Structure/ Stratification/ Batch Effects

### 1.3.1 PCA
```bash
plink --bfile chipX_qc4 --pca 10 --out chipX_pca
```

###### After QC, merge chips again ######
---
## 1. Impute Sex
```bash
plink --bfile gwas_data_indiv_filtered
      --impute-sex 0.8 
      --make-bed --out gwas_data_sex_imputed
      --allow-no-sex
```
| command                 | what it does                                                                      |
|-------------------------|-----------------------------------------------------------------------------------|
| ```--bfile```           | reads in all .bed, .bam & .fam files                                              |
| ```--impute-sex 0.8```  | command to impute sex, F>0.8 likely male, F<0.2 likely female, 0.8 standard value |
| ```--make-bed```        | needed to produce new output files                                                |
| ```--allow-no-sex```    | allows PLINK to process those individuals at all                                  |

- impute sex = infer sex: change sex assignments to imputed values
- get sex of individuals based on X chromosome heterozygosity rates
  - SNPs determine whether an individual is likely male (low heterozygosity on X, since XY) or female (normal heterozygosity on X, since XX)
- outputs a .sexcheck file

```bash
# The terminal output I get
515439 MB RAM detected; reserving 257719 MB for main workspace.
1376653 variants loaded from .bim file.
2009 people (992 males, 916 females, 101 ambiguous) loaded from .fam.
Ambiguous sex IDs written to
../../students/amlg/project/data/gwas_data_sex_imputed.nosex .
Using 1 thread (no multithreaded calculations invoked).
Before main variant filters, 2009 founders and 0 nonfounders present.
Calculating allele frequencies... done.
Warning: 106246 het. haploid genotypes present (see
../../students/amlg/project/data/gwas_data_sex_imputed.hh ); many commands
treat these as missing.
Warning: Nonmissing nonmale Y chromosome genotype(s) present; many commands
treat these as missing.
Total genotyping rate is 0.478922.
1376653 variants and 2009 people pass filters and QC.
Note: No phenotypes present.
--impute-sex: 29856 Xchr and 0 Ychr variant(s) scanned, 1889/2009 sexes
imputed. Report written to
../../students/amlg/project/data/gwas_data_sex_imputed.sexcheck .
--make-bed to ../../students/amlg/project/data/gwas_data_sex_imputed.bed +
../../students/amlg/project/data/gwas_data_sex_imputed.bim +
../../students/amlg/project/data/gwas_data_sex_imputed.fam ... done.
```

## 2. Check Sex
---

## How do do .md files

### code blocks
Use triple backticks (```) for code blocks:

```bash
plink --bfile data --mind 0.1 --geno 0.05 --maf 0.01 --make-bed --out qc_filtered
```
You can specify language after the triple backticks (e.g., bash, r, python) for syntax highlighting.

Also possible like this
```bash
plink --bfile data/raw/genotypes \
  --mind 0.1 \
  --geno 0.05 \
  --maf 0.01 \
  --make-bed \
  --out data/clean/genotypes_qc
```

### bullet lists
with - or 1.

### tables
| Step         | Command                              | Description                  |
|--------------|--------------------------------------|------------------------------|
| Initial QC   | plink --mind 0.1                     | Remove individuals with >10% missing data |
| MAF Filter   | plink --maf 0.01                     | Remove rare variants         |

### links
[link name](https://github.com/antoniamlg/GT_GWASproject)

### images
![image name](image.png)

### text formatting
- use (*) or (_) for *italic* or _italic_ text
- use (**) or (__) for __bold__ text