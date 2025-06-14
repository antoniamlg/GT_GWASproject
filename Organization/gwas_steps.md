# GWAS Pipeline Documentation

- [GWAS Pipeline Documentation](#gwas-pipeline-documentation)
  - [Steps Overview](#steps-overview)
  - [Datanumbers](#datanumbers)
- [0. Filter by phenotype \[X\]](#0-filter-by-phenotype-x)
- [1. Initial per-chip QC (split) \[X\]](#1-initial-per-chip-qc-split-x)
      - [As we do not have that many individuals (now ~`1071`), I will start with the SNP-level QC first, to preserve as many individuals as I can.](#as-we-do-not-have-that-many-individuals-now-1071-i-will-start-with-the-snp-level-qc-first-to-preserve-as-many-individuals-as-i-can)
  - [1.1 SNP-level QC \[X\]](#11-snp-level-qc-x)
    - [1.1.1 Call rate/Missingness \[X\]](#111-call-ratemissingness-x)
      - [Results](#results)
    - [1.1.2 Filter SNPs that are missing in 100% of individuals + with 0.05 threshold \[X\]](#112-filter-snps-that-are-missing-in-100-of-individuals--with-005-threshold-x)
    - [1.1.3 Replot missingness to see if i actually did something. \[X\]](#113-replot-missingness-to-see-if-i-actually-did-something-x)
    - [1.1.4 Merge chips \[X\]](#114-merge-chips-x)
    - [1.1.5 Sex-imputation (TRY) \[X\]](#115-sex-imputation-try-x)
    - [1.1.5 MAF \& HWE \[X\]](#115-maf--hwe-x)
  - [1.2 Sample-level QC \[X\]](#12-sample-level-qc-x)
    - [1.2.1 Sex checks \[X\]](#121-sex-checks-x)
    - [1.2.2 Missingness \[X\]](#122-missingness-x)
    - [1.2.3 Heterozygosity outliers \[X\]](#123-heterozygosity-outliers-x)
      - [Decide on thresholds](#decide-on-thresholds)
    - [1.2.4 Relatedness / Duplicates \[X\]](#124-relatedness--duplicates-x)
      - [Why check relatedness](#why-check-relatedness)
      - [Merging the datasets #### \[X\]](#merging-the-datasets--x)
- [2. PCA \[X\]](#2-pca-x)
      - [I also redid the plot with sex-colouring](#i-also-redid-the-plot-with-sex-colouring)
- [3. GWAS on height](#3-gwas-on-height)
  - [3.1 Prepare phenotype file \[X\]](#31-prepare-phenotype-file-x)
  - [3.2 Run the GWAS](#32-run-the-gwas)
  - [3.3 Count significant SNPs](#33-count-significant-snps)
  - [3.4 Manhattan/QQ Plots](#34-manhattanqq-plots)
  - [3.5 QQ-plot](#35-qq-plot)
  - [3.6 Genomic Inflation Factor](#36-genomic-inflation-factor)
  - [3.7 Genomic Control \& Most Significant p After GC](#37-genomic-control--most-significant-p-after-gc)
  - [3.8 Replot](#38-replot)
- [4. PCA II + Adjust for PCs](#4-pca-ii--adjust-for-pcs)
  - [4.1 Use PCs as covariates in GWAS](#41-use-pcs-as-covariates-in-gwas)
- [5. Do additional analyses using Plink, GCTA, R or any other tool you might find relevant](#5-do-additional-analyses-using-plink-gcta-r-or-any-other-tool-you-might-find-relevant)
  - [5.1 Association test conditioning on the most significant SNP (Option 8)](#51-association-test-conditioning-on-the-most-significant-snp-option-8)
    - [5.1.1 Identify top SNP](#511-identify-top-snp)
    - [5.1.2 Run conditional analysis](#512-run-conditional-analysis)
  - [5.2 Distribution of phenotypes by genotype at the most significant SNP (Option 6)](#52-distribution-of-phenotypes-by-genotype-at-the-most-significant-snp-option-6)
  - [5.2.1 Extract genotype for one SNP with PLINK](#521-extract-genotype-for-one-snp-with-plink)

## Steps Overview
| Step         | Command                              | Description                  |
|--------------|--------------------------------------|------------------------------|
| Split IID by chip   | splitbychip.ipynb (Python)    | extract IIDs specific to one chip |
|    |                      |          |

Also check out this link [QC README](https://github.com/kaspermunch/PopulationGenomicsCourse/blob/master/Exercises/GWAS_QC/step_by_step_tutorial.md).

## Datanumbers
| Chip Name             | # Individuals | # Variants          | #Ind. after phenotype filtering  | #Var. after SNP QC (hwe + maf) | # Ind. after Sample QC
|-----------------------|---------------|---------------------|----------------------------------|------|-----|
| whole dataset         | 2009          | 1376653             | 1071                             | 768 |
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

### 1.1.1 Call rate/Missingness [X]
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

## 1.2 Sample-level QC [X]

### 1.2.1 Sex checks [X]

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

### 1.2.2 Missingness [X]

Look at section 1.1.1, we already did that there.
Calculated it again anyways to be sure.
```bash
plink --bfile GWA-QC --missing --out GWA-QC
```

### 1.2.3 Heterozygosity outliers [X]

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

### 1.2.4 Relatedness / Duplicates [X]
Yay, after a loooooong way through the QC we finally reached the first question (*pat on the back) <br>
#### Why check relatedness ####
Related individuals share segments of genome - identical by descent (IBD). If there are too many related individuals present, it can bias your association results. So you want to detect pairs with high relatedness and remove one from each pair. <br>
In the [exercises](https://github.com/kaspermunch/PopulationGenomicsCourse/tree/master/Exercises/GWAS_QC) we got the following pipeline.

#### Merging the datasets #### [X]
1. Choose base dataset
2. Create merge list (txt)
3. Merge: `plink --bfile OmniExpress_plus_snpqc --merge-list merge_list.txt --make-bed --out merged_chips` <br>
  => 1117723 variants and 768 people pass filters and QC.

---
1. Flagging pairs with $\hat{\pi}$ > 0.185
```bash
plink --bfile merged_afterqc --genome --min 0.185 --out related/merged_afterqc
```
* `--genome` : calculate pairwise genome-wide IBD estimates for all pairs of individuals
* `--min 0.185` : only output pairs with PI_HAT > 0.185 (a measure of relatedness)
  * $\hat{\pi}$ ~ 0.5: dupicates or twins
  * $\hat{\pi}$ ~ 0.25: 1st degree relatives (parent-child, siblings)
  * $\hat{\pi}$ ~ 0.125: 2nd degree relatives (grandparent, uncle)
  * `0.185` is a common threshold to flag close relatives
* Result: a .genome file listing pairs of individuals with estimated relatedness > 0.185.

2. Pruning SNPs to get independent markers
```bash
plink --bfile merged_afterqc --indep-pairwise 500kb 5 0.2 --out related/merged_afterqc
```
* `--indep-pairwise 500kb 5 0.2` : perform linkage disequilibrium (LD) pruning by scanning SNPs in 500 kb windows, shifting 5 SNPs at a time, removing one SNP from pairs with r² > 0.2
* This removes correlated SNPs so relatedness estimation isn’t biased by SNPs in strong LD.
* Result: 
  * GWA-QC.prune.in (SNPs kept)
  * GWA-QC.prune.out (SNPs removed)

3. Calculate IBD (pairwise relatedness) between each pair of individuals using pruned SNPs
```bash
plink --bfile merged_afterqc --extract related/merged_afterqc.prune.in --genome --min 0.185 --out merged_afterqc
```
* `--extract GWA-QC.prune.in : use only the pruned set of SNPs for relatedness.`
* `--genome & --min 0.185 as above.`
* This reduces bias and improves accuracy of relatedness estimates.
* This command generated a .genome file that only includes pairs of individuals with π̂  > 0.185, which is the standard threshold above which individuals are considered closely related.

* I also did a plot to visualize the different types of relatedness we have :) It only shows the related individuals.

4. Remove related individuals
Remove a member from each of the pairs that are too closely related from the data set. To keep it simple you can just always remove the individual mentioned first.

```bash
plink --bfile  GWA-QC --remove wrong_ibd.txt --make-bed --out GWA-QC
```
* `--remove wrong_ibd.txt : remove individuals listed in this text file (one individual ID per line).`
* This file is typically created by you after inspecting the .genome file, picking one individual from each related pair to exclude.
* `--make-bed : create a new binary PLINK file with those individuals removed.`
* `--out GWA-QC : new cleaned dataset prefix.`

Chat CPT also wants to mention:
* Choosing which individuals to remove. The .genome file lists pairs, so you need a strategy to pick who to remove:
  * Remove the individual with lower call rate.
  * Remove duplicates or related individuals with more missingness.
  * Alternatively, remove the first individual listed as a simple rule (as your comment suggests).
* Automating removal file creation: You can write a script to parse .genome and generate wrong_ibd.txt accordingly.
* Visualize relatedness: Plot histogram of PI_HAT values to see the distribution of relatedness.
* Check for duplicates: Look for pairs with PI_HAT > 0.9 (almost identical) which can be duplicates or sample swaps.
* Consider threshold tuning: Threshold 0.185 is common but can be adjusted depending on study design.
* Document all removals: Keep track of removed individuals for transparency.

```
**Task 1: Do a QC. Are there any closely related individuals?**
After running step 3 above ("Claculating IBD") we get a .genome file as an output. It contains 7149 pairs of related individuals with relatedness >= 0.125.
`awk 'NR>1' ibdcalc_merged.genome | sort -k10,10nr | head`
Since we only have ~700 individuals, we assume, that many are related in multiple pairings. This could be the case if they are siblings or cousins, which would create many pairwise combinations.
With running `awk 'NR>1 && $10 >= 0.125 {print $1; print $2}' ibdcalc_merged.genome | sort | uniq | wc -l` we get the result, that our dataset has 123 closely related individuals. We will remove them in step 4 (above).
```

# 2. PCA [X]
```bash
plink --bfile afterqc --pca 10 --out afterqc_pca
```
This calculates Eigenvalues and Eigenvectors. I used a [Python script](https://github.com/antoniamlg/GT_GWASproject/blob/main/Scripts/PCAplot.ipynb) to visualize them graphically. <br>

**Question 2: Do a PCA plot. What does it tell you about the samples?**
* PCA plot shows, how samples are distributed and how they relate to each other based on the first two principal components (PC1 and PC2).
* both components account for 57.3% of the total variance which is a substantial portion
* samples:
  * **Clustering**: samples of same chip type cluster together -> chip has a strong influence on data structure (they are all in one place, doesn't seem to be the case here)
  * **Separation**: the colours are not well separated, so samples from different chips do not seem to have chip-wise distinct profiles in PCA space
  * **Overlap**: since the colors overlap significally, it suggests that the chip types may not be the primary source of variation.
* there are no batch effects (nice)
* the biological variation seems to be stronger than the technical variation
* the v-shaped pattern in the PCA can be indicative of population structure like different ancestries among the samples (CITE)
  * it is common that genetic differences between populations dominate in the first few principal components
  * **Admixture**: the point of the V may represent individuals with mixed ancestry, while the arms of the V represent individuals from more genetically distinct populations
  * **Gradien of ancestry**: the arms of the V may reflect a cline/gradient of ancestry, such as from Northern to Southern Europe
* How to confirm this?
  * Overlay with ancestry labels (which we don't have)
  * use reference populations (beyond scope of this project)

#### I also redid the plot with sex-colouring
You cannot see any good clustering between female and male samples (sad).

# 3. GWAS on height

## 3.1 Prepare phenotype file [X]
So that it corresponds to PLINK standards aka FID IID Phenotype.
```bash
# duplicate IID to get FID column
awk '{print $1, $1, $2}' height.txt > height_plink.txt
```
=> height = quantitative trait

## 3.2 Run the GWAS
```bash
plink --bfile your_prefix \
      --pheno height_plink.txt \
      --assoc \
      --out gwas_height
```
* what is `--assoc`:  performs basic association testing, for quantitative traits (like height) it's a linear regression

I could also use this:
```bash
plink --bfile gwa --pheno height.txt --linear --adjust --out gwas_height
```
The `--adjust´ option adds multiple testing correction results (like Bonferroni and FDR) to your output
  
## 3.3 Count significant SNPs
Count SNPs with p < 5e-8 (genome-wide significance)
```bash
awk '$9 < 5e-8' gwas_height.assoc | wc -l
```
* **gwas_height**: 7
  * basic lin. regression without any covariates (like PCs)
  * assumes that all individuals are from a single homogeneous population
  * problem: if your dataset has population structure (e.g. people from different ancestries), results are heavily infiltrated because allele frequency differences can correlate with height due to confounding - not true biological effects.
* **gwas_height_adj**: 0
  * this runs a linear regression adjusted for covariates (PCs), that controls for population structure
  * had 20.993 SNPs -> are too much though, something went wrong
* so I guess I am sticking with the 7 SNPs

## 3.4 Manhattan/QQ Plots
Find the script for the Manhattan and QQ plot under [plot script]().
-> the cool stuff
* X-axis: Chromosomes (1–22), showing the genomic position of each SNP.
* Y-axis: −log10(p-value), so higher points are more statistically significant.
* Dots: Each represents a SNP. The higher the dot, the stronger the association with the trait being studied.
* Red dot on chromosome 22: This SNP is highly significant — likely a strong candidate for further investigation.

**How many significant loci do you find?**
6
1092480     rs1349787
1092509     rs6588911
1092518    rs35178888
1116614      rs306894
1116649     rs3093474
1117147     rs1084422

=> there is no known height-association of these SNPs :(

## 3.5 QQ-plot
That will help answer:
* Are other variants also associated? → Look for a cluster or "bump" in Manhattan plot.
* General inflation? → QQ plot: early departure from expected = inflation.

Interpretation:
- assess, whether the distribution of observed p-values deviates from what would be expected under H_0 (aka no association between genetic variants and traits)
- red line: line of equality - if all your observed p-values matched the expected ones perfectly (no true associations), all points would lie on this line
- blue dots: 
  - each dot = SNP
  - dots above red line = lower p-values than expected, suggesting potential true associations
  - strong upward deviation at tail (top right) suggests significant associations (yay)
- my data:
  - most points follow red line - data is well-behaved overall
  - deviation at the end of the line likely represents a genome-wide significant hit

## 3.6 Genomic Inflation Factor
Inflation factor helps detect population structure or other confounding effects.
It’s a measure of how much your test statistics are inflated compared to what you expect under the null. λ > 1 suggests inflation, possibly from:
* population stratification (hidden population structure)
* technical artifacts
λ ≈ 1 means your p-values look well-calibrated. Calculating λ helps you check if your GWAS test statistics are unbiased or inflated. <br>
How does it work? <br>
* It adjusts your test statistics to correct for inflation by dividing by λ.
* You then recalculate p-values from the adjusted chi-squared values.
* This controls false positives if λ > 1. <br>

That answers:
* What is λ?
* What is corrected p-value of top SNP? <br>

My result => Genomic inflation factor (lambda): 0.956. This means, there is no inflation and population stratification or other confounding effects don’t look like a big problem here.

## 3.7 Genomic Control & Most Significant p After GC
My most significant GC-adjusted p-value: 5.1859116007633647e-17.

## 3.8 Replot
How to interpret differences:
* If the original and GC-adjusted Manhattan plots look almost the same (which often happens when lambda is near 1 like yours), it means:
  * Your association results are robust.
  * No inflation or confounding is affecting your findings.
* If the GC-adjusted plot had fewer SNPs passing the genome-wide significance threshold, it would mean some original signals might have been false positives due to inflation.
* Conversely, if your lambda is below 1 (like 0.956), sometimes the GC correction can make p-values a bit more significant, but usually differences are minor. <br>

**Also, on a side note, I did QC separately per chip but the GWAS not. So that could also be something to do in future analyses.**

# 4. PCA II + Adjust for PCs
There may still be population structure (ancestry differences) in your sample. People from different ancestries have different allele frequencies, and if uncorrected, that can fake associations. That’s where PCA helps:
* It captures hidden structure (e.g., European vs. African ancestry).
* You then adjust for top PCs in your GWAS to remove this confounding.

```bash
plink --bfile gwa --indep-pairwise 500kb 5 0.2 --out gwa
plink --bfile gwa --extract gwa.prune.in --pca 20 --out gwa
```

=> it looks pretty the same, only mirrored and with less variance.

## 4.1 Use PCs as covariates in GWAS
```bash
plink --bfile afterqc \
  --pheno height_adj.txt \
  --covar pca_ldpruned/pca_ldpruned.eigenvec \
  --covar-number 1-10 \
  --linear
  ```

  => and that did not work at all

# 5. Do additional analyses using Plink, GCTA, R or any other tool you might find relevant

## 5.1 Association test conditioning on the most significant SNP (Option 8)

### 5.1.1 Identify top SNP
Run `sort -g -k9 plink.assoc | head`, assuming that your association results are in `plink.assoc`.

| CHR | SNP        | BP        | A1 | TEST | NMISS | OR | STAT | P      |
|-----|------------|-----------|----|------|-------|----|------|--------|
| 10  | i6006521   | 101595996 | 1  | NA   | NA    | NA | NA   | NA     |
| 10  | i6009456   | 96521657  | 1  | NA   | NA    | NA | NA   | NA     |
| 10  | i6058695   | 88717154  | 1  | NA   | NA    | NA | NA   | NA     |
| 10  | i6059147   | 72195439  | 1  | NA   | NA    | NA | NA   | NA     |
| 10  | rs1000135  | 10854819  | 1  | NA   | NA    | NA | NA   | NA     |
| 10  | rs1003441  | 116044364 | 1  | NA   | NA    | NA | NA   | NA     |
| 10  | rs1004296  | 14815459  | 1  | NA   | NA    | NA | NA   | NA     |
| 10  | rs1006193  | 81159255  | 1  | NA   | NA    | NA | NA   | NA     |
| 10  | rs1006218  | 129149586 | 1  | NA   | NA    | NA | NA   | NA     |
| 10  | rs1006791  | 133839610 | 1  | NA   | NA    | NA | NA   | NA     |

**Why don't we take i SNPs into account?**
The i SNPs (like i6006521, i6009456, etc.) are typically imputed SNPs or non-reference variants introduced during certain types of preprocessing or genotyping. Why you should exclude them:
1. Lack of rsID: They usually don't have a known rsID, so it's harder to match them to known variants or published GWAS hits.
2. Lower confidence: They may come from lower-confidence calls or imputation with low quality.
3. Not present in external PRS tools or replication studies.

So we will use `rs1000135`.

### 5.1.2 Run conditional analysis
```bash
plink \
  --bfile data \
  --pheno phenotype.txt \               
  --covar covariates.txt \              
  --condition rs1000135 \                 
  --linear \
  --out conditioned_assoc
```
Use python script - createCovariancestxt - to get a file witht the sex info from the metadata. This will be used to control for covariance in the conditional analysis.

Comment:
* PLINK is skipping the association test because it thinks you have too many variables (covariates) compared to samples, or the phenotype is missing or invalid.
* apparently, I have extra samples in phenotype/covariate files: the phenotype and covariate files contain samples not genotyped or not included in your genotype data .fam. Maybe the pheno/covar files were created from a larger cohort or different subset.
* I am stopping, it's a lot of data wrangling and I feel it is a bit pointless

## 5.2 Distribution of phenotypes by genotype at the most significant SNP (Option 6)

## 5.2.1 Extract genotype for one SNP with PLINK
```bash
plink --bfile ../data/final_dataset/afterqc \
      --snp rs1000135 \
      --recode A --out ../data/final_dataset/rs1000135_geno
```

* the plot I get after executing the python script, shows just a flat line, so I am only seeing one genotype group (0,0). So 
  * either rs1000135_C is monomorphic (no variation), only one individual has it
  * SNP wasn't genotyped properly or filtered

**THE END**
