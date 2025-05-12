# GWAS Pipeline Documentation

## Steps Overview
| Step         | Command                              | Description                  |
|--------------|--------------------------------------|------------------------------|
| Split IID by chip   | splitbychip.ipynb (Python)    | extract IIDs specific to one chip |
|    |                      |          |

Also check out this link [QC README](https://github.com/kaspermunch/PopulationGenomicsCourse/blob/master/Exercises/GWAS_QC/step_by_step_tutorial.md). 

# 1. Initial per-chip QC (split) [x]
- why: different chips = different SNPs, genotyping errors etc.

Python script to:
1. Read the metadata and split IIDs by chip.
2. Create `.keep` files for each chip group.

Use `.keep` to filter bed file into chip specific files:
```bash
plink --bfile gwas_data --keep /faststorage/project/populationgenomics/students/amlg/project/data/splitbychip/chipX.keep --make-bed --out /faststorage/project/populationgenomics/students/amlg/project/data/splitbychip/chipX_data
```

## 1.1 Sample-level QC []
I just checked, that's also what we did in the exercise. There, missingness and sex check was switched (but it should't matter too much right?).

### 1.1.1 Missingness []

```bash
plink --bfile GWA-QC --missing --out GWA-QC
```

### 1.1.2 Sex checks []

```bash
plink --bfile GWA-data --check-sex --out GWA-QC
```
Then, remove proplematic sex with
```bash
plink --bfile GWA-data --remove wrong_sex.txt --make-bed --out GWA-QC
```

### 1.1.3 Heterozygosity outliers []

```bash
plink --bfile GWA-QC --het --out GWA-QC 
```
calculate the observed heterozygosity rate per individual using the formula:

Het = (N(NM) − O(Hom))/N(NM) -> look at the exercises!

* then, do a plot and find outliers

Make a file with the FID and IID of all individuals that have a genotype missing rate >=0.03 or a heterozygosity rate that is more than 3 s.d. from the mean. Then use plink to remove these individuals from the data set.
```bash
plink --bfile GWA-QC --remove wrong_het_missing_values.txt --make-bed --out GWA-QC
```

### 1.1.4 Relatedness / Duplicates []

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

## 1.2 SNP-level QC []

### 1.2.1 Call rate/Missingness []
Run the --missing command again to generate the GWA-data.lmiss with the missing data rate for each SNP. <br>
Use R to make a histogram of the missing data rates (F_MISS).
Run the test-missing command and make a list of all the names of all SNPs where the differential missingness p-value is less than 1e-5. Save the list as fail-diffmiss-qc.txt.

### 1.2.2 MAF []

### 1.2.3 HWE []

after doing all this, remove all low-quality SNPs
```bash
plink --bfile GWA-QC --exclude fail-diffmiss-qc.txt --geno 0.05 --hwe 0.00001 --maf 0.01 --make-bed --out GWA-QC
```

In addition to removing SNPs identified with differential call rates between cases and controls, this command removes SNPs with call rate less than 95% with --geno option and deviation from HWE (p<1e-5) with the --hwe option. It also removes all SNPs with minor allele frequency less than a specified threshold using the --maf option.
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