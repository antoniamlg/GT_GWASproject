# Plink

Definition: whole genome association analysis toolset

Information from 

[Input filtering - PLINK 1.9](https://www.cog-genomics.org/plink/1.9/filter)

- parses each command line as a collection of flags (—) plus parameters (which immediately follow a flag)
- part of almost every plink run:

| `--file <file>` | input data |
| --- | --- |
| `--calculations <...>` | calculations to perform on input data |
| `--out <outputpath/file>` | data output prefix, chose different for every run |

## **Citation instructions**

If you use PLINK 1.9 in any published work, please cite both the software (as an electronic resource/URL):

Package : PLINK [version]

Authors : Shaun Purcell, Christopher Chang

URL     : www.cog-genomics.org/plink/1.9/

and the manuscript(s) describing the methods you used. Our primary methods paper is:

[Chang CC, Chow CC, Tellier LCAM, Vattikuti S, Purcell SM, Lee JJ (2015) Second-generation PLINK: rising to the challenge of larger and richer datasets.](https://doi.org/10.1186/s13742-015-0047-8) GigaScience, 4.

# Standard data input

- most calculations operate on tables of samples and variant calls
- these commands are important for the project
    - `--bfile [prefix]` references the binary fileset of .bed + .bim + .fam
    - `--bed <filename>
    --bim <filename>
    --fam <filename>` replace `--bfile` and read in files separately
    - `--no-sex
    --no-fid
    --no-parents
    --no-pheno` allow you to use .fam or .ped files which lack family ID, parental ID, sex, and/or phenotype columns.

# Input filtering

exclude samples and/or variants from analysis batch based on variety of criteria

### ID lists

| --keep <filename> | space/tab delimited text file with FID’s (1st col), IID’s (2nd col) → removes all unlisted samples from current analysis |
| --- | --- |
| --remove <filename> | same as `keep` but for all listed samples |
| --keep-fam <filename>
--remove-fam <filename> | keep or remove entire family |

### Minor allele frequencies/counts

| --maf [minimum freq] | filters out all variants with minor allele frequency below provided threshold (default 0.01) |
| --- | --- |
| --max-maf <maximum freq> | upper maf bound |
| --mac <minimum count>
--max-mac <maximum count> | lower and upper minor allele count bounds |

### HWE

| --hwe <p-value> ['midp'] ['include-nonctrl’] | filters out variants which have a HWE exact test p-value below the provided threshold
→ low threshold recommended (serious genotyping erros often yield extreme p-values)
→ there is a HWE p-value calculator 😁 |
| --- | --- |

### Miscellaneous

`--allow-no-sex` samples with ambigous sex have their phenotypes set to missing → flag prevents this