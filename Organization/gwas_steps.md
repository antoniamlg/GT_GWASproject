# GWAS Pipeline Documentation

## How do do .md files

### code blocks
Use triple backticks (```) for code blocks:

```bash plink --bfile data --mind 0.1 --geno 0.05 --maf 0.01 --make-bed --out qc_filtered ```
You can specify language after the triple backticks (e.g., bash, r, python) for syntax highlighting.

Also possible like this
```bash
plink --bfile data/raw/genotypes \
  --mind 0.1 \
  --geno 0.05 \
  --maf 0.01 \
  --make-bed \
  --out data/clean/genotypes_qc

### bullet lists
with - or 1.

### tables
| Step         | Command                              | Description                  |
|--------------|--------------------------------------|------------------------------|
| Initial QC   | plink --mind 0.1                     | Remove individuals with >10% missing data |
| MAF Filter   | plink --maf 0.01                     | Remove rare variants         |

### links
[<link name>](https://github.com/antoniamlg/GT_GWASproject)

### images
![image name](image.png)

### text formatting
- use (*) or (_) for *italic* or _italic_ text
- use (**) or (__) for __bold__ text