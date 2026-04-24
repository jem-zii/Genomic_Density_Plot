# Genomic Feature Density Plot

A Python tool for computing and visualizing the spatial distribution of genomic features across chromosomes from a GFF3 annotation file.

---

## What It Does

Parses a GFF3 genome annotation file, computes feature density using a sliding window approach, and produces a multi-chromosome density plot — one panel per chromosome, color-coded, with peak density markers.

<img width="477" height="914" alt="image" src="https://github.com/user-attachments/assets/75b0088a-d97e-4e31-af2f-daf4f85f4fb3" />

---

## Output

Running the script produces a multi-panel figure where each chromosome gets its own density curve showing how features are distributed along its length. A dashed vertical line marks the position of peak density on each chromosome.

Stats are printed to the terminal:
```
Input file   : Drosophila_melanogaster.BDGP6.54.115.gff3
Feature type : gene
Window size  : 100,000 bp
Features loaded  : 13,986
Chromosomes found: ['2L', '2R', '3L', '3R', '4', 'X', 'Y']
Windows computed : 1,371
Plot saved to: Drosophila_melanogaster.BDGP6.54.115_density.png
```

---

## Setup

```bash
git clone https://github.com/YOUR_USERNAME/genomic-density-plot.git
cd genomic-density-plot
pip install -r requirements.txt
```

Download a GFF3 annotation file from [Ensembl](https://ftp.ensembl.org/pub/current_gff3/) or [NCBI](https://www.ncbi.nlm.nih.gov/genome/), then run:

```bash
python genomic_density.py --file annotation.gff3 --feature gene
```

---

## Usage

```
python genomic_density.py --file FILE [--feature FEATURE] [--window WINDOW] [--output OUTPUT]
```

| Argument | Short | Default | Description |
|----------|-------|---------|-------------|
| `--file` | `-f` | required | Path to input GFF3 annotation file |
| `--feature` | `-t` | `gene` | Feature type to analyze |
| `--window` | `-w` | `100000` | Sliding window size in base pairs |
| `--output` | `-o` | auto | Output PNG path |

### Examples

```bash
# Gene density with default 100 kbp windows
python genomic_density.py --file annotation.gff3 --feature gene

# Transposable element density with 500 kbp windows
python genomic_density.py --file annotation.gff3 --feature transposable_element --window 500000

# Save output to a specific location
python genomic_density.py --file annotation.gff3 --feature gene --output ~/Desktop/density.png
```

---

## Tested On

- *Drosophila melanogaster* BDGP6.54 (Ensembl release 115)
- Compatible with any GFF3 file following the standard specification

