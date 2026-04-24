"""
genomic_density.py
------------------
A tool for computing and visualizing genomic feature density
across chromosomes from a GFF3 annotation file.

1. Load and parse a GFF3 file into a pandas dataframe
2. Compute sliding window feature density per chromosome
3. Visualize as a multi-chromosome density plot
4. CLI interface via argparse

Libraries used: numpy, pandas, matplotlib, seaborn, argparse
"""


import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import seaborn as sns
from pathlib import Path


# 1. Load and parse the GFF3 file

gff_cols = [
    'seqname', 'source', 'feature', 'start', 'end', 'score',
    'strand', 'frame', 'attribute'
]

def load_gff(filepath, feature_type='gene'):
    
    """
    Read a GFF3 annotation file and return a filtered dataframe.
    """

    df = pd.read_csv(
        filepath, 
        sep='\t', 
        comment='#', 
        names=gff_cols, 
        dtype={
            'seqname': 'str',
            'start': 'Int64',
            'end': 'Int64'
        }
    )

    df['length'] = df['end'] - df['start']

    if feature_type:
        df = df[df['feature'] == feature_type]

    main_chroms = {'2L', '2R', '3L', '3R', '4', 'X', 'Y'}
    df = df[df['seqname'].isin(main_chroms)]

    return df.reset_index(drop=True)


# 2. Compute sliding window density

def compute_density(df, chrom, window_size=100_000):
    
    """
    Count how many features fall inside each sliding window
    along a single chromosome.
    """

    chrom_df = df[df['seqname'] == chrom]
    chrom_len = chrom_df['end'].max()
    windows = np.arange(0, chrom_len, window_size)

    counts = [
        len(chrom_df[
            (chrom_df['start'] >= w) &
            (chrom_df['start'] < w + window_size)
        ])
        for w in windows
    ]

    return pd.DataFrame({
        'window_start': windows,
        'count': counts,
        'chrom': chrom
    })

def compute_all_density(df, window_size=100_000):
    
    """
    Run compute_density() across every chromosome and combine results.
    """
    
    chroms = df['seqname'].unique()
 
    return pd.concat([
        compute_density(df, chrom, window_size)
        for chrom in chroms
    ]).reset_index(drop=True)


# 3. Visualize

def plot_density(density_df, feature_type, window_size, output_path):
    
    """
    Plot a multi-chromosome genomic feature density figure.
    Each chromosome gets its own panel.
    """

    chroms = density_df['chrom'].unique()
    n = len(chroms)
    palette = sns.color_palette('husl', n)

    fig = plt.figure(figsize=(12, 3 * n), facecolor='#fafafa')
    fig.suptitle(
        f'Genomic Feature Density - {feature_type} | window = {window_size:,} bp',
        fontsize=14, fontweight='bold', y=1.01
    )

    gs = gridspec.GridSpec(n, 1, figure=fig, hspace=0.6)

    for i, chrom in enumerate(chroms):
        ax = fig.add_subplot(gs[i])
        chrom_data = density_df[density_df['chrom'] == chrom]

        x = chrom_data['window_start'] / 1_000_000   # convert bp to Mbp
        y = chrom_data['count']

        ax.fill_between(x, y, alpha=0.4, color=palette[i])
        ax.plot(x, y, linewidth=1.0, color=palette[i])

        peak_x = x.iloc[y.values.argmax()]
        ax.axvline(peak_x, color=palette[i], linestyle='--',
                   linewidth=0.8, alpha=0.7, label=f'Peak: {peak_x:.1f} Mbp')
        
        ax.set_title(chrom, fontsize=11, fontweight='bold', loc='left')
        ax.set_ylabel('Features / window', fontsize=9)
        ax.set_xlabel('Position (Mbp)', fontsize=9)
        ax.legend(fontsize=8, loc='upper right')
        ax.grid(True, alpha=0.2)
        ax.set_xlim(0, x.max())

    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight', facecolor='#fafafa')
    print(f"Plot saved to: {output_path}")
    plt.show()


# 4. CLI interface via argparce

def parse_args():

    """
    Define and parse command-line arguments so the script can be
    run from the terminal with options like --file, --feature, --window.
    """

    parser = argparse.ArgumentParser(
        description='Compute and visualize genomic feature density from a GFF3 file.'
    )

    parser.add_argument(
        '--file', '-f',
        required=True,
        help='Path to input GFF3 annotation file'
    )
    parser.add_argument(
        '--feature', '-t',
        default='gene',
        help='Feature type to analyze (default: gene)'
    )
    parser.add_argument(
        '--window', '-w',
        type=int,
        default=100_000,
        help='Sliding window size in base pairs (default: 100000)'
    )
    parser.add_argument(
        '--output', '-o',
        default=None,
        help='Output PNG path (default: saves next to input file)'
    )

    return parser.parse_args()


# Main function
def main():
    args = parse_args()
    gff_path = Path(args.file)

    output_path = (
        Path(args.output) if args.output 
        else gff_path.parent / f'{gff_path.stem}_density.png'
    )

    print(f"Input file: {gff_path}")
    print(f"Feature type: {args.feature}")
    print(f"Window size: {args.window:,} bp")
    print(f"Output path: {output_path}")

    print("\nLoading annotation file...")
    df = load_gff(gff_path, feature_type=args.feature)
    print(f"Features loaded  : {len(df):,}")
    print(f"Chromosomes found: {list(df['seqname'].unique())}")

    print("\nComputing density...")
    density_df = compute_all_density(df, window_size=args.window)
    print(f"Windows computed : {len(density_df):,}")

    print("\nGenerating plot...")
    plot_density(density_df, args.feature, args.window, output_path)
 
 
if __name__ == "__main__":
    main()
