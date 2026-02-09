#!/usr/bin/env python3
"""
VCF Comparison Tool
Compares DeepSomatic output VCF against NIST ground truth VCF
Calculates precision, recall, and F1 scores
"""

import argparse
import gzip
from pathlib import Path
from collections import defaultdict
import csv


def normalize_chrom(chrom):
    """Normalize chromosome names (handle b37 vs hg19 differences)"""
    # Remove 'chr' prefix if present
    chrom = chrom.replace('chr', '')
    return chrom


def parse_vcf(vcf_path):
    """
    Parse VCF file and extract variants
    Returns: set of (chrom, pos, ref, alt) tuples
    """
    variants = set()
    vcf_path = Path(vcf_path)
    
    # Handle both compressed and uncompressed VCF
    if vcf_path.suffix == '.gz':
        opener = gzip.open
        mode = 'rt'
    else:
        opener = open
        mode = 'r'
    
    with opener(vcf_path, mode) as f:
        for line in f:
            # Skip header lines
            if line.startswith('#'):
                continue
            
            fields = line.strip().split('\t')
            if len(fields) < 5:
                continue
            
            chrom = normalize_chrom(fields[0])
            pos = int(fields[1])
            ref = fields[3]
            alt = fields[4]
            
            # Handle multiple alternate alleles
            for alt_allele in alt.split(','):
                variants.add((chrom, pos, ref, alt_allele))
    
    return variants


def compare_variants(truth_variants, query_variants):
    """
    Compare truth and query variant sets
    Returns: (TP, FP, FN) sets
    """
    # True Positives: in both truth and query
    tp = truth_variants & query_variants
    
    # False Positives: in query but not in truth
    fp = query_variants - truth_variants
    
    # False Negatives: in truth but not in query
    fn = truth_variants - query_variants
    
    return tp, fp, fn


def calculate_metrics(tp_count, fp_count, fn_count):
    """Calculate precision, recall, and F1 score"""
    # Precision: TP / (TP + FP)
    precision = tp_count / (tp_count + fp_count) if (tp_count + fp_count) > 0 else 0.0
    
    # Recall: TP / (TP + FN)
    recall = tp_count / (tp_count + fn_count) if (tp_count + fn_count) > 0 else 0.0
    
    # F1 Score: 2 * (Precision * Recall) / (Precision + Recall)
    f1 = 2 * (precision * recall) / (precision + recall) if (precision + recall) > 0 else 0.0
    
    return precision, recall, f1


def save_detailed_results(tp, fp, fn, output_path):
    """Save detailed variant comparison results to CSV"""
    with open(output_path, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['Category', 'Chromosome', 'Position', 'Reference', 'Alternate'])
        
        for variant in sorted(tp):
            writer.writerow(['True Positive', *variant])
        
        for variant in sorted(fp):
            writer.writerow(['False Positive', *variant])
        
        for variant in sorted(fn):
            writer.writerow(['False Negative', *variant])
    
    print(f"✓ Detailed results saved to: {output_path}")


def visualize_results(tp_count, fp_count, fn_count, precision, recall, f1, output_dir):
    """Generate visualization of comparison results"""
    try:
        import matplotlib.pyplot as plt
        import seaborn as sns
        
        # Set style
        sns.set_style("whitegrid")
        
        # Create figure with subplots
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))
        
        # Plot 1: Variant counts
        categories = ['True\nPositives', 'False\nPositives', 'False\nNegatives']
        counts = [tp_count, fp_count, fn_count]
        colors = ['#2ecc71', '#e74c3c', '#f39c12']
        
        bars = ax1.bar(categories, counts, color=colors, alpha=0.7, edgecolor='black')
        ax1.set_ylabel('Count', fontsize=12, fontweight='bold')
        ax1.set_title('Variant Classification', fontsize=14, fontweight='bold')
        ax1.grid(axis='y', alpha=0.3)
        
        # Add count labels on bars
        for bar, count in zip(bars, counts):
            height = bar.get_height()
            ax1.text(bar.get_x() + bar.get_width()/2., height,
                    f'{int(count)}',
                    ha='center', va='bottom', fontweight='bold', fontsize=11)
        
        # Plot 2: Performance metrics
        metrics = ['Precision', 'Recall', 'F1 Score']
        values = [precision, recall, f1]
        colors2 = ['#3498db', '#9b59b6', '#1abc9c']
        
        bars2 = ax2.bar(metrics, values, color=colors2, alpha=0.7, edgecolor='black')
        ax2.set_ylabel('Score', fontsize=12, fontweight='bold')
        ax2.set_ylim([0, 1.0])
        ax2.set_title('Performance Metrics', fontsize=14, fontweight='bold')
        ax2.grid(axis='y', alpha=0.3)
        
        # Add value labels on bars
        for bar, value in zip(bars2, values):
            height = bar.get_height()
            ax2.text(bar.get_x() + bar.get_width()/2., height,
                    f'{value:.3f}',
                    ha='center', va='bottom', fontweight='bold', fontsize=11)
        
        plt.tight_layout()
        
        # Save figure
        output_path = Path(output_dir) / "comparison_report.png"
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        print(f"✓ Visualization saved to: {output_path}")
        
    except ImportError:
        print("⚠ matplotlib/seaborn not available, skipping visualization")


def main():
    parser = argparse.ArgumentParser(description="Compare VCF files and calculate accuracy metrics")
    parser.add_argument("--truth", required=True, help="Path to ground truth VCF file")
    parser.add_argument("--query", required=True, help="Path to query VCF file (DeepSomatic output)")
    parser.add_argument("--output-dir", default="deepsomatic_output", help="Output directory for results")
    
    args = parser.parse_args()
    
    # Create output directory
    output_dir = Path(args.output_dir)
    output_dir.mkdir(exist_ok=True)
    
    print("=" * 60)
    print("VCF Comparison Tool")
    print("=" * 60)
    
    # Load VCF files
    print(f"\n📂 Loading ground truth: {args.truth}")
    truth_variants = parse_vcf(args.truth)
    print(f"   Found {len(truth_variants)} variants in ground truth")
    
    print(f"\n📂 Loading query VCF: {args.query}")
    query_variants = parse_vcf(args.query)
    print(f"   Found {len(query_variants)} variants in query")
    
    # Compare variants
    print("\n🔍 Comparing variants...")
    tp, fp, fn = compare_variants(truth_variants, query_variants)
    
    tp_count = len(tp)
    fp_count = len(fp)
    fn_count = len(fn)
    
    # Calculate metrics
    precision, recall, f1 = calculate_metrics(tp_count, fp_count, fn_count)
    
    # Print results
    print("\n" + "=" * 60)
    print("RESULTS")
    print("=" * 60)
    print(f"\n📊 Variant Counts:")
    print(f"   True Positives (TP):  {tp_count:6d}")
    print(f"   False Positives (FP): {fp_count:6d}")
    print(f"   False Negatives (FN): {fn_count:6d}")
    
    print(f"\n📈 Performance Metrics:")
    print(f"   Precision: {precision:.4f} ({precision*100:.2f}%)")
    print(f"   Recall:    {recall:.4f} ({recall*100:.2f}%)")
    print(f"   F1 Score:  {f1:.4f} ({f1*100:.2f}%)")
    
    # Save detailed results
    print(f"\n💾 Saving results...")
    csv_path = output_dir / "comparison_results.csv"
    save_detailed_results(tp, fp, fn, csv_path)
    
    # Generate visualization
    visualize_results(tp_count, fp_count, fn_count, precision, recall, f1, output_dir)
    
    print("\n" + "=" * 60)
    print("✅ Comparison complete!")
    print("=" * 60)


if __name__ == "__main__":
    main()
