#!/usr/bin/env python3
"""
DeepSomatic Variant Visualizer - Simple Version
Mostra i grafici delle varianti direttamente.
"""

import gzip
from pathlib import Path
from collections import defaultdict
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

# Configura lo stile
sns.set_style("whitegrid")
plt.rcParams['font.size'] = 10


def parse_vcf(vcf_path):
    """Legge e analizza il file VCF."""
    print(f"📖 Lettura file VCF: {vcf_path}")
    
    variants = []
    stats = defaultdict(int)
    
    open_func = gzip.open if str(vcf_path).endswith('.gz') else open
    
    with open_func(vcf_path, 'rt') as f:
        for line in f:
            if line.startswith('#'):
                continue
            
            fields = line.strip().split('\t')
            if len(fields) < 8:
                continue
            
            chrom = fields[0]
            pos = int(fields[1])
            ref = fields[3]
            alt = fields[4]
            qual = float(fields[5]) if fields[5] != '.' else 0
            
            # Determina il tipo di variante
            if len(ref) == 1 and len(alt) == 1:
                var_type = "SNP"
            elif len(ref) > len(alt):
                var_type = "DEL"
            elif len(ref) < len(alt):
                var_type = "INS"
            else:
                var_type = "COMPLEX"
            
            variant = {
                'chrom': chrom,
                'pos': pos,
                'ref': ref,
                'alt': alt,
                'qual': qual,
                'type': var_type
            }
            
            variants.append(variant)
            stats[var_type] += 1
    
    print(f"✓ Trovate {len(variants)} varianti")
    for var_type, count in stats.items():
        print(f"  • {var_type}: {count}")
    
    return variants, stats


def visualize_variants(vcf_path='deepsomatic_output/output.vcf.gz'):
    """Visualizza i grafici delle varianti."""
    vcf_path = Path(vcf_path)
    
    if not vcf_path.exists():
        print(f"✗ File VCF non trovato: {vcf_path}")
        return
    
    print("\n" + "=" * 70)
    print("🎨 DeepSomatic Variant Visualizer")
    print("=" * 70 + "\n")
    
    # Parse VCF
    variants, stats = parse_vcf(vcf_path)
    
    if not variants:
        print("\n⚠️  Nessuna variante trovata!")
        return
    
    # Crea figura con 3 subplot
    fig = plt.figure(figsize=(16, 5))
    
    # 1. Grafico a torta dei tipi di varianti
    ax1 = plt.subplot(1, 3, 1)
    colors = ['#FF6B6B', '#4ECDC4', '#45B7D1', '#FFA07A']
    explode = [0.05] * len(stats)
    
    wedges, texts, autotexts = ax1.pie(
        stats.values(),
        labels=stats.keys(),
        autopct='%1.1f%%',
        colors=colors[:len(stats)],
        explode=explode,
        shadow=True,
        startangle=90
    )
    
    for text in texts:
        text.set_fontsize(11)
        text.set_weight('bold')
    for autotext in autotexts:
        autotext.set_color('white')
        autotext.set_fontsize(9)
        autotext.set_weight('bold')
    
    ax1.set_title('Distribuzione Tipi di Varianti', fontsize=13, weight='bold', pad=15)
    
    # 2. Istogramma quality score
    ax2 = plt.subplot(1, 3, 2)
    qualities = [v['qual'] for v in variants if v['qual'] > 0]
    
    if qualities:
        ax2.hist(qualities, bins=30, color='#4ECDC4', edgecolor='black', alpha=0.7)
        ax2.axvline(np.median(qualities), color='red', linestyle='--', linewidth=2, 
                   label=f'Mediana: {np.median(qualities):.1f}')
        ax2.axvline(np.mean(qualities), color='orange', linestyle='--', linewidth=2, 
                   label=f'Media: {np.mean(qualities):.1f}')
        ax2.legend(fontsize=9)
    
    ax2.set_xlabel('Quality Score (QUAL)', fontsize=11, weight='bold')
    ax2.set_ylabel('Numero di Varianti', fontsize=11, weight='bold')
    ax2.set_title('Distribuzione Quality Score', fontsize=13, weight='bold', pad=15)
    ax2.grid(True, alpha=0.3)
    
    # 3. Scatter plot posizioni genomiche
    ax3 = plt.subplot(1, 3, 3)
    colors_map = {'SNP': '#FF6B6B', 'DEL': '#4ECDC4', 'INS': '#45B7D1', 'COMPLEX': '#FFA07A'}
    
    for v in variants:
        color = colors_map.get(v['type'], 'gray')
        ax3.scatter(v['pos'], v['qual'], c=color, alpha=0.6, s=30, edgecolors='black', linewidth=0.3)
    
    ax3.set_xlabel('Posizione Genomica (bp)', fontsize=11, weight='bold')
    ax3.set_ylabel('Quality Score', fontsize=11, weight='bold')
    ax3.set_title('Varianti lungo il Cromosoma', fontsize=13, weight='bold', pad=15)
    
    # Legenda
    from matplotlib.patches import Patch
    legend_elements = [Patch(facecolor=colors_map[t], label=f'{t} ({stats[t]})') 
                      for t in stats.keys()]
    ax3.legend(handles=legend_elements, loc='upper right', fontsize=9)
    ax3.grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    print("\n📊 Grafici generati!")
    print("💡 Chiudi la finestra per terminare.\n")
    
    plt.show()


if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description="Visualizza varianti DeepSomatic")
    parser.add_argument('--vcf', type=str, default='deepsomatic_output/output.vcf.gz',
                       help='Percorso al file VCF')
    
    args = parser.parse_args()
    
    try:
        visualize_variants(args.vcf)
    except FileNotFoundError as e:
        print(f"\n✗ Errore: {e}")
        exit(1)
    except Exception as e:
        print(f"\n✗ Errore: {e}")
        import traceback
        traceback.print_exc()
        exit(1)
