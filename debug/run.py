#!/usr/bin/env python3
"""
DeepSomatic Runner - Minimal Version
Esegue DeepSomatic con interfaccia semplice.
"""

import argparse
import subprocess
from pathlib import Path

def run_deepsomatic(model_type="WGS", visualize=False):
    """Esegue DeepSomatic."""
    
    # Percorsi
    base_dir = Path(__file__).parent
    test_dir = base_dir / "deepsomatic_test"
    output_dir = base_dir / "deepsomatic_output"
    output_dir.mkdir(exist_ok=True)
    
    # File di input
    bam = test_dir / "NA12878_S1.chr20.10_10p1mb.bam"
    ref = test_dir / "ucsc.hg19.chr20.unittest.fasta"
    
    # Verifica file
    if not bam.exists() or not ref.exists():
        print(f"✗ File di input mancanti in {test_dir}")
        return False
    
    print(f"🧬 DeepSomatic - {model_type}")
    print(f"📥 Input: {bam.name}")
    print(f"� Output: {output_dir}")
    print(f"🚀 Avvio analisi...\n")
    
    # Comando Docker
    base_dir_unix = str(base_dir).replace("\\", "/")
    cmd = [
        "docker", "run", "--rm",
        "-v", f"{base_dir_unix}:/data",
        "google/deepsomatic:1.8.0",
        "run_deepsomatic",
        "--model_type", model_type,
        "--ref", "/data/deepsomatic_test/ucsc.hg19.chr20.unittest.fasta",
        "--reads_tumor", "/data/deepsomatic_test/NA12878_S1.chr20.10_10p1mb.bam",
        "--output_vcf", "/data/deepsomatic_output/output.vcf.gz",
        "--output_gvcf", "/data/deepsomatic_output/output.g.vcf.gz",
        "--sample_name_tumor", "tumor_sample",
        "--num_shards", "1",
        "--logging_dir", "/data/deepsomatic_output/logs",
        "--intermediate_results_dir", "/data/deepsomatic_output/intermediate"
    ]
    
    try:
        subprocess.run(cmd, check=True)
        
        print("\n✅ Completato!")
        print(f" VCF: {output_dir / 'output.vcf.gz'}")
        
        # Visualizza se richiesto
        if visualize:
            print("\nGenerazione grafici...")
            subprocess.run(["python", "visualize_variants.py"], check=True)
        
        return True
        
    except subprocess.CalledProcessError:
        print("\n✗ Errore durante l'esecuzione")
        return False
    except KeyboardInterrupt:
        print("\n\n⚠️  Interrotto")
        return False


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Esegue DeepSomatic")
    parser.add_argument("--model", default="WGS", choices=["WGS", "WES", "PACBIO", "HYBRID"],
                       help="Tipo di modello (default: WGS)")
    parser.add_argument("--visualize", action="store_true",
                       help="Mostra grafici dopo l'analisi")
    
    args = parser.parse_args()
    
    success = run_deepsomatic(args.model, args.visualize)
    exit(0 if success else 1)
