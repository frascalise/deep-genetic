import pysam
import matplotlib.pyplot as plt
from matplotlib_venn import venn2

def get_variants(vcf_path):
    """Estrae un set di tuple (cromosoma, posizione) dal VCF."""
    vcf = pysam.VariantFile(vcf_path)
    variants = set()
    for record in vcf:
        # Filtra solo le varianti che passano i filtri (se presenti)
        if 'PASS' in record.filter or len(record.filter) == 0:
            # Salviamo (CHROM, POS, REF, ALT) per essere precisi
            variants.add((record.chrom, record.pos, record.ref, record.alts[0]))
    return variants

# --- CONFIGURAZIONE ---
ground_truth_file = "deepsomatic_test/test_nist.b37_chr20_100kbp_at_10mb.vcf.gz" # Il tuo file di verità
deepsomatic_file = "deepsomatic_output/output.g.vcf.gz"       # Il tuo output

print("Caricamento varianti...")
truth_set = get_variants(ground_truth_file)
predicted_set = get_variants(deepsomatic_file)

# --- CALCOLO INTERSEZIONI ---
# True Positives (TP): Trovate da entrambi
tp = truth_set.intersection(predicted_set)
# False Positives (FP): Trovate dal modello ma non esistono nella verità
fp = predicted_set - truth_set
# False Negatives (FN): Esistono nella verità ma il modello le ha perse
fn = truth_set - predicted_set

print(f"\n--- RISULTATI ---")
print(f"Ground Truth Totali: {len(truth_set)}")
print(f"DeepSomatic Totali:  {len(predicted_set)}")
print(f"-------------------")
print(f"True Positives (Corrette): {len(tp)}")
print(f"False Positives (Rumore):  {len(fp)}")
print(f"False Negatives (Perse):   {len(fn)}")

# --- CALCOLO METRICHE ---
if len(predicted_set) > 0:
    precision = len(tp) / (len(tp) + len(fp))
    print(f"Precision: {precision:.4f} (Quanto sono affidabili le predizioni?)")
else:
    print("Precision: 0.0")

if len(truth_set) > 0:
    recall = len(tp) / (len(tp) + len(fn))
    print(f"Recall:    {recall:.4f} (Quante varianti vere ho trovato?)")
else:
    print("Recall: 0.0")

if precision + recall > 0:
    f1 = 2 * (precision * recall) / (precision + recall)
    print(f"F1-Score:  {f1:.4f}")

# --- GRAFICO DI VENN ---
plt.figure(figsize=(8, 8))
venn2([truth_set, predicted_set], set_labels=('Ground Truth', 'DeepSomatic'))
plt.title("Confronto Varianti Somatiche")
plt.show()