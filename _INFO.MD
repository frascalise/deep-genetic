# Deepsomatic Info

Sto usando i dati del paziente HG008.
HG008-N -> DNA normale alla nascita
HG008-T -> DNA tumorale creato al computer per vedere se l'algoritmo riesce a trovare le mutazioni inserite

## Input
Input di deepsomatic:
1. File dei dati (4 file)
    - Tumor BAM (.bam) -> file con le letture del tumore
    - Tumor Index (.bam.bai) -> indice del file BAM del tumore
    - Normal BAM (.bam) -> file con le letture del normale
    - Normal Index (.bam.bai) -> indice del file BAM del normale
2. File di riferimento (2 file)
    - Reference Genome (.fa o .fasta) -> genoma di riferimento
    - Reference Index (.fa.fai o .fasta.fai) -> indice del genoma di riferimento

Quando deepsomatic viene eseguito bisogna specificare il percorso di questi file:
- --reads_tumor -> percorso del file BAM del tumore
- --reads_normal -> percorso del file BAM del normale
- --ref -> percorso del file FASTA del genoma di riferimento
Gli indici vengono cercati automaticamente da deepsomatic, basta che si trovino nella stessa cartella dei file BAM e FASTA.

## Output
Output di deepsomatic:
1. File dei risultati VCF (.vcf o .vcf.gz) -> file con le mutazioni somatiche identificate, ogni riga rappresenta una possibile mutazione trovata nel tumore rispetto al tessuto sano
2. Indice dei risultati (.vcf.gz.tbi) -> indice del file dei risultati
3. Logs -> file di log
4. (opzionale se si usa --output_gvcf) File gVCF (.g.vcf.gz) -> contiene tutto il genoma e per ogni punto del DNA ti dice se c'e una mutazione o no

## Post-processing
Il file VCF contiene le mutazioni somatiche identificate da deepsomatic. Bisogna confrontare questo file con il file VCF del ground truth per valutare le performance di deepsomatic. Per fare questo bisogna:
1. DeepSomatic ha trovato la mutazione nel posto giusto?
    - Confronto POS e CHROM (devono essere uguali in GT e VCF)
    - Confronto REF e ALT (devono essere uguali in GT e VCF)
    - Se tutto coincide allora deepsomatic ha trovato una mutazione corretta, altrimenti sono falsi positivi o falsi negativi

2. DeepSomatic ha stimato correttamente quanta parte del tumore ha questa mutazione?
    - Confronto i VAF (Variant Allele Frequency) tra i due file VCF, se sono simili allora deepsomatic ha stimato correttamente quanta parte del tumore ha questa mutazione
    - Per trovare il VAF in VCF guardo la colonna FORMAT alla ricerca di VAF, poi prendo il valore corrispondente nella colonna TUMOR
    - Per trovare il VAF in GT guardo la colonna FORMAT alla ricerca di AF (Allele Frequency), poi prendo il valore corrispondente nella colonna HG008-T

