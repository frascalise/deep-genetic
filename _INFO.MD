# Deepsomatic Info

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

Output di deepsomatic:
1. File dei risultati (.vcf o .vcf.gz) -> file con le mutazioni somatiche identificate, ogni riga rappresenta una possibile mutazione trovata nel tumore rispetto al tessuto sano
2. Indice dei risultati (.vcf.gz.tbi) -> indice del file dei risultati
3. Logs -> file di log
4. (opzionale se si usa --output_gvcf) File gVCF (.g.vcf.gz) -> contiene tutto il genoma e per ogni punto del DNA ti dice se c'e una mutazione o no