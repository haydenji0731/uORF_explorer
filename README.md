# uORF_explorer

**uORF_explorer** is a suite of tools used to explore potentially protein coding exons in upstream open reading frame (uORF) regions in a target genome. The suite is largely divided into two parts: (1) checking for sequence-level conservation of uORF sequences and (2) constructing protein coding exon candidates using some splicing evidence.

WIP: I'm working on polishing the readability and usability of this suite. If you have a specific request, please open an issue. I'll try to accomodate it as much as possible.

### Installation ###

We recommend that you create a Python virtual environment (venv) or a conda environment and then install the dependencies specified in requirements.txt as following:
```
$ pip install -r requirements
```
If you encounter an error in installing one of the dependencies, you can directly modify this file to better suit your environment.

### Running uORFinder ###

**uORFinder** is a Python program that maps (typically) short uORF nucleotide sequences to a target genome. We use mappy (wrapper for [minimap2](https://github.com/lh3/minimap2)) for mapping / alignment. You'll need GTF/GFF and FASTA files for both the query uORFs and target genome. You can easily extract nucleotide sequences for a given GTF/GFF file using [GffRead](https://ccb.jhu.edu/software/stringtie/gff.shtml). Run uORFinder as follows:
```
$ uORFinder.py -ogtf query_uorfs.gtf -ofa query_uorfs.fa -tgtf target_annotation.gff3 -tfa target_genome.fa -tdb target_genome_db --tmp-dir tmp_dir -o out_dir -t num_threads
```
Creating the target genome database is the biggest speed bottleneck in this program. We recommend that you save the target_genome_db file somewhere safe if you plan on reusing it. You can also increase the number of threads to expedite the alignment process. The outputs of uORFinder include conservation levels (maximum if multiple) for each query uORF in TSV format and GFF files annotating whether sequence matches occurred in the target genome.

### Running uORFconnector ###

**uORFconnector** is a Python program that constructs novel protein coding exons from uORF regions by connecting them to downstream CDS regions (see Figure below), generating what we call uORFconnected transcripts. You'll need GTF/GFF and FASTA files for the query uORFs and reference genome. The GTF/GFF file for the reference genome should contain the annotation for all transcripts and their CDS regions that you want consider connecting uORFs to. Additionally, you must provide a BED file containing the junctions to use when uORFs are spliced into a downstream CDS region. These junctions can come from any source (e.g., RNA-seq experiment).

<img width="453" alt="image" src="https://github.com/haydenji0731/uORF_explorer/assets/65577557/17ee2af6-fb92-4620-9a3c-fd5379801c33">

You can also provide a TSV file containing gene synonyms. Given a query uORF, the program searches and extracts all transcripts at the same gene locus the ORF was previously annotated on, based on gene names. Some genes have multiple aliases, and if you know that your reference annotation contains genes frequently referred to as multiple names, it's safer to provide the list of gene synonoms (see examples/gene_syns.tsv for an example). This way, you can minimize the risk of skipping over query uORFs because their corresponding gene loci weren't found in the reference annotation. uORFconnector generates a TSV file containing gene name-to-gene ID mappings, which tells you which genes were considering "missing." Run uORFconnector as follows:
```
$ uORFconnector.py -orf query_uorfs.gtf -junc junctions.bed --tmp tmp_dir \
--ref-gtf reference_annotation.gtf --ref-db reference_db -g reference_genome.fa \
-op output_prefix --collapse -syn gene_syn.tsv
```
Adding the --collapse flag ensures that your final output contains only transcripts encoding unique protein sequences. The outputs of uORFconnector include the GTF file annotating uORFconnected transcripts as well as their nucleotide and protein sequences in FASTA format, and more. uORFconnector generates GTF file without the UTR regions annotated, but if needed, you can use the add_utrs.py program located under scripts directory to add them separately.

### Running uORFconnector in de novo mode ###

This suite also includes a variation of the uORFconnector program that incorporates a deep learning-based splice site predictor called [Splam](https://github.com/Kuanhao-Chao/splam). Splam-generated scores are representative of how likely a given splice junction is real, and you can adjust the threshold to filter out low scoring splice junctions. Run uORFconnector in its de novo mode as follows:
```
$ uORFconnector_denovo.py -orf query_uorfs.gtf --ref-gtf .reference_annotation.gtf \
--ref-db  reference_db -syn gene_syn.tsv -g reference_genome.fa -op output_prefix \
-tmp tmp_dir -map gene_mappings.tsv --collapse
```
This generates pretty much identical outputs to when uORFconnector is run in its default mode, except that you'll also get the junctions.bed file outputted by Splam. 

### References ###

> Ji, H. J., & Salzberg, S. (2024). Upstream open reading frames may contain hundreds of novel human exons. bioRxiv, 2024.2003.2022.586333. [doi:10.1101/2024.03.22.586333]

> Pertea, G., & Pertea, M. (2020). GFF utilities: GffRead and GffCompare [version 2; peer review: 3 approved].
> *F1000Research*, **9**:304. [doi:10.12688/f1000research.23297.2]

> Chao, K.-H., Mao, A., Salzberg, S. L., & Pertea, M. (2023). Splam: a deep-learning-based splice site predictor that improves spliced alignments. bioRxiv, 2023.2007.2027.550754. [doi:10.1101/2023.07.27.550754]

> Danecek, P., Bonfield, J. K., Liddle, J., Marshall, J., Ohan, V., Pollard, M. O., ... & Li, H. (2021). Twelve years of SAMtools and BCFtools. *Gigascience*, > **10(2)**, giab008. [doi:10.1093/gigascience/giab008]

> Li, H. (2018). Minimap2: pairwise alignment for nucleotide sequences.
> *Bioinformatics*, **34**:3094-3100. [doi:10.1093/bioinformatics/bty191]

> Li, H. (2021). New strategies to improve minimap2 alignment accuracy.
> *Bioinformatics*, **37**:4572-4574. [doi:10.1093/bioinformatics/btab705]
