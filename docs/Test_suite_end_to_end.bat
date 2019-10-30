.\FluentDNA.py --layout=alignment --fasta="FluentDNA/example_data/alignments" --outname="Test 7 Gene Families from Fraxinus"
.\FluentDNA.py --fasta="FluentDNA/example_data/hg38_chr19_sample.fa" --outname="Test Simple"
.\FluentDNA.py "FluentDNA/example_data/hg38_chr19_sample.fa"
.\FluentDNA.py --fasta="FluentDNA/example_data/gnetum_sample.fa" --outname="Test Annotation Track" --ref_annotation="FluentDNA/example_data/Gnetum_sample_genes.gff" --annotation_width=18 --layout=annotation_track --contigs scaffold989535 scaffold103297
.\FluentDNA.py --fasta="FluentDNA/example_data/gnetum_sample.fa" --outname="Test Annotation Highlight and Outline" --ref_annotation="FluentDNA/example_data/Gnetum_sample_genes.gff" --query_annotation="FluentDNA/example_data/Gnetum_query_genes.gff" --contigs scaffold989535 scaffold103297 --outname="Test Annotation Highlight and Outline"
.\FluentDNA.py --fasta="FluentDNA/example_data/gnetum_sample.fa" --outname="Test Annotated Genome" --ref_annotation="FluentDNA/example_data/Gnetum_sample_genes.gff" --contigs scaffold989535 scaffold103297
.\FluentDNA.py --fasta="FluentDNA/example_data/gnetum_sample.fa" --outname="Test Ideogram Small" --contigs scaffold830595 --radix="([3,3,3,3,3,9], [5,3,3,3,3 ,53],1,1)"
.\FluentDNA.py --fasta="FluentDNA/example_data/Human selenoproteins.fa" --outname="Test Multipart file" --sort_contigs
.\FluentDNA.py --fasta="FluentDNA/example_data/gnetum_sample.fa" --outname="Test Gnetum Ideogram" --ref_annotation="FluentDNA/example_data/Gnetum_sample_genes.gff" --query_annotation="FluentDNA/example_data/Gnetum_query_genes.gff" --contigs scaffold830595 --radix="([3,3,3,3,3,17], [5,3,3,3,3 ,53],1,1)"
.\FluentDNA.py --fasta="FluentDNA/example_data/hg38_chr20_sample.fa" --radix="([3,3,3,3,3, 27], [5,3,3,3,3,3,53],1,1)" --ref_annotation="FluentDNA/example_data/Hg38_genes_chr20_sample.gff" --outname="Test Ideogram" --layout=ideogram
.\FluentDNA.py --runserver

# Test Suite Used in Pycharm for regression testing: (may rely on some large files not included in the distribution)
fluentdna.py --layout=alignment --fasta="D:\josiah\Documents\Research\Colleagues\VCF MSA\one_alignment" --outname="Test Long MSA" --sort_contigs
fluentdna.py --fasta="example_data/gnetum_sample.fa" --outname="Test Annotated Genome" --ref_annotation="example_data/Gnetum_sample_genes.gff" --contigs scaffold989535 scaffold103297
fluentdna.py --fasta="example_data/gnetum_sample.fa" --outname="Test Annotation Highlight and Outline" --ref_annotation="example_data/Gnetum_sample_genes.gff" --query_annotation="example_data/Gnetum_query_genes.gff" --contigs scaffold989535 scaffold103297 --outname="Test Annotation Highlight and Outline"
fluentdna.py --fasta="example_data/gnetum_sample.fa" --outname="Test Annotation Track" --ref_annotation="example_data/Gnetum_sample_genes.gff" --annotation_width="18" --layout=annotation_track --contigs scaffold989535 scaffold103297
fluentdna.py --custom_layout="([2,3,5,7,11,13,17,999], [0,0,0,0,0,0,1,6])"  --fasta="example_data/hg38_chr19_sample.fa" --outname="Test Custom Layout - Prime Numbers"
fluentdna.py --layout=alignment --fasta="example_data\alignments" --outname="Test 7 Gene Families from Fraxinus" --sort_contigs
fluentdna.py --fasta="example_data/gnetum_sample.fa" --outname="Test Gnetum Ideogram" --ref_annotation="example_data/Gnetum_sample_genes.gff" --query_annotation="example_data/Gnetum_query_genes.gff" --contigs scaffold830595 --radix="([3,3,3,3,3,17], [5,3,3,3,3 ,53],1,1)" --natural_colors
fluentdna.py --fasta="D:\josiah\Projects\DDV\FluentDNA\example_data\hg38_chr20_sample.fa" --radix="([3,3,3,3,3, 27], [5,3,3,3,3,3,53],1,1)" --ref_annotation="D:\josiah\Projects\DDV\FluentDNA\example_data\Hg38_genes_chr20_sample.gff" --outname="Test Ideogram" --layout=ideogram
fluentdna.py --fasta="example_data/Human selenoproteins.fa" --outname="Test Multipart file" --sort_contigs --natural_colors
fluentdna.py --fasta=example_data/whole_genome_alignment/chr21_hg38_gapped.fa --extrafastas example_data/whole_genome_alignment/chr21_hg38_unique.fa example_data/whole_genome_alignment/panTro5_to_hg38_chr21_unique.fa example_data/whole_genome_alignment/panTro5_to_hg38_chr21_gapped.fa --outname="Test Parallel Layout"
fluentdna.py --fasta="example_data/hg38_chr19_sample.fa" --outname="Test Simple"
fluentdna.py --fasta=D:\Genomes/Human/hg38.fa --chainfile=D:\Genomes/Human/alignments/hg38ToPanTro6.over.chain --extrafastas "D:\Genomes\Chimpanzee panTro\panTro6.fa" --contigs chr21 --outname="Test translocations"
fluentdna.py "example_data\hg38_chr19_sample.fa"

