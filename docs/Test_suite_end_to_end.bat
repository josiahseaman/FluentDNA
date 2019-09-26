.\FluentDNA.exe --layout=alignment --fasta="DDV/example_data/alignments" --outname="Test 7 Gene Families from Fraxinus"
.\FluentDNA.exe --fasta="DDV/example_data/hg38_chr19_sample.fa" --outname="Test Simple"
.\FluentDNA.exe "DDV/example_data/hg38_chr19_sample.fa"
.\FluentDNA.exe --fasta="DDV/example_data/gnetum_sample.fa" --outname="Test Annotation Track" --ref_annotation="DDV/example_data/Gnetum_sample_genes.gff" --annotation_width=18 --layout=annotation_track --contigs scaffold989535 scaffold103297
.\FluentDNA.exe --fasta="DDV/example_data/gnetum_sample.fa" --outname="Test Annotation Highlight and Outline" --ref_annotation="DDV/example_data/Gnetum_sample_genes.gff" --query_annotation="DDV/example_data/Gnetum_query_genes.gff" --contigs scaffold989535 scaffold103297 --outname="Test Annotation Highlight and Outline"
.\FluentDNA.exe --fasta="DDV/example_data/gnetum_sample.fa" --outname="Test Annotated Genome" --ref_annotation="DDV/example_data/Gnetum_sample_genes.gff" --contigs scaffold989535 scaffold103297
.\FluentDNA.exe --fasta="DDV/example_data/gnetum_sample.fa" --outname="Test Ideogram Small" --contigs scaffold830595 --radix="([3,3,3,3,3,9], [5,3,3,3,3 ,53],1,1)"
.\FluentDNA.exe --fasta="DDV/example_data/Human selenoproteins.fa" --outname="Test Multipart file" --sort_contigs
.\FluentDNA.exe --fasta="DDV/example_data/gnetum_sample.fa" --outname="Test Gnetum Ideogram" --ref_annotation="DDV/example_data/Gnetum_sample_genes.gff" --query_annotation="DDV/example_data/Gnetum_query_genes.gff" --contigs scaffold830595 --radix="([3,3,3,3,3,17], [5,3,3,3,3 ,53],1,1)"
.\FluentDNA.exe --fasta="DDV/example_data/hg38_chr20_sample.fa" --radix="([3,3,3,3,3, 27], [5,3,3,3,3,3,53],1,1)" --ref_annotation="DDV/example_data/Hg38_genes_chr20_sample.gff" --outname="Test Ideogram" --layout=ideogram
.\FluentDNA.exe --runserver

