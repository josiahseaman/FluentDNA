# FluentDNA Data Visualization Tool 

This application creates visualizations of FASTA formatted DNA nucleotide data.
FluentDNA generates a DeepZoomImage visualizations similar to Google Maps for FASTA files.

From 2Mbp of Genomic Sequence, FluentDNA generates this image.  Changes in nucleotide bias make
individual genome elements visible even without an annotation.  Add your annotation files to see how
they align with the sequence features.
![Example FluentDNA output of Human Chr19 2MBp](https://dnaskittle.com/ddvresults/dnadata/Test%20Simple/sources/Test%20Simple.png)

***
## FluentDNA Quick Start

You can start using FluentDNA:
 1. Downloading and unzipping the [Latest Release](https://github.com/josiahseaman/FluentDNA/releases) (Mac and Windows only).
 2. Open a terminal (command line) in the same folder you unzipped FluentDNA.
 3. Run the command `./fluentdna --fasta="DDV/example_data/Human selenoproteins.fa" --runserver`
 4. Once your private server has started, all your results available at [http://127.0.0.1:8000](http://127.0.0.1:8000).  Note that this server is not actually using the internet at all, it's just using your browser as a file viewer.

To use FluentDNA as a python module (required for Linux), follow the [pip install instructions](https://github.com/josiahseaman/FluentDNA/blob/python-master/docs/installation.md).

***

## Example Use Cases
* [Simple FASTA file (DNA)](#simple-fasta-file-dna)
* [Multi-part FASTA file (DNA)](#multi-part-fasta-file-dna)
* [Annotated Genomes](#annotated-genomes)
* [Multiple Sequence Alignment Families](#multiple-sequence-alignment-families)
* [Alignment of two Genomes](#alignment-of-two-genomes)
* [History](#history)


***

### Simple FASTA file (DNA)

Generating a basic visualization of a FASTA file downloaded from NCBI or another source is accomplished with the following commands:

**Command:** `./fluentdna --fasta="DDV/example_data/hg38_chr19_sample.fa" --outname="Test Simple"`

This generates an image pyramid with the standard legend (insert image of legend) and nucleotide number display.

**Input Data Example:**

* FA File: [DDV/example_data/hg38_chr19_sample.fa](https://github.com/josiahseaman/FluentDNA/blob/python-master/DDV/example_data/hg38_chr19_sample.fa)

**Result:** [Hg38 chr19 sample](https://dnaskittle.com/ddvresults/dnadata/Test%20Simple/)
![Example FluentDNA output of Human Chr19 2MBp](https://dnaskittle.com/ddvresults/dnadata/Test%20Simple/Test%20Simple.png)

It is also possible to generate an image file only that can be accessed with an image viewer using `--no_webpage`.

**Command:** `./fluentdna --fasta="DDV/example_data/hg38_chr19_sample.fa" --outname="Test Simple" --no_webpage`


***

### Multi-part FASTA file (DNA)

Multi-part FASTA format includes multiple sequences in the same file. A sequence record in a FASTA format consists of a single-line description (sequence name), followed by line(s) of sequence data. The first character of the description line is a greater-than (">") symbol. Multi-part format Example:
```
>seq0
AATGCCA
>seq1
GCCCTAT
```

The following command generates a multi-part FASTA file visualization:

**Input Data Example:**
* FA File: https://github.com/josiahseaman/FluentDNA/blob/python-master/DDV/example_data/Human%20selenoproteins.fa

**Command:** `./fluentdna --fasta="DDV/example_data/Human selenoproteins.fa"`
**Result:** [Human Selenoproteins](https://dnaskittle.com/ddvresults/dnadata/Human%20selenoproteins/)

![](https://dnaskittle.com/ddvresults/dnadata/Human%20selenoproteins/Human%20selenoproteins.png)

This generates a multi-scale image of the multi-part FASTA file.  Note that if you don't specify `--outname=` it will default to the name of the FASTA file.

Using this simple command, FluentDNA can visualize an entire draft genome at once.
**Result:** [Ash Tree Genome (_Fraxinus excelsior_)](https://dnaskittle.com/ddvresults/dnadata/Ash%20Tree%20Genome%20-%20BATG-0_5/)
![Fraxinus excelsior genome](https://github.com/josiahseaman/FluentDNA/raw/python-master/DDV/example_data/British%20Ash%20Tree%20Genome.png)

Additional options - see also:
- `--outname`
- `--sort_contigs`
- `--no_titles`
- `--natural_colors`
- `--base_width`

***
### Annotated Genomes
By specifying `--ref_annotation=` you can include a gene annotation to be rendered alongside your sequence.  This is currently setup to show gene introns and exons.  But the features rendered and colors used can be changed in `DDV/Annotations.py`

**Command:** `./fluentdna --fasta="DDV/example_data/gnetum_sample.fa" --ref_annotation="DDV/example_data/Gnetum_sample_genes.gff"`

**Input Data Example:**
* GFF File: https://github.com/josiahseaman/FluentDNA/blob/python-master/DDV/example_data/Gnetum_sample_genes.gff
* FA file: https://github.com/josiahseaman/FluentDNA/blob/python-master/DDV/example_data/gnetum_sample.fa


**Result:** [Gnetum montanum Annotation](https://dnaskittle.com/ddvresults/dnadata/Gnetum%20montanum%20Annotation%20-%20blue%20gene%20-%20yellow%20exon%20-%20green%20CDS/)
![Gnetum montanum Annotation](https://github.com/josiahseaman/FluentDNA/raw/python-master/DDV/example_data/Gnetum%20montanum%20Annotation%20-%20blue%20gene%20-%20yellow%20exon%20-%20green%20CDS.png)

You can download the full _Gnetum montanum_ files from [Data Dryad](https://datadryad.org//resource/doi:10.5061/dryad.0vm37).

***

### Multiple Sequence Alignment Families

To visualize a multiple sequence alignment you need to use the `--layout=alignment` option to tell FluentDNA to treat each entry in a multipart fasta file as being one row of an alignment.  To show many MSAs at once, just point `--fasta=` to a folder instead of a file.

**Input Data Example:**
* Folder with FA files: https://github.com/josiahseaman/FluentDNA/tree/python-master/DDV/example_data/alignments

**Important Note!** Make sure there are no other files in your folder besides sequence files.  If FluentDNA decides on an unreasonably long "max width" it is because it picked up a non-sequence file in the folder.

Each fasta file represents one aligned protein family.  The input fasta file should already be aligned with another program.  Each protein is represented by one entry in the fasta file with `-` inserted for gap characters like this:
```
>GeneA_in_species1
ACTCA--ACGATC------GGGT
>GeneA_in_species2
ACTCAAAACGATCTCTCTAGGGT
```

**Command:** `./fluentdna --layout=alignment --fasta="DDV/example_data\alignments" --outname="Example 7 Gene Families from Fraxinus"`

This generates a multi-scale image of the multiple alignment.  The multiple alignment results are sorted by gene name.  For a smoother layout use `--sort_contigs` which will sort them by row count (copy number).

**Result:** [Example 7 Gene Families from Fraxinus](https://dnaskittle.com/ddvresults/dnadata/Example%207%20Gene%20Families%20from%20Fraxinus/)
![](https://raw.githubusercontent.com/josiahseaman/FluentDNA/python-master/DDV/example_data/Example%207%20Gene%20Families%20from%20Fraxinus.png)

This layout allows users to check thousands of MSAs.  Here we used FluentDNA to quality check the merging software for 2,961 putative gene families: [Fraxinus Homologous Gene Groups](https://dnaskittle.com/ddvresults/dnadata/Fraxinus%20Homologous%20Gene%20Groups/)

***

### Alignment of two Genomes

This generates an alignment visualization of two genomes (A and B).  Since this is a whole genome alignment the files are very large and not included in the FluentDNA installation.  Download each of these files and unzip them into a local folder called `data/`.

#### Input Data Download Links (unzip these into data folder):
* [Human Genome FA File (hg38.fa)](http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz)
* [Alignment LiftOver Chain File (hg38ToPanTro5.over.chain)](http://hgdownload.cse.ucsc.edu/goldenPath/hg38/liftOver/hg38ToPanTro5.over.chain.gz)
* [Chimpanzee Genome FA file (panTro5.fa)](http://hgdownload.cse.ucsc.edu/goldenPath/panTro5/bigZips/panTro5.fa.gz)

**Command:** `./fluentdna --fasta=data/hg38.fa --chainfile=data/hg38ToPanTro5.over.chain --extrafastas data/panTro5.fa --chromosomes chr19 --outname="Human vs Chimpanzee"`

This generates a multi-scale image of the alignment.  There are 4 columns in this visualization:

* Column 1. Genome A (gapped entire DNA of genome A)
* Column 2. Genome A (unique DNA of A)
* Column 3. Genome B (unique DNA of B)
* Column 4. Genome B (gapped entire DNA of genome B)

**Result:** [Human vs Chimpanzee_chr19 (natural colors)](https://dnaskittle.com/ddvresults/dnadata/Parallel_hg38_and_panTro5_chr19/)
![Rows of sequence side by side Human and chimp.  Gaps where they have unique sequence.](https://github.com/josiahseaman/FluentDNA/raw/python-master/DDV/example_data/Human%20vs%20Chimpanzee_chr19.png)

Figure 1 in the paper can be seen at Nucleotide Position 548,505 which corresponds to HG38 chr19:458,731.  The difference in coordinates is due to the gaps inserted for sake of alignment.

**Note:** Whole genome alignment visualizations can be processed in batches, one visualization per chromosome.  Simply specify each of the reference chromosomes that you would like to generate.  `--outname` will be used as a prefix for the name of the folder and the visualization. For example, the above command generates a folder called "Human vs Chimpanzee_chr19".

***

## History
This project is a fork of the C# DDV developed at Concordia University.
https://github.com/photomedia/DDV/

DDV Licence:
https://github.com/photomedia/DDV/blob/master/DDV-license.txt


## Support Contact
If you run into any problems or would like to use FluentDNA in research, contact me at **josiah@newline.us**.  I'm happy to support my own software and always interested in new collaborations.