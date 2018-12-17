import shlex
import unittest
import sys
import fluentdna
from unittest.mock import patch  # Python 3 only.  Used for injecting command line arguments
#import testtools  # TODO: run tests in parallel.

def fluent(cmd_string):
    """This function mimics the command line setup and uses the same argument string to invoke
    the main FluentDNA executable.  These are designed to test end to end use cases.
    NOTE: Tests need visual inspection of the results.  Not everything will throw an error.
    For example, a use case may run but not output any annotation labels."""
    s = ['fluentdna.py'] + shlex.split(cmd_string)
    with patch.object(sys, 'argv', s):
        fluentdna.main()


class FluentDNACase(unittest.TestCase):

    def test_gene_families(self):
        fluent('--layout=alignment --fasta="example_data/alignments" --outname="Test 7 Gene Families from Fraxinus"')
    def test_simple(self):
        fluent('--fasta="example_data/hg38_chr19_sample.fa" --outname="Test Simple"')
    def test_quick(self):
        fluent('"example_data/hg38_chr19_sample.fa"')
    def test_annotation_track(self):
        fluent('--fasta="example_data/gnetum_sample.fa" --outname="Test Annotation Track" --ref_annotation="example_data/Gnetum_sample_genes.gff" --annotation_width=18 --layout=annotation_track --contigs scaffold989535 scaffold103297')
    def test_annotation_highlight_and_outline(self):
        fluent('--fasta="example_data/gnetum_sample.fa" --outname="Test Annotation Highlight and Outline" --ref_annotation="example_data/Gnetum_sample_genes.gff" --query_annotation="example_data/Gnetum_query_genes.gff" --contigs scaffold989535 scaffold103297 --outname="Test Annotation Highlight and Outline"')
    def test_annotated_genome(self):
        fluent('--fasta="example_data/gnetum_sample.fa" --outname="Test Annotated Genome" --ref_annotation="example_data/Gnetum_sample_genes.gff" --contigs scaffold989535 scaffold103297')
    def test_ideogram_small(self):
        fluent('--fasta="example_data/gnetum_sample.fa" --outname="Test Ideogram Small" --contigs scaffold830595 --radix="([3,3,3,3,3,9], [5,3,3,3,3 ,53],1,1)"')
    def test_multipart_file(self):
        fluent('--fasta="example_data/Human selenoproteins.fa" --outname="Test Multipart file" --sort_contigs')
    def test_gnetum_ideogram(self):
        fluent('--fasta="example_data/gnetum_sample.fa" --outname="Test Gnetum Ideogram" --ref_annotation="example_data/Gnetum_sample_genes.gff" --query_annotation="example_data/Gnetum_query_genes.gff" --contigs scaffold830595 --radix="([3,3,3,3,3,17], [5,3,3,3,3 ,53],1,1)"')
    def test_ideogram(self):
        fluent('--fasta="example_data/hg38_chr20_sample.fa" --radix="([3,3,3,3,3, 27], [5,3,3,3,3,3,53],1,1)" --ref_annotation="example_data/Hg38_genes_chr20_sample.gff" --outname="Test Ideogram" --layout=ideogram')
    def test_custom_layout(self):
        fluent('--custom_layout="([2,3,5,7,11,13,17,999], [0,0,0,0,0,0,1,6])"  --fasta="example_data/hg38_chr19_sample.fa" --outname="Test Custom Layout"')
    def test_multiple_file_retrieval(self):
        fluent('--fasta=example_data/whole_genome_alignment/chr21_hg38_gapped.fa --extrafastas example_data/whole_genome_alignment/chr21_hg38_unique.fa example_data/whole_genome_alignment/panTro5_to_hg38_chr21_unique.fa example_data/whole_genome_alignment/panTro5_to_hg38_chr21_gapped.fa --outname="Test Multiple File Retrieval"')

    @unittest.skip("Skipped: server never closes")
    def test_server(self):
        fluent('--runserver')

    @unittest.skip("Skipped: Test is very slow")
    def test_large_chromosome(self):
        fluent('--fasta="D:\josiah\Documents\Research\Thesis - Genome Symmetry\data\Hymenoscyphus_fraxineus_EIv2.23.fa"  --outname="Test Large Chromosome" --fasta=D:\Genomes\hg38.fa --chainfile=data/hg38ToPanTro5.over.chain --extrafastas D:\Genomes\panTro5.fa --contigs chr21 --outname="Test translocations"')
    @unittest.skip("Skipped: Test is very slow")
    def test_translocation(self):
        fluent('--fasta=D:\Genomes\hg38.fa --chainfile=data/hg38ToPanTro5.over.chain --extrafastas D:\Genomes\panTro5.fa --contigs chr21 --outname="Test translocations"')


if __name__ == '__main__':
    unittest.main()
