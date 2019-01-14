from datetime import datetime
from os.path import basename, join

from DNASkittleUtils.Contigs import read_contigs
from Annotations import GFF, create_fasta_from_annotation, FeatureRep

class Persist(object):
    def __init__(self, extract_contigs, lengths):
        self.start = datetime.now()
        self.gff_filename = r"D:\josiah\Projects\DDV\DDV\data\Hg38_genes.gff"
        self.output_folder = r"D:\josiah\Documents\Research\Colleagues\Yan Wong\Annotation Fastas\Original Scale"
        self.annotation = GFF(self.gff_filename)
        self.extract_contigs = extract_contigs
        self.lengths = lengths
        # self.chromosomes = ['' + str(a) for a in [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13,
        #                                           14, 15, 16, 17, 18, 19, 20, 21, 22, 'X', 'Y']]
        print("Done reading everything", datetime.now() - self.start)


    def do_work(self, index):
        chr_name = self.extract_contigs[index]
        annotation_fasta = join(self.output_folder, '%s_chr%s.fa' % (basename(self.gff_filename), chr_name))
        print(annotation_fasta)
        create_fasta_from_annotation(self.annotation, [chr_name],
                                     scaffold_lengths=[self.lengths[index]],
                                     output_path=annotation_fasta,
                                     annotation_width=60,
                                     base_width=60,
                                     features={'gene':FeatureRep('C', 3),
                                     'exon':FeatureRep('T', 2),})
        return True

if __name__ == '__main__':
    contigs = read_contigs(  # r"D:\josiah\Projects\DDV\DDV\data\chr20_hg38.fa")
        r"D:\Genomes\Human\hg38_chromosome.fa")
    extract_contigs = [x.name.split()[0].replace('chr', '') for x in contigs]
    lengths = [len(x.seq) for x in contigs]
    dispatcher = Persist(extract_contigs, lengths)

    jobs = list(range(len(dispatcher.extract_contigs)))

    import multiprocessing

    with multiprocessing.Pool(processes=4) as pool:  # number of simultaneous processes.  Watch your RAM usage
        pool.map(dispatcher.do_work, jobs)
    print("Done with everything", datetime.now() - dispatcher.start)
    # dispatcher.do_work(0)