import multiprocessing
from datetime import datetime
from os.path import basename, join
from glob import glob

import sys
from DNASkittleUtils.Contigs import read_contigs, write_contigs_to_file

import fluentdna
from Annotations import GFF, create_fasta_from_annotation, FeatureRep, squish_fasta


class Persist(object):
    def __init__(self, extract_contigs, lengths):
        self.output_folder = r"D:\josiah\Documents\Research\Colleagues\Yan Wong\Annotation Fastas\Original Scale"
        self.extract_contigs = extract_contigs
        self.lengths = lengths
        # self.chromosomes = ['' + str(a) for a in [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13,
        #                                           14, 15, 16, 17, 18, 19, 20, 21, 22, 'X', 'Y']]

    def do_work(self, args):
        annotation, chr_name, length = args
        annotation_fasta = join(self.output_folder, '%s_chr%s.fa' % ('hg38', chr_name))
        print(annotation_fasta)
        create_fasta_from_annotation({chr_name:annotation}, [chr_name],
                                     scaffold_lengths=[length],
                                     output_path=annotation_fasta,
                                     annotation_width=60,
                                     base_width=60,
                                     features={'gene':FeatureRep('C', 3),
                                     'exon':FeatureRep('T', 2),})
        return True

def convert_annotation():
    contigs = read_contigs(  # r"D:\josiah\Projects\DDV\DDV\data\chr20_hg38.fa")
        r"D:\Genomes\Human\hg38_chromosome.fa")
    extract_contigs = [x.name.split()[0].replace('chr', '') for x in contigs]
    lengths = [len(x.seq) for x in contigs]
    annotation = GFF(r"D:\josiah\Projects\DDV\DDV\data\Hg38_genes.gff").annotations
    dispatcher = Persist(extract_contigs, lengths)

    jobs = [(annotation[i], i, lengths[extract_contigs.index(i)]) for i in extract_contigs]

    print("Done reading everything", datetime.now() - start)
    with multiprocessing.Pool(processes=4) as pool:  # number of simultaneous processes.  Watch your RAM usage
        pool.map(dispatcher.do_work, jobs)


def squish_folder():
    input_folder = r"D:\josiah\Documents\Research\Colleagues\Yan Wong\Annotation Fastas\Original Scale"
    output_folder = r"D:\josiah\Documents\Research\Colleagues\Yan Wong\Annotation Fastas\Scale 1-60"
    for fasta in glob(join(input_folder, '*.fa')):
        contigs = read_contigs(fasta)  # there's just one
        s_contigs = squish_fasta(contigs, 1, 60)
        write_contigs_to_file(join(output_folder, basename(fasta)), s_contigs)


def set_sys_argv(args):
    sys.argv = ['fluentdna.py'] + args
    fluentdna.main()
    print('Done with ', basename(args[0][8:]))

if __name__ == '__main__':
    start = datetime.now()
    output_folder = r"D:\josiah\Documents\Research\Colleagues\Yan Wong\Annotation Fastas\Scale 1-60"
    args = []
    for fasta in glob(join(output_folder, '*.fa')):
        args.append(['--fasta=%s' % fasta,
                    '--radix=([3,7,13,3], [5,11,7,21],1,1)',
                    '--layout=ideogram', '--quick', '--no_titles', '--no_beep'])
        # set_sys_argv(args[-1])
    with multiprocessing.Pool(processes=4) as pool:  # number of simultaneous processes.  Watch your RAM usage
        pool.map(set_sys_argv, args)
    print("Done with everything", datetime.now() - start)
    # self.palette['T'] = (0, 0, 0)  # exon
    # self.palette['C'] = (0, 0, 0)  # gene