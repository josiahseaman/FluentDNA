import os
from array import array


class GFF(object):
    def __init__(self, annotation_file):
        self.specimen, self.gff_version, \
        self.genome_version, self.date, \
        self.file_name, self.annotations \
            = self._import_gff(annotation_file)

    def _import_gff(self, annotation_file):
        assert os.path.isfile(annotation_file)

        specimen = None
        gff_version = 2
        genome_version = None
        date = None
        file_name = os.path.splitext(os.path.basename(annotation_file))[0]
        annotations = {}

        open_annotation_file = open(annotation_file, 'r')
        counter = 0
        print("Opening Annotation file...")
        for line in open_annotation_file.readlines():
            if line.startswith("#"):
                if "gff-version" in line:
                    gff_version = line.split(' ')[1]
                    if gff_version != 2:
                        raise ValueError("GFF version %s is not currently supported!" % gff_version)
                elif "genome-build" in line:
                    specimen = line.split(' ')[1]
                elif "genome-version " in line:  # NOTE: Keep the space after genome-version!!!
                    genome_version = line.split(' ')[1]
                elif "genome-date" in line:
                    date = line.split(' ')[1]
            else:
                if counter == 0:
                    print("Beginning reading of version 2 file:", annotation_file)

                counter += 1

                elements = line.split('\t')

                chromosome = elements[0]

                if chromosome not in annotations:
                    annotations[chromosome] = []
                    if len(annotations) < 10:
                        print(chromosome, end=", ")
                    elif len(annotations) == 10:
                        print('...')

                ID = counter
                source = elements[1]
                feature = elements[2]
                start = int(elements[3])
                end = int(elements[4])

                if elements[5] == '.':
                    score = None
                else:
                    score = float(elements[5])

                if elements[6] == '.':
                    strand = None
                else:
                    strand = elements[6]

                if elements[7] == '.':
                    frame = None
                else:
                    frame = int(elements[7])

                if len(elements) >= 9:
                    attribute = elements[8:]
                else:
                    attribute = None

                annotation = self.Annotation(chromosome, ID,
                                             source, feature,
                                             start, end,
                                             score, strand,
                                             frame, attribute, line)

                annotations[chromosome].append(annotation)

                # if counter % 10000 == 0:
                #     sys.stdout.write('.')

        open_annotation_file.close()

        return specimen, gff_version, genome_version, date, file_name, annotations

    class Annotation(object):
        def __init__(self, chromosome, ID, source, feature, start, end, score, strand, frame, attribute, line):
            assert isinstance(chromosome, str)
            assert isinstance(ID, int)
            assert isinstance(source, str)
            assert isinstance(feature, str)
            assert isinstance(start, int)
            assert isinstance(end, int)
            assert score is None or isinstance(score, float)
            assert isinstance(strand, str)
            assert frame is None or isinstance(frame, int)
            assert isinstance(attribute, list)

            self.chromosome = chromosome
            self.ID = ID
            self.source = source
            self.feature = feature
            self.start = start
            self.end = end
            self.score = score
            self.strand = strand
            self.frame = frame
            self.attribute = attribute
            self.line = line


def create_fasta_from_annotation(gff_filename, target_chromosome, out_name, chromosome_length):
    from DDVUtils import write_complete_fasta
    gff = GFF(gff_filename)
    filler = 'X'
    count = 0
    seq = array('u', filler * chromosome_length)
    for chromosome in gff.annotations.keys():
        if chromosome == target_chromosome or chromosome == target_chromosome.replace('chr', ''):  # only one
            for entry in gff.annotations[chromosome]:
                assert isinstance(entry, GFF.Annotation), "I'm confused"
                if entry.feature == 'exon':
                    count += 1
                    for i in range(entry.start, entry.end + 1):
                        seq[i] = 'G'
                if entry.feature == 'gene':
                    # TODO: output header JSON every time we find a gene
                    for i in range(entry.start, entry.end + 1):
                        if seq[i] == filler:
                            seq[i] = 'T'

    write_complete_fasta(out_name, seq, header='>%s\n' % target_chromosome)
    print("Done", gff.file_name)
    print("Found %i exons" % count)


def purge_annotation(gff_filename):
    features_of_interest = ['exon', 'gene']
    gff = GFF(gff_filename)
    total = 0
    kept = 0
    survivors = []
    for chromosome in gff.annotations.keys():
        for entry in gff.annotations[chromosome]:
            assert isinstance(entry, GFF.Annotation), "I'm confused"
            total += 1
            if entry.feature in features_of_interest:
                if survivors:
                    last = survivors[-1]
                    if last.start == entry.start and last.end == entry.end and entry.feature == last.feature:
                        continue  # skip this because it's a duplicate of what we already have
                kept += 1
                survivors.append(entry)

    with open(gff_filename[:-4] + '_trimmed.gtf', 'w') as out:
        for entry in survivors:
            out.write(entry.line)

    print("Done", gff.file_name)
    print("Kept %.2f percent = %i / %i" % (kept / total * 100, kept, total))


if __name__ == '__main__':
    # annotation = r'HongKong\Pan_Troglodytes_refseq2.1.4.gtf'
    # target_chromosome = 'chr20'
    # create_fasta_from_annotation(annotation, target_chromosome, 'Chimp_test_' + target_chromosome + '.fa', 63 * 1000 * 1000)

    # annotation = r'HongKong\Pan_Troglodytes_refseq2.1.4.gtf'
    annotation = r'HongKong\Homo_Sapiens_GRCH38_trimmed.gtf'
    purge_annotation(annotation)
