from __future__ import print_function, division, absolute_import, \
    with_statement, generators, nested_scopes
import os
from collections import namedtuple, defaultdict
from itertools import chain

from DNASkittleUtils.Contigs import Contig

from DDV import gap_char
from DNASkittleUtils.DDVUtils import editable_str


class GFF(object):
    def __init__(self, annotation_file):
        self.specimen, self.gff_version, \
        self.genome_version, self.date, \
        self.file_name, self.annotations, self.chromosome_lengths \
            = self._import_gff(annotation_file)

    def _import_gff(self, annotation_file):
        assert os.path.isfile(annotation_file)

        specimen = None
        gff_version = '2'
        genome_version = None
        date = None
        file_name = os.path.splitext(os.path.basename(annotation_file))[0]
        annotations = {}
        chromosome_lengths = {}

        with open(annotation_file, 'r') as open_annotation_file:
            counter = 0
            print("Opening Annotation file...")
            for line in open_annotation_file.readlines():
                if line.startswith("#"):
                    if "gff-version" in line:
                        gff_version = line.split()[1]
                        if int(gff_version) != '2':
                            print("WARNING: Expecting GFF Version 2, not  %s!" % gff_version)
                    elif "genome-build" in line:
                        specimen = line.split()[1]
                    elif "genome-version " in line:  # NOTE: Keep the space after genome-version!!!
                        genome_version = line.split()[1]
                    elif "genome-date" in line:
                        date = line.split()[1]
                elif line.strip():
                    if counter == 0:
                        print("Version 2 GFF file:", annotation_file)

                    counter += 1
                    try:
                        elements = line.split('\t')

                        chromosome = elements[0]

                        if chromosome not in annotations:
                            annotations[chromosome] = []
                            chromosome_lengths[chromosome] = 0
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
                            pairs = [pair.strip() for pair in elements[8].split(';') if pair]
                            attributes = {pair.split('=')[0]: pair.split('=')[1].replace('"', '') for pair in pairs if
                                          len(pair.split('=')) == 2}
                        else:
                            attributes = {}

                        annotation = self.Annotation(chromosome, ID,
                                                     source, feature,
                                                     start, end,
                                                     score, strand,
                                                     frame, attributes, line)

                        annotations[chromosome].append(annotation)
                        chromosome_lengths[chromosome] = max(chromosome_lengths[chromosome], annotation.end)
                    except IndexError as e:
                        print(e, line)

        return specimen, gff_version, genome_version, date, file_name, annotations, chromosome_lengths

    class Annotation(object):
        def __init__(self, chromosome, ID, source, feature, start, end, score, strand, frame, attributes, line):
            assert isinstance(chromosome, str)
            assert isinstance(ID, int)
            assert isinstance(source, str)
            assert isinstance(feature, str)
            assert isinstance(start, int)
            assert isinstance(end, int)
            assert score is None or isinstance(score, float)
            assert isinstance(strand, str)
            assert frame is None or isinstance(frame, int)
            assert isinstance(attributes, dict)

            self.chromosome = chromosome
            self.ID = ID
            self.source = source
            self.feature = feature
            self.start = start
            self.end = end
            self.score = score
            self.strand = strand
            self.frame = frame
            self.attributes = attributes
            self.line = line


def handle_tail(seq_array, scaffold_lengths, sc_index):
    if scaffold_lengths is not None:
        remaining = scaffold_lengths[sc_index] - len(seq_array)
        seq_array.extend( gap_char * remaining)


def squish_fasta(scaffolds, annotation_width, base_width):
    print("Squishing annotation by %i / %i" % (base_width, annotation_width))
    squished_versions = []
    skip_size = base_width // annotation_width
    remainder = base_width - (skip_size * annotation_width)
    skips = list(chain([skip_size] * (annotation_width - 1), [skip_size + remainder]))
    for contig in scaffolds:
        work = editable_str('')
        i = 0; x = 0
        while i < len(contig.seq):
            work.append(contig.seq[i])
            i += skips[x % annotation_width]
            x += 1
        squished_versions.append(Contig(contig.name, ''.join(work)))
    return squished_versions


def create_fasta_from_annotation(gff, scaffold_names, scaffold_lengths=None, output_path=None, features=None,
                                 annotation_width=100, base_width=100):
    from DNASkittleUtils.Contigs import write_contigs_to_file, Contig
    FeatureRep = namedtuple('FeatureRep', ['symbol', 'priority'])
    if features is None:
        features = {'CDS':FeatureRep('G', 1),  # 1 priority is the most important
                    'exon':FeatureRep('T', 2),
                    'gene':FeatureRep('C', 3),
                    'mRNA':FeatureRep('A', 4),
                    'transcript':FeatureRep('N', 5)}
    symbol_priority = defaultdict(lambda: 20, {f.symbol: f.priority for f in features.values()})
    if isinstance(gff, str):
        gff = GFF(gff)  # gff parameter was a filename
    count = 0
    scaffolds = []
    for sc_index, scaff_name in enumerate(scaffold_names):  # Exact match required (case sensitive)
        if scaff_name in gff.annotations.keys():
            seq_array = editable_str(gap_char * (gff.chromosome_lengths[scaff_name] + 1))
            for entry in gff.annotations[scaff_name]:
                assert isinstance(entry, GFF.Annotation), "This isn't a proper GFF object"
                if entry.feature in features.keys():
                    count += 1
                    my = features[entry.feature]
                    for i in range(entry.start, entry.end + 1):
                        if symbol_priority[seq_array[i]] > my.priority :
                            seq_array[i] = my.symbol
                if entry.feature == 'gene':
                    # TODO: output header JSON every time we find a gene
                    pass
            handle_tail(seq_array, scaffold_lengths, sc_index)
            scaffolds.append(Contig(scaff_name, ''.join(seq_array)))
        else:
            print("No matches for '%s'" % scaff_name)
    if scaffolds:
        print("Done", gff.file_name, "Found %i features" % count, "on %i scaffolds" % len(scaffolds))
    else:
        print("WARNING: No matching scaffold names were found between the annotation and the request.")
    if annotation_width != base_width:
        scaffolds = squish_fasta(scaffolds, annotation_width, base_width)
    if output_path is not None:
        write_contigs_to_file(output_path, scaffolds)
    return scaffolds


def purge_annotation(gff_filename, features_of_interest=('exon', 'gene')):
    gff = GFF(gff_filename)
    total = 0
    kept = 0
    survivors = []
    for chromosome in gff.annotations.keys():
        for entry in gff.annotations[chromosome]:
            assert isinstance(entry, GFF.Annotation), "This isn't a GFF annotation."
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



def find_universal_prefix(annotation_list):
    """ :type annotation_list: list(GFF.Annotation) """
    names = []
    if len(annotation_list) < 2:
        return ''
    for entry in annotation_list:
        assert hasattr(entry, 'attributes'), "This isn't a proper GFF object %s" % type(entry)
        names.append(extract_gene_name(entry))  # flattening the structure
    start = 0
    for column in zip(*names):
        if all([c == column[0] for c in column]):
            start += 1
        else:
            break
    # shortened_names = [name[start:] for name in names]
    prefix = names[0][:start]
    while len(prefix) and prefix[-1].isdigit() and prefix[-1] != '0':
        prefix = prefix[:-1]  # chop off last letter
    return prefix


def extract_gene_name(entry, remove_prefix=''):
    if not entry.attributes:
        name = entry.line.split('\t')[-1]  # last part
        if '"' in name:
            name = name.split('"')[1].replace('Motif:', '')  # repeatmasker format: name inside quotes
    elif 'Name' in entry.attributes:
        name = entry.attributes['Name']
    elif 'ID' in entry.attributes:  # TODO case sensitive?
        name = entry.attributes['ID']
    elif 'gene_name'in entry.attributes:
        name = entry.attributes['gene_name']
    else:
        name = ';'.join(['%s=%s' % (key, val) for key, val in entry.attributes.items()])
    return name.replace(remove_prefix, '', 1)


if __name__ == '__main__':
    # annotation = r'DDV\data\Pan_Troglodytes_refseq2.1.4.gtf'
    # target_chromosome = 'chr20'
    # create_fasta_from_annotation(annotation, target_chromosome, 'Chimp_test_' + target_chromosome + '.fa')

    # annotation = r'DDV\data\Pan_Troglodytes_refseq2.1.4.gtf'
    annotation = r'DDV\data\Homo_Sapiens_GRCH38_trimmed.gtf'
    purge_annotation(annotation)
