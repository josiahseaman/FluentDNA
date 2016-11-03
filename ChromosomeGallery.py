import os

from DDVUtils import make_output_dir_with_suffix


class ChromosomeGallery:
    def __init__(self, output_prefix):
        self.output_prefix = output_prefix
        self.output_folder = None

    # def setup_for_reference_chromosome(self, ref_chr):
    #     ending = ref_chr + '__squished' * self.squish_gaps + \
    #         '__separate_translocations' * self.separate_translocations + \
    #         '__translocations' * self.show_translocations_only + \
    #         '__aligned_only' * self.aligned_only
    #     self.output_folder = make_output_dir_with_suffix(self.output_prefix, ending)
    #     names = {'ref': ref_chr + '_%s.fa' % first_word(self.ref_source),
    #              'query': '%s_to_%s_%s.fa' % (first_word(self.query_source), first_word(self.ref_source), ref_chr)
    #              }  # for collecting all the files names in a modifiable way
    #     self.ref_sequence = pluck_contig(ref_chr, self.ref_source)  # only need the reference chromosome read, skip the others
    #
    #     return names, ref_chr


    def do_chromosome_jobs(self, chromosomes) -> list:
        """This method builds out all the assets needed to complete a job and then distributes
        them to threads which produce gapped FASTA files.
        It returns a list of Batches which contain the paths to the FASTA files that were produced."""
        self.create_parent_directory()
        # TODO: self.write_all_annotation_fastas()
        # assets = self.assets_needed_per_chromosome(chromosomes)

        assert isinstance(chromosomes, list), "'Chromosomes' must be a list! A single element list is okay."
        batches = []
        for chromosome in chromosomes:
            batches.append(self.process_chromosome(chromosome))
        return batches
        # workers = multiprocessing.Pool(6)  # number of simultaneous processes. Watch your RAM usage
        # workers.map(self._parse_chromosome_in_chain, chromosomes)


    def process_chromosome(self, chromosome):
        return NotImplementedError("This method must be overriden in Child Classes.")


    def create_parent_directory(self):
        """This uses the output_prefix to create a parent folder, then doubles up on the path so that
        dnadata/Human_  turns into dnadata/Human/Human_chr20.  This is for the purpose of grouping
        all related visualizations in a parent folder."""
        parent_folder = self.output_prefix[:-1] if self.output_prefix.endswith('_') else self.output_prefix
        parent_path = make_output_dir_with_suffix(parent_folder, '')
        self.output_prefix = os.path.join(parent_path, os.path.basename(self.output_prefix))