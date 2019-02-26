from glob import glob
from os.path import join, basename, isdir

from Ideogram import Ideogram


def batch_render(args, folder, output_name, renderer):
    renderer.render_multiple_files(folder, args.output_dir, output_name)
    #Initial setup
    #Point_mapping grows as needed, just copy from N-1 each time
    #preserve same image canvas the whole time
    # longest = sum(len(c.seq) for c in self.contigs)

    # longest = 2 * 1000 * 1000
    # renderer.each_layout[0].build_coordinate_mapping(longest)
    # accumulated_point_mapping = self.each_layout[0].point_mapping
    # for file_index, fasta in enumerate(glob(join(folder, '*.fa*'))):
    #     print("Processing", basename(fasta))
    #
    #     # Calculate longest layout first, then reuse
    #     # Layout contigs one at a time
    #     for layout_index, contig in enumerate(self.contigs):
    #         self.i_layout = layout_index  # for intercompatibility
    #         renderer = self.each_layout[self.i_layout]
    #         if file_index:  # not the first
    #             renderer.point_mapping = accumulated_point_mapping
    #
    #
    #     renderer.process_file(fasta, args.output_dir, output_name, args.no_webpage, [])

    from fluentdna import finish_webpage
    finish_webpage(args, renderer, output_name)


class IdeogramManager:
    def __init__(self, args):
        assert args.radix, "You must provide a --radix argument for Ideograms."
        radix_settings = eval(args.radix)
        if len(radix_settings) == 4 and \
            type(radix_settings[0]) == type(radix_settings[1]) == type([]) and \
            type(radix_settings[2]) == type(radix_settings[3]) == type(1):
            layout = Ideogram(radix_settings,
                              ref_annotation=args.ref_annotation, query_annotation=args.query_annotation,
                              repeat_annotation=args.repeat_annotation,
                              low_contrast=args.low_contrast, use_titles=args.use_titles,
                              use_labels=args.use_labels)
            if isdir(args.fasta):
                batch_render(args, args.fasta, args.output_name, layout)
            else:
                from fluentdna import create_tile_layout_viz_from_fasta
                create_tile_layout_viz_from_fasta(args, args.fasta, args.output_name, layout)
        else:
            raise ValueError("Invalid radix settings.  Follow the example.")

