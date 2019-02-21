from Ideogram import Ideogram





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
            self.create_tile_layout_viz_from_fasta(args, args.fasta, args.output_name, layout)
        else:
            raise ValueError("Invalid radix settings.  Follow the example.")

    def create_tile_layout_viz_from_fasta(self, args, fasta, output_name, layout):
        print("Creating Large Image from Input Fasta...")
        layout.process_file(fasta, args.output_dir, output_name, args.no_webpage, args.contigs)

        from fluentdna import finish_webpage
        finish_webpage(args, layout, output_name)