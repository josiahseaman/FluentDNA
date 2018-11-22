#!/usr/bin/env python
"""
FluentDNA DDV 2.0 is a new version of DDV written in Python that allows you to generate a single image
for an entire genome.  It was necessary to switch platforms and languages because of intrinsic
limitations in the size of image that could be handled by: C#, DirectX, Win2D, GDI+, WIC, SharpDX,
or Direct2D. We tried a lot of options.

The python version has matured significantly past the previous feature set.

"""
from __future__ import print_function, division, absolute_import, \
    with_statement, generators, nested_scopes

import os
import sys
#############################################################################
# IMPORTANT!  Make sure there are import here for non-builtin packages.  Those go below.
#############################################################################
# print("Setting up Python...")

if getattr(sys, 'frozen', False):
    BASE_DIR = os.path.dirname(sys.executable)
    os.environ["PATH"] += os.pathsep + os.path.join(BASE_DIR, 'bin')
    os.environ["PATH"] += os.pathsep + os.path.join(BASE_DIR, 'bin', 'env')
else:
    try:
        BASE_DIR = os.path.dirname(os.path.abspath(__file__))
    except:  # just in case __file__ isn't defined in some contexts
        import DDV
        BASE_DIR = os.path.dirname(DDV.__file__)
print('Running in:', BASE_DIR)

sys.path.append(BASE_DIR)
sys.path.append(os.path.join(BASE_DIR, 'bin'))
sys.path.append(os.path.join(BASE_DIR, 'bin', 'env'))

os.chdir(BASE_DIR)

import multiprocessing
multiprocessing.freeze_support()

# ----------BEGIN MAIN PROGRAM----------
from DDV import VERSION

import argparse

from DNASkittleUtils.CommandLineUtils import just_the_name
from DDV.DDVUtils import create_deepzoom_stack, make_output_dir_with_suffix, base_directories, \
    hold_console_for_windows, beep
from DDV.ParallelGenomeLayout import ParallelLayout
from DDV.AnnotatedTrackLayout import  AnnotatedTrackLayout
from DDV.Ideogram import Ideogram
from DDV.HighlightedAnnotation import HighlightedAnnotation
from DDV.ChainParser import ChainParser
from DDV.UniqueOnlyChainParser import UniqueOnlyChainParser
from DDV.AnnotatedAlignment import AnnotatedAlignment
from DDV.TileLayout import TileLayout
from DDV.TransposonLayout import TransposonLayout
from DDV.MultipleAlignmentLayout import MultipleAlignmentLayout
from DNASkittleUtils.Contigs import write_contigs_to_file, read_contigs

if sys.platform == 'win32':
    OS_DIR = 'windows'
    EXTENSION = '.exe'
    SCRIPT = '.cmd'
else:
    OS_DIR = 'linux'
    EXTENSION = ''
    SCRIPT = ''


def query_yes_no(question, default='yes'):
    valid = {'yes': True, 'y': True, "no": False, 'n': False}

    if default is None:
        prompt = " [y/n] "
    elif default in ['yes', 'y']:
        prompt = " [Y/n] "
    elif default in ['no', 'n']:
        prompt = " [y/N] "
    else:
        raise ValueError("Invalid default answer!")

    while True:
        sys.stdout.write('\n' + question + prompt)

        choice = input().lower()

        if default is not None and choice == '':
            return valid[default]
        elif choice in valid:
            return valid[choice]
        else:
            sys.stdout.write("Please respond with 'yes' or 'no'.\n")


def run_server(home_directory):
    try:
        from http import server
        from socketserver import TCPServer
    except ImportError:  # Python 2 imports
        import SimpleHTTPServer as server
        from SocketServer import TCPServer

    print("Setting up HTTP Server based from", home_directory)
    os.chdir(home_directory)

    ADDRESS = "127.0.0.1"
    PORT = 8000

    handler = server.SimpleHTTPRequestHandler
    httpd = TCPServer((ADDRESS, PORT), handler)

    print("Open a browser at http://%s:%s" %(ADDRESS, str(PORT)))
    httpd.serve_forever()


def done(args, output_dir):
    """Ensure that server always starts when requested.
    Otherwise system exit."""
    if args.run_server:
        run_server(output_dir)
    beep()
    hold_console_for_windows()
    sys.exit(0)


def ddv(args):
    SERVER_HOME, base_path = base_directories(args.output_name)

    if not args.layout and args.run_server:
        done(args, SERVER_HOME)



    if args.layout == "NONE":  # Complete webpage generation from existing image
        layout = TileLayout(use_titles=args.use_titles, sort_contigs=args.sort_contigs,
                            low_contrast=args.low_contrast, base_width=args.base_width,
                            custom_layout=args.custom_layout)
        layout.generate_html(args.output_dir, args.output_name)
        print("Creating Deep Zoom Structure for Existing Image...")
        create_deepzoom_stack(args.image, os.path.join(args.output_dir, 'GeneratedImages', "dzc_output.xml"))
        print("Done creating Deep Zoom Structure.")
        done(args, args.output_dir)

    elif args.layout == "tiled":  # Typical Use Case
        # TODO: allow batch of tiling layout by chromosome
        create_tile_layout_viz_from_fasta(args, args.fasta, args.output_name)
        done(args, args.output_dir)

    # ==========TODO: separate views that support batches of contigs============= #
    elif args.layout == 'transposon':
        layout = TransposonLayout()
        # if len(args.contigs) != 1:
        #     raise NotImplementedError("Chromosome Argument requires exactly one chromosome e.g. '--contigs chr12'")
        layout.process_all_repeats(args.fasta, args.output_dir, just_the_name(args.output_dir), args.ref_annotation, args.contigs)
        print("Done with Transposons")
        done(args, args.output_dir)

    elif args.layout == 'alignment':
        layout = MultipleAlignmentLayout(sort_contigs=args.sort_contigs)
        layout.process_all_alignments(args.fasta,
                                      args.output_dir,
                                      args.output_name)
        finish_webpage(args, layout, args.output_name)
        print("Done with Alignments")
        done(args, args.output_dir)

    elif args.layout == "parallel":  # Parallel genome column layout OR quad comparison columns
        if not args.chain_file:  # life is simple
            # TODO: support drag and drop of multiple files
            create_parallel_viz_from_fastas(args, len(args.extra_fastas) + 1, args.output_dir,
                                            args.output_name, [args.fasta] + args.extra_fastas)
            done(args, args.output_dir)
        else:  # parse chain files, possibly in batch
            chain_parser = ChainParser(chain_name=args.chain_file,
                                       first_source=args.fasta,
                                       second_source=args.extra_fastas[0],
                                       output_prefix=base_path,
                                       trial_run=args.trial_run,
                                       separate_translocations=args.separate_translocations,
                                       no_titles=args.no_titles,
                                       squish_gaps=args.squish_gaps,
                                       show_translocations_only=args.show_translocations_only,
                                       aligned_only=args.aligned_only,
                                       extract_contigs=args.contigs)
            print("Creating Gapped and Unique Fastas from Chain File...")
            batches = chain_parser.parse_chain(args.contigs)
            del chain_parser
            print("Done creating Gapped and Unique.")
            for batch in batches:  # multiple contigs, multiple views
                create_parallel_viz_from_fastas(args, len(batch.fastas),
                                                batch.output_folder, batch.output_folder,
                                                batch.fastas)
            done(args, SERVER_HOME)
    elif args.layout == "annotation_track":
        layout = AnnotatedTrackLayout(args.fasta, args.ref_annotation, args.annotation_width)
        layout.render_genome(args.output_dir, args.output_name, args.contigs)
        finish_webpage(args, layout, args.output_name)
        done(args, args.output_dir)
    elif args.layout == "annotated":
        layout = HighlightedAnnotation(args.ref_annotation, args.query_annotation, args.repeat_annotation,
                                       use_titles=args.use_titles, sort_contigs=args.sort_contigs,
                                       low_contrast=args.low_contrast, base_width=args.base_width,
                                       custom_layout=args.custom_layout)
        layout.process_file(args.fasta, args.output_dir, args.output_name,
                            args.no_webpage, args.contigs)
        finish_webpage(args, layout, args.output_name)
        done(args, args.output_dir)

    elif args.layout == "unique":
        """UniqueOnlyChainParser(chain_name='data\\hg38ToPanTro4.over.chain',
                               first_source='data\\hg38.fa',
                               second_source='',
                               output_folder_prefix='Hg38_unique_vs_panTro4_')"""
        unique_chain_parser = UniqueOnlyChainParser(chain_name=args.chain_file,
                                                    first_source=args.fasta,
                                                    second_source=args.fasta,
                                                    output_prefix=base_path,
                                                    trial_run=args.trial_run,
                                                    separate_translocations=args.separate_translocations)
        batches = unique_chain_parser.parse_chain(args.contigs)
        print("Done creating Gapped and Unique Fastas.")
        del unique_chain_parser
        combine_files(batches, args, args.output_name)
        # for batch in batches:
        #     create_tile_layout_viz_from_fasta(args, batch.fastas[0], batch.output_folder, args.output_name)
        done(args, SERVER_HOME)

    elif args.layout == 'ideogram':
        assert args.radix, "You must provide a --radix argument for Ideograms."
        radix_settings = eval(args.radix)
        if len(radix_settings) == 4 and \
            type(radix_settings[0]) == type(radix_settings[1]) == type([]) and \
            type(radix_settings[2]) == type(radix_settings[3]) == type(1):
            layout = Ideogram(radix_settings,
                              ref_annotation=args.ref_annotation, query_annotation=args.query_annotation,
                              repeat_annotation=args.repeat_annotation,
                              low_contrast=args.low_contrast, use_titles=args.use_titles)
            create_tile_layout_viz_from_fasta(args, args.fasta, args.output_name, layout)
        else:
            print("Invalid radix settings.  Follow the example.")
        done(args, args.output_dir)


    elif args.ref_annotation and args.layout != 'transposon':  # parse chain files, possibly in batch
        anno_align = AnnotatedAlignment(chain_name=args.chain_file,
                                        first_source=args.fasta,
                                        first_annotation=args.ref_annotation,
                                        second_source=args.extra_fastas[0],
                                        second_annotation=args.query_annotation,
                                        output_prefix=base_path,
                                        trial_run=args.trial_run,
                                        separate_translocations=args.separate_translocations,
                                        squish_gaps=args.squish_gaps,
                                        show_translocations_only=args.show_translocations_only,
                                        aligned_only=args.aligned_only)
        print("Creating Aligned Annotations using Chain File...")
        batches = anno_align.parse_chain(args.contigs)
        del anno_align
        print("Done creating Gapped Annotations.")
        for batch in batches:  # multiple contigs, multiple views
            create_parallel_viz_from_fastas(args, len(batch.fastas), args.output_dir, args.output_name,
                                            batch.fastas)
        done(args, SERVER_HOME)
    else:
        raise NotImplementedError("What you are trying to do is not currently implemented!")


def create_parallel_viz_from_fastas(args, n_genomes, output_dir, output_name, fastas):
    print("Creating Large Comparison Image from Input Fastas...")
    layout = ParallelLayout(n_genomes=n_genomes, low_contrast=args.low_contrast, base_width=args.base_width)
    layout.process_file(output_dir, output_name, fastas, args.no_webpage, args.contigs)

    finish_webpage(args, layout, output_name)



def create_tile_layout_viz_from_fasta(args, fasta, output_name, layout=None):
    print("Creating Large Image from Input Fasta...")
    if layout is None:
        layout = TileLayout(use_titles=args.use_titles, sort_contigs=args.sort_contigs,
                            low_contrast=args.low_contrast, base_width=args.base_width,
                            custom_layout=args.custom_layout)
    layout.process_file(fasta, args.output_dir, output_name, args.no_webpage, args.contigs)

    finish_webpage(args, layout, output_name)


def combine_files(batches, args, output_name):
    from itertools import chain
    contigs = list(chain(*[read_contigs(batch.fastas[0]) for batch in batches]))
    fasta_output = output_name + '.fa'
    write_contigs_to_file(fasta_output, contigs)
    create_tile_layout_viz_from_fasta(args, fasta_output, output_name)


def finish_webpage(args, layout, output_name):
    final_location = layout.final_output_location
    print("Done creating Large Image at ", final_location)
    if not args.no_webpage:
        layout.generate_html(args.output_dir, output_name)
        del layout
        print("Creating Deep Zoom Structure from Generated Image...")
        create_deepzoom_stack(os.path.join(args.output_dir, final_location),
                              os.path.join(args.output_dir, 'GeneratedImages', "dzc_output.xml"))
        print("Done creating Deep Zoom Structure.")
    else:
        del layout


def main():
    if len(sys.argv) == 2 and not sys.argv[1].startswith('-'):  # there's only one input and it does have a flag
        print("--Starting in Quick Mode--")
        print("This will convert the one FASTA file directly to an image and place it in the same "
              "folder as the image for easy access.  "
              # "The scaffolds will be sorted by length for best layout."
              "Recommend you open large files with 'Windows Photo Viewer'.")
        sys.argv[1] = '--fasta=' + sys.argv[1]
        sys.argv.append("--quick")
    if "--quick" in sys.argv:
        sys.argv.append("--no_webpage")  # don't generate a full webpage (deepzoom is time consuming)

        # sys.argv.append("--sort_contigs")

    parser = argparse.ArgumentParser(usage="%(prog)s [options]",
                                     description="Creates visualizations of FASTA formatted DNA nucleotide data.",
                                     add_help=True)

    parser = argparse.ArgumentParser(prog='fluentdna')
    parser.add_argument('--quick',
                        action='store_true',
                        help="Shortcut for dropping the file on fluentdna.exe.  Only an image will be generated "
                             "in the same directory as the FASTA.  This is the default behavior if you drop "
                             "a file onto the program or a filepath is the only argument.",
                        dest="quick")

    parser.add_argument("-f", "--fasta",
                        type=str,
                        help="Path to main FASTA file to process into new visualization.",
                        dest="fasta")
    parser.add_argument("-o", "--outname",
                        type=str,
                        help="What to name the output folder (not a path). Defaults to name of the fasta file.",
                        dest="output_name")
    parser.add_argument("-r", "--runserver",
                        action='store_true',
                        help="Run Web Server after computing.",
                        dest="run_server")
    parser.add_argument("-c", "--contigs",
                        nargs='+',
                        type=str,
                        help="List contigs you'd like visualized from the file separated by spaces. "
                             "This can be used to pluck out your contig of interest from a large file. "
                             "REQUIRED for Chain File alignments.",
                        dest="contigs")
    parser.add_argument("-cc", "--chromosomes", nargs='+', type=str,
                        help="Synonym for --contigs for backwards compatibility.", dest="contigs")

    parser.add_argument('-s', '--sort_contigs',
                        action='store_true',
                        help="Sort the entries of the fasta file by length.  This option will kick in "
                             "automatically if your file has more than 10,000 separate FASTA entries.",
                        dest="sort_contigs")
    parser.add_argument('-nc', '--natural_colors',
                        action='store_true',
                        help="Use low contrast, natural colors that are easier on the eyes",
                        dest="low_contrast")
    parser.add_argument("-l", "--layout",
                        type=str,
                        help="The type of layout to perform. Will autodetect between Tiled and "
                            "Parallel. Only needed if you want non-default option like 'alignment', "
                             "'unique' or 'annotation_track'.",
                        choices=["tiled", "annotated", "ideogram", "alignment", "annotation_track",
                                 "parallel", "unique", "transposon"],
                        dest="layout")  # Don't set a default so we can do error checking on it later
    parser.add_argument("-x", "--extrafastas",
                        nargs='+',
                        type=str,
                        help="Path to secondary FASTA files to process when doing Parallel layout.",
                        dest="extra_fastas")
    parser.add_argument("-bw", "--base_width",
                        help="Overrides the default 100bp column width in standard --layout=tiled. "
                        "Use this only if you are trying to accomplish something custom. "
                        "The rest of the layout will ratio adjust, so base_width=200 will produce "
                        "columns that are 2,000 lines tall and rows containing 40 Mbp, etc.",
                        dest="base_width")

    parser.add_argument("-nt", "--no_titles",
                        action='store_true',
                        help="No gaps for a title.  Useful when combined with separate_translocations",
                        dest="no_titles")
    parser.add_argument("-nw", "--no_webpage",
                        action='store_true',
                        help="Use if you only want an image.  No webpage or zoomstack will be calculated.  "
                        "You can use --image option later to resume the process to get a deepzoom stack.",
                        dest="no_webpage")
    parser.add_argument("-q", "--trial_run",
                        action='store_true',
                        help="Only show the first 1 Mbp.  This is a fast run for testing.",
                        dest="trial_run")
    ### Chain Files
    parser.add_argument("-cf", "--chainfile",
                        type=str,
                        help="Path to Chain File when doing Parallel Comparisons layout.",
                        dest="chain_file")
    parser.add_argument("-t", "--separate_translocations",
                        action='store_true',
                        help="Don't edit in translocations, list them at the end.",
                        dest="separate_translocations")
    parser.add_argument("-g", "--squish_gaps",
                        action='store_true',
                        help="If two gaps are approximately the same size, subtract the intersection.",
                        dest="squish_gaps")
    parser.add_argument("-k", "--show_translocations_only",
                        action='store_true',
                        help="Used to highlight the locations of translocations (temporary)",
                        dest='show_translocations_only')
    parser.add_argument("-a", "--aligned_only",
                        action='store_true',
                        help="Don't show the unaligned pieces of ref or query sequences.",
                        dest='aligned_only')

    ### Annotations
    parser.add_argument("-ra", "--ref_annotation",
                        type=str,
                        help="Path to Annotation File for Reference Genome (first).",
                        dest="ref_annotation")
    parser.add_argument("-qa", "--query_annotation",
                        type=str,
                        help="Path to Annotation File for Query Genome (second).",
                        dest="query_annotation")
    parser.add_argument("-rp", "--repeat_annotation",
                        type=str,
                        help="Path to Annotation File for Repeats which will be shaded.",
                        dest="repeat_annotation")

    parser.add_argument("-aw", "--annotation_width",
                        help="Overrides the default 100 pixel column width for annotations. "
                        "annotation_width=1 will only sample one pixel per display line, "
                        "skipping intermediate intervals.  If annotated features are less than"
                        "base_width / annotation_width bp in length it's possible they won't be visible.",
                        dest="annotation_width")

    ### Other
    parser.add_argument("-i", "--image",
                        type=str,
                        help="Path to already computed big image to process with DeepZoom. "
                             "No layout will be performed if an image is passed in.",
                        dest="image")
    parser.add_argument("-rx", "--radix",
                        type=str,
                        help="String that is a python literal for the radix settings. "
                             "x and y radices, and scale\n"
                             "Example: '([5,5,5,5,11], [5,5,5,5,5 ,53], 1, 1)'",
                        dest="radix")
    parser.add_argument("-cl", "--custom_layout",
                        type=str,
                        help='Changes the layout based on ([number of repeating units], [padding between units])'
                             'Custom layout must be formatted as two integer lists of euqal length.\n'
                             'For example: --custom_layout="([10,100,100,10,3,999], [0,0,0,3,18,108])"',
                        dest="custom_layout")
    parser.add_argument('-n', '--update_name', dest='update_name', help='Query for the name of this program as known to the update server', action='store_true')
    parser.add_argument('-v', '--version', dest='version', help='Get current version of program.', action='store_true')

    args = parser.parse_args()

    # Respond to an updater query
    if args.update_name:
        print("DDV")
        sys.exit(0)
    elif args.version:
        print(VERSION)
        sys.exit(0)

    # Errors
    if args.base_width:
        args.base_width = int(args.base_width)
    if args.annotation_width:
        args.annotation_width = int(args.annotation_width)

    #Layout Defaults
    if not args.layout:
        if args.extra_fastas:  # separate because unique can use a chain file without extra_fastas
            args.layout = 'parallel'
        if args.fasta:
            if args.radix:
                args.layout = 'ideogram'
            elif args.chain_file:
                if not args.extra_fastas:
                    args.layout = 'unique'
            elif args.ref_annotation or args.query_annotation or args.repeat_annotation:
                args.layout = "annotated"
            else:
                args.layout = 'tiled'
    if args.image and not args.layout:
        args.layout = "NONE"


    if args.image and (args.fasta or args.layout or args.extra_fastas or args.chain_file):
        parser.error("No layout will be performed if an existing image is passed in! "
                     "Please only define an existing 'image' and the desired 'outfile'.")
    if not args.image and not args.fasta and not args.run_server:
        parser.error('Please define a a file to process.  Ex: fluentdna.py --fasta="example_data/phiX.fa"')
    if args.image and args.no_webpage:
        parser.error("This parameter combination doesn't make sense.  You've provided a precalculated image"
                     "and asked DDV to only generate an image with no DeepZoom stack or webpage.")

    if args.extra_fastas and not args.layout:
        args.layout = "parallel"
    if args.layout and args.layout == "parallel" and not args.extra_fastas:
        parser.error("When doing a Parallel, you must at least define 'extrafastas'!")
    # if args.layout and args.layout == 'unique' and args.extra_fastas:
    #     parser.error("For Unique view, you don't need to specify 'extrafastas'.")
    # if args.contigs and not (args.chain_file or args.layout == 'transposon'):
    #     parser.error("Listing 'contigs' is only relevant when parsing Chain Files or Repeats!")
    # if args.extra_fastas and "parallel" not in args.layout:
    #     parser.error("The 'extrafastas' argument is only used when doing a Parallel layout!")
    if args.chain_file and args.layout not in ["parallel", "unique"]:
        parser.error("The 'chainfile' argument is only used when doing a Parallel or Unique layout!")
    if args.chain_file and args.extra_fastas and len(args.extra_fastas) > 1:
        parser.error("Chaining more than two samples is currently not supported! Please only specify one 'extrafastas' when using a Chain input.")
    if args.layout == "unique" and not args.chain_file:
        parser.error("You must have a 'chainfile' to make a Unique layout!")
    if args.show_translocations_only and args.separate_translocations:
        parser.error("It just doesn't make sense to ask to show translocations in context while separating them.  You've got to pick one or the other.")

    # Set post error checking defaults
    if not args.contigs and args.chain_file and args.layout != 'unique':
        print("Error: you must list the name of a contig you wish to display for an alignment.\n"
              "Example: --contigs chrM chrX --chain_file=input.chain.liftover", file=sys.stderr)

    if args.output_name and args.chain_file and args.output_name[-1] != '_':
        args.output_name += '_'  # prefix should always end with an underscore

    # Set dependent defaults
    if not args.output_name and args.layout:
        if args.chain_file:
            args.output_name = 'Parallel_%s_and_%s_' % (just_the_name(args.fasta), just_the_name(args.extra_fastas[0]))
            if args.layout == "unique":
                args.output_name = '%s_unique_vs_%s_' % (just_the_name(args.fasta), just_the_name(args.extra_fastas[0]))
        else:
            either_name = args.fasta or args.image
            args.output_name = os.path.basename(os.path.splitext(either_name)[0])
    if args.output_name:
        args.output_name = args.output_name.strip()
    args.use_titles = not args.no_titles

    #Output directory: after args.output_name is set
    SERVER_HOME, base_path = base_directories(args.output_name)
    if args.quick:
        args.output_dir = os.path.dirname(
            os.path.abspath(args.fasta))  # just place the image next to the fasta
    else:
        args.output_dir = make_output_dir_with_suffix(base_path, '')

    ddv(args)


if __name__ == "__main__":
    main()
