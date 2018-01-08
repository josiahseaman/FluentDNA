#!/usr/bin/env python
"""
DDV 2.0 is a new version of DDV written in Python that allows you to generate a single image
for an entire genome.  It was necessary to switch platforms and languages because of intrinsic
limitations in the size of image that could be handled by: C#, DirectX, Win2D, GDI+, WIC, SharpDX,
or Direct2D. We tried a lot of options.

The python version has matured significantly past the previous feature set.

"""
from __future__ import print_function
from __future__ import division
from __future__ import unicode_literals
from __future__ import absolute_import
from __future__ import with_statement
from __future__ import generators
from __future__ import nested_scopes

import multiprocessing
import os
import sys


print("Setting up Python...")

if getattr(sys, 'frozen', False):
    BASE_DIR = os.path.dirname(sys.executable)
    os.environ["PATH"] += os.pathsep + os.path.join(BASE_DIR, 'bin')
    os.environ["PATH"] += os.pathsep + os.path.join(BASE_DIR, 'bin', 'env')
else:
    BASE_DIR = os.path.dirname(os.path.abspath(__file__))
print('Running in:', BASE_DIR)

sys.path.append(BASE_DIR)
sys.path.append(os.path.join(BASE_DIR, 'bin'))
sys.path.append(os.path.join(BASE_DIR, 'bin', 'env'))

os.chdir(BASE_DIR)

multiprocessing.freeze_support()

# ----------BEGIN MAIN PROGRAM----------
from DDV import VERSION

import shutil
import argparse

from DDV.DDVUtils import create_deepzoom_stack, make_output_dir_with_suffix, base_directories
from DNASkittleUtils.CommandLineUtils import just_the_name
from DDV.ParallelGenomeLayout import ParallelLayout
from DDV.ChainParser import ChainParser
from DDV.UniqueOnlyChainParser import UniqueOnlyChainParser
from DDV.AnnotatedAlignment import AnnotatedAlignment
from DDV.TileLayout import TileLayout
from DDV.TransposonLayout import TransposonLayout
from DDV.MultipleAlignmentLayout import MultipleAlignmentLayout


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

    print("Serving at %s on port %s" %(ADDRESS, str(PORT)))
    httpd.serve_forever()


def ddv(args):
    SERVER_HOME, base_path = base_directories(args)

    if not args.layout_type and args.run_server:
        run_server(SERVER_HOME)
        sys.exit(0)

    if args.ref_annotation and args.layout_type != 'transposon':  # parse chain files, possibly in batch
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
        batches = anno_align.parse_chain(args.chromosomes)
        del anno_align
        print("Done creating Gapped Annotations.")
        for batch in batches:  # multiple chromosomes, multiple views
            create_parallel_viz_from_fastas(args, len(batch.fastas), batch.output_folder, args.output_name,
                                            batch.fastas)
        sys.exit(0)

    if args.layout_type == "NONE":  # Complete webpage generation from existing image
        output_dir = make_output_dir_with_suffix(base_path, '')
        layout = TileLayout(use_titles=not args.no_titles, sort_contigs=args.sort_contigs,
                            high_contrast=args.high_contrast)
        layout.generate_html(base_path, output_dir, args.output_name)
        print("Creating Deep Zoom Structure for Existing Image...")
        create_deepzoom_stack(args.image, os.path.join(output_dir, 'GeneratedImages', "dzc_output.xml"))
        print("Done creating Deep Zoom Structure.")
        # TODO: Copy over html structure
        sys.exit(0)
    elif args.layout_type == "tiled":  # Typical Use Case
        # TODO: allow batch of tiling layout by chromosome
        if args.quick:
            output_dir = os.path.dirname(os.path.abspath(args.fasta))  # just place the image next to the fasta
        else:
            output_dir = make_output_dir_with_suffix(base_path, '')

        create_tile_layout_viz_from_fasta(args, args.fasta, output_dir, args.output_name)
        if args.run_server:
            run_server(output_dir)
        sys.exit(0)

    # views that support batches of chromosomes

    if args.layout_type == 'transposon':
        layout = TransposonLayout()
        output_dir = make_output_dir_with_suffix(base_path, '')
        # if len(args.chromosomes) != 1:
        #     raise NotImplementedError("Chromosome Argument requires exactly one chromosome e.g. '--chromosomes chr12'")
        layout.process_all_repeats(args.fasta, output_dir, just_the_name(output_dir), args.ref_annotation, args.chromosomes)
        print("Done with Transposons")
        sys.exit(0)
    if args.layout_type == 'alignment':
        output_dir = make_output_dir_with_suffix(base_path, '')
        layout = MultipleAlignmentLayout()
        layout.process_all_alignments(args.fasta,
                                      output_dir,
                                      args.output_name)
        print("Done with Alignments")
        sys.exit(0)
    if args.layout_type == "parallel":  # Parallel genome column layout OR quad comparison columns
        if not args.chain_file:  # life is simple
            # TODO: allow batch of tiling layout by chromosome
            # TODO: support drag and drop
            output_dir = make_output_dir_with_suffix(base_path, '')
            create_parallel_viz_from_fastas(args, len(args.extra_fastas) + 1, output_dir, args.output_name,
                                            [args.fasta] + args.extra_fastas)
            sys.exit(0)
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
                                       aligned_only=args.aligned_only)
            print("Creating Gapped and Unique Fastas from Chain File...")
            batches = chain_parser.parse_chain(args.chromosomes)
            del chain_parser
            print("Done creating Gapped and Unique.")
            for batch in batches:  # multiple chromosomes, multiple views
                create_parallel_viz_from_fastas(args, len(batch.fastas), batch.output_folder,
                                                args.output_name, batch.fastas)
            sys.exit(0)
    elif args.layout_type == "unique":
        """UniqueOnlyChainParser(chain_name='hg38ToPanTro4.over.chain',
                               first_source='HongKong\\hg38.fa',
                               second_source='',
                               output_folder_prefix='Hg38_unique_vs_panTro4_')"""
        unique_chain_parser = UniqueOnlyChainParser(chain_name=args.chain_file,
                                                    first_source=args.fasta,
                                                    second_source=args.fasta,
                                                    output_prefix=base_path,
                                                    trial_run=args.trial_run,
                                                    separate_translocations=args.separate_translocations)
        batches = unique_chain_parser.parse_chain(args.chromosomes)
        print("Done creating Gapped and Unique Fastas.")
        del unique_chain_parser
        for batch in batches:
            create_tile_layout_viz_from_fasta(args, batch.fastas[0], batch.output_folder, args.output_name)
        sys.exit(0)

    elif args.layout_type == "original":
        raise NotImplementedError("Original layout is not implemented!")
    else:
        raise NotImplementedError("What you are trying to do is not currently implemented!")


def create_parallel_viz_from_fastas(args, n_genomes, output_dir, output_name, fastas):
    print("Creating Large Comparison Image from Input Fastas...")
    layout = ParallelLayout(n_genomes=n_genomes)
    layout.process_file(output_dir, output_name, fastas)
    layout_final_output_location = layout.final_output_location
    del layout
    try:
        for extra_fasta in fastas:
            shutil.copy(extra_fasta, os.path.join(output_dir, os.path.basename(extra_fasta)))
    except shutil.Error:
        pass  # Same file is not a problem.  shutil.SameFileError is not defined in 2.7
    print("Done creating Large Image and HTML.")
    print("Creating Deep Zoom Structure from Generated Image...")
    create_deepzoom_stack(os.path.join(output_dir, layout_final_output_location),
                          os.path.join(output_dir, 'GeneratedImages', "dzc_output.xml"))
    print("Done creating Deep Zoom Structure.")

    if args.run_server:
        run_server(output_dir)



def create_tile_layout_viz_from_fasta(args, fasta, output_dir, output_name, layout=None):
    print("Creating Large Image from Input Fasta...")
    if layout is None:
        layout = TileLayout(use_titles=not args.no_titles, sort_contigs=args.sort_contigs,
                            high_contrast=args.high_contrast)
    layout.process_file(fasta, output_dir, output_name)
    layout_final_output_location = layout.final_output_location
    # try:
    #     shutil.copy(fasta, os.path.join(output_dir, os.path.basename(fasta)))
    # except shutil.SameFileError:
    #     pass  # not a problem
    print("Done creating Large Image at ", layout_final_output_location)
    if not args.no_webpage:
        layout.generate_html(fasta, output_dir, output_name)
        del layout
        print("Creating Deep Zoom Structure from Generated Image...")
        create_deepzoom_stack(os.path.join(output_dir, layout_final_output_location), os.path.join(output_dir, 'GeneratedImages', "dzc_output.xml"))
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
        sys.argv.append("--no_webpage")  # don't generate a full webpage (deepzoom is time consuming)
        sys.argv.append("--quick")
        # sys.argv.append("--sort_contigs")

    parser = argparse.ArgumentParser(usage="%(prog)s [options]",
                                     description="Creates visualizations of FASTA formatted DNA nucleotide data.",
                                     add_help=True)

    parser = argparse.ArgumentParser(prog='DDV.exe')
    parser.add_argument('--quick',
                        action='store_true',
                        help="Shortcut for dropping the file on DDV.exe.  Only an image will be generated "
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
    parser.add_argument('-s', '--sort_contigs',
                        action='store_true',
                        help="Sort the entries of the fasta file by length.  This option will kick in "
                             "automatically if your file has more than 10,000 separate FASTA entries.",
                        dest="sort_contigs")
    parser.add_argument('-hc', '--high_contrast',
                        action='store_true',
                        help="Use high contrast colors",
                        dest="high_contrast")
    parser.add_argument("-l", "--layout",
                        type=str,
                        help="The type of layout to perform. Will autodetect between Tiled and "
                        "Parallel. Really only need if you want the Original DDV layout or Unique only layout.",
                        choices=["tiled", "parallel", "alignment", "unique", "transposon", "original"],
                        dest="layout_type")  # Don't set a default so we can do error checking on it later
    parser.add_argument("-x", "--extrafastas",
                        nargs='+',
                        type=str,
                        help="Path to secondary FASTA files to process when doing Parallel layout.",
                        dest="extra_fastas")

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
    parser.add_argument("-c", "--chainfile",
                        type=str,
                        help="Path to Chain File when doing Parallel Comparisons layout.",
                        dest="chain_file")
    parser.add_argument("-ch", "--chromosomes",
                        nargs='+',
                        type=str,
                        help="Chromosome to parse from Chain File. NOTE: Defaults to 'chr21' for testing.",
                        dest="chromosomes")
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

    ### Other
    parser.add_argument("-i", "--image",
                        type=str,
                        help="Path to already computed big image to process with DeepZoom. "
                             "No layout will be performed if an image is passed in.",
                        dest="image")
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
    if args.layout_type == "original":
        parser.error("The 'original' layout is not yet implemented in Python!")  # TODO: Implement the original layout

    if args.image and (args.fasta or args.layout_type or args.extra_fastas or args.chain_file):
        parser.error("No layout will be performed if an existing image is passed in! Please only define an existing 'image' and the desired 'outfile'.")
    if not args.image and not args.fasta and not args.run_server:
        parser.error("Please either define a 'fasta' file or an 'image' file!")
    if args.image and args.no_webpage:
        parser.error("This parameter combination doesn't make sense.  You've provided a precalculated image"
                     "and asked DDV to only generate an image with no DeepZoom stack or webpage.")

    if args.extra_fastas and not args.layout_type:
        args.layout_type = "parallel"
    if args.layout_type and args.layout_type == "parallel" and not args.extra_fastas:
        parser.error("When doing a Parallel, you must at least define 'extrafastas'!")
    # if args.layout_type and args.layout_type == 'unique' and args.extra_fastas:
    #     parser.error("For Unique view, you don't need to specify 'extrafastas'.")
    if args.chromosomes and not (args.chain_file or args.layout_type == 'transposon'):
        parser.error("Listing 'Chromosomes' is only relevant when parsing Chain Files or Repeats!")
    # if args.extra_fastas and "parallel" not in args.layout_type:
    #     parser.error("The 'extrafastas' argument is only used when doing a Parallel layout!")
    if args.chain_file and args.layout_type not in ["parallel", "unique"]:
        parser.error("The 'chainfile' argument is only used when doing a Parallel or Unique layout!")
    if args.chain_file and len(args.extra_fastas) > 1:
        parser.error("Chaining more than two samples is currently not supported! Please only specify one 'extrafastas' when using a Chain input.")
    if args.layout_type == "unique" and not args.chain_file:
        parser.error("You must have a 'chainfile' to make a Unique layout!")
    if args.show_translocations_only and args.separate_translocations:
        parser.error("It just doesn't make sense to ask to show translocations in context while separating them.  You've got to pick one or the other.")

    # Set post error checking defaults
    if not args.image and not args.layout_type and args.fasta:
        args.layout_type = "tiled"
    if args.image and not args.layout_type:
        args.layout_type = "NONE"

    if args.chain_file and not args.chromosomes:
        args.chromosomes = ['chr21']

    if args.output_name and args.chain_file and args.output_name[-1] != '_':
        args.output_name += '_'  # prefix should always end with an underscore

    # Set dependent defaults
    if not args.output_name and args.layout_type:
        if args.chain_file:
            args.output_name = 'Parallel_%s_and_%s_' % (just_the_name(args.fasta), just_the_name(args.extra_fastas[0]))
            if args.layout_type == "unique":
                args.output_name = '%s_unique_vs_%s_' % (just_the_name(args.fasta), just_the_name(args.extra_fastas[0]))
        else:
            args.output_name = just_the_name(args.fasta or args.image)

    ddv(args)


if __name__ == "__main__":
    main()
