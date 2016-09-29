"""
self.is an optional addon to DDV written in Python that allows you to generate a single image
for an entire genome.  It was necessary to switch platforms and languages because of intrinsic
limitations in the size of image that could be handled by: C#, DirectX, Win2D, GDI+, WIC, SharpDX,
or Direct2D. We tried a lot of options.

self.python file contains basic image handling methods.  It also contains a re-implementation of
Josiah's "Tiled Layout" algorithm which is also in DDVLayoutManager.cs.
"""
import os
import sys
import multiprocessing


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
__version__ = '1.0.0'

import shutil
import argparse
import subprocess

# from http import server

from DDVUtils import create_deepzoom_stack
from TileLayout import TileLayout
from ParallelGenomeLayout import ParallelLayout
from ChainParser import ChainParser
from UniqueOnlyChainParser import UniqueOnlyChainParser


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


def ddv(args):
    if getattr(sys, 'frozen', False):
        BASE_DIR = os.path.dirname(sys.executable)
    else:
        BASE_DIR = os.path.dirname(os.path.abspath(__file__))

    output_dir = os.path.join(BASE_DIR, "www-data", "dnadata", args.output_name)

    print("Creating Chromosome Output Directory...")
    os.makedirs(output_dir, exist_ok=True)
    print("Done creating Directories.")

    if not args.layout_type:  # Shortcut for old visualizations to create dz stack from existing large image
        print("Creating Deep Zoom Structure for Existing Image...")
        create_deepzoom_stack(args.image, os.path.join(output_dir, 'GeneratedImages', "dzc_output.xml"))
        shutil.copy(args.image, os.path.join(output_dir, os.path.basename(args.image)))
        print("Done creating Deep Zoom Structure.")
        # TODO: Copy over html structure
        sys.exit(0)
    elif args.layout_type == "tiled":  # Typical Use Case
        print("Creating Large Image from Input Fasta...")
        layout = TileLayout()
        layout.process_file(args.input_fasta, output_dir, args.output_name)
        shutil.copy(args.input_fasta, os.path.join(output_dir, os.path.basename(args.input_fasta)))
        print("Done creating Large Image and HTML.")

        print("Creating Deep Zoom Structure from Generated Image...")
        create_deepzoom_stack(os.path.join(output_dir, layout.final_output_location), os.path.join(output_dir, 'GeneratedImages', "dzc_output.xml"))
        print("Done creating Deep Zoom Structure.")

        if args.run_server:
            pass
            # handler_class = server.BaseHTTPRequestHandler
            # server.test(HandlerClass=handler_class, port=8000, bind='')

        sys.exit(0)
    elif args.layout_type == "parallel":  # Parallel genome column layout OR quad comparison columns
        n_genomes = len(args.extra_fastas) + 1

        if args.chain_file:
            print("Creating Gapped and Unique Fastas from Chain File...")
            if args.layout_type == "parallel":
                chain_parser = ChainParser(chain_name=args.chain_file,
                                           first_source=args.input_fasta,
                                           second_source=args.extra_fastas[0],
                                           output_folder=output_dir,
                                           trial_run=False)
                chain_parser.parse_chain(args.chromosomes)
                n_genomes = 4
                swap_columns = True  # TODO: Make this variable
                if swap_columns:
                    args.extra_fastas = chain_parser.output_fastas.reverse()
                else:
                    args.extra_fastas = chain_parser.output_fastas
                args.fasta = args.extra_fastas.pop()
            elif args.layout_type == "unique-parallel":
                unique_chain_parser = UniqueOnlyChainParser(chain_name=args.chain_file,
                                                            first_source=args.input_fasta,
                                                            second_source=args.extra_fastas[0],
                                                            output_folder=output_dir,
                                                            trial_run=False)
                unique_chain_parser.parse_chain(args.chromosomes)
                n_genomes = 2
                args.extra_fastas = unique_chain_parser.output_fastas
                args.fasta = args.extras.pop()
                # TODO: Does this need to be sent to NOT Parallel Layout?
            print("Done creating Gapped and Unique Fastas.")

        print("Creating Large Comparison Image from Input Fastas...")
        layout = ParallelLayout(n_genomes=n_genomes)
        layout.process_file(args.input_fasta, output_dir, args.output_name, args.extra_fastas)
        shutil.copy(args.input_fasta, os.path.join(output_dir, os.path.basename(args.input_fasta)))
        for extra_fasta in args.extra_fastas:
            shutil.copy(extra_fasta, os.path.join(output_dir, os.path.basename(extra_fasta)))
        print("Done creating Large Image and HTML.")

        print("Creating Deep Zoom Structure from Generated Image...")
        create_deepzoom_stack(os.path.join(output_dir, layout.final_output_location), os.path.join(output_dir, 'GeneratedImages', "dzc_output.xml"))
        print("Done creating Deep Zoom Structure.")

        if args.run_server:
            pass
            # handler_class = server.BaseHTTPRequestHandler
            # server.test(HandlerClass=handler_class, port=8000, bind='')

        sys.exit(0)
    elif args.layout_type == "original":
        raise NotImplementedError("Original layout is not implemented!")
    else:
        raise NotImplementedError("What you are trying to do is not currently implemented!")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(usage="%(prog)s [options]",
                                     description="Creates visualizations of FASTA formatted DNA nucleotide data.",
                                     add_help=True)

    parser = argparse.ArgumentParser(prog='DDV.exe')
    parser.add_argument('-n', '--update_name', dest='update_name', help='Query for the name of this program as known to the update server', action='store_true')
    parser.add_argument('-v', '--version', dest='version', help='Get current version of program.', action='store_true')

    parser.add_argument("-i", "--image",
                        type=str,
                        help="Path to already laid out big image to process with DeepZoom. No layout will be performed if an image is passed in.",
                        dest="image")
    parser.add_argument("-f", "--fasta",
                        type=str,
                        help="Path to main FASTA file to process into new visualization.",
                        dest="input_fasta")
    parser.add_argument("-o", "--outname",
                        type=str,
                        help="What to name the output folder (not a path). Defaults to name of the fasta file.",
                        dest="output_name")
    parser.add_argument("-l", "--layout",
                        type=str,
                        help="The type of layout to perform. Will autodetect between Tiled and Parallel. Really only need if you want the Original DDV layout.",
                        choices=["original", "tiled", "parallel", "unique-parallel"],
                        dest="layout_type")  # Don't set a default so we can do error checking on it later
    parser.add_argument("-x", "--extrafastas",
                        nargs='+',
                        type=str,
                        help="Path to secondary FASTA files to process when doing Parallel layout.",
                        dest="extra_fastas")
    parser.add_argument("-c", "--chainfile",
                        type=str,
                        help="Path to Chain File when doing Parallel Comparisons layout.",
                        dest="chain_file")
    parser.add_argument("-ch", "--chromosomes",
                        nargs='+',
                        type=str,
                        help="Chromosome to parse from Chain File. NOTE: Defaults to 'chrY' for testing.",
                        dest="chromosomes")
    parser.add_argument("-s", "--server",
                        action='store_true',
                        help="Run Web Server after computing.",
                        dest="run_server")

    args = parser.parse_args()

    # Respond to an updater query
    if args.update_name:
        print("DDV")
        sys.exit(0)
    elif args.version:
        print(__version__)
        sys.exit(0)

    # Errors
    if args.layout_type == "original":
        parser.error("The 'original' layout is not yet implemented in Python!")  # TOOD: Implement the original layout

    if args.image and (args.input_fasta or args.layout_type or args.extra_fastas or args.chain_file):
        parser.error("No layout will be performed if an existing image is passed in! Please only define an existing 'image' and the desired 'outfile'.")
    if not args.image and not args.input_fasta:
        parser.error("Please either define a 'fasta' file or an 'image' file!")

    if args.extra_fastas and not args.layout_type:
        args.layout_type = "parallel"
    if "parallel" in args.layout_type and not args.extra_fastas:
        parser.error("When doing a Parallel layout, you must at least define 'extrafastas' if not 'extrafastas' and a 'chainfile'!")
    if args.chromosomes and not args.chain_file:
        parser.error("Listing 'Chromosomes' is only relevant when parsing Chain Files!")
    if args.extra_fastas and "parallel" not in args.layout_type:
        parser.error("The 'extrafastas' argument is only used when doing a Parallel layout!")
    if args.chain_file and "parallel" not in args.layout_type:
        parser.error("The 'chainfile' argument is only used when doing a Parallel layout!")
    if args.chain_file and len(args.extra_fastas) > 1:
        parser.error("Chaining more than two samples is currently not supported! Please only specify one 'extrafastas' when using a Chain input.")

    # Set post error checking defaults
    if not args.image and not args.layout_type:
        args.layout_type = "tiled"

    if args.chain_file and not args.chromosomes:
        args.chromosomes = ['chrY']

    # Set dependent defaults
    if not args.output_name:
        if args.chain_file:
            args.output_name = os.path.splitext(os.path.basename(args.input_fasta))[0] + '_AND_' + os.path.splitext(os.path.basename(args.extra_fastas[0]))[0]
            if args.layout_type == "unique-parallel":
                args.output_name += "_UNIQUE"
        else:
            args.output_name = os.path.splitext(os.path.basename(args.input_fasta or args.image))[0]

    ddv(args)
