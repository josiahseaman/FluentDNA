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
import shutil
import argparse

from DDVUtils import create_deepzoom_stack
from TileLayout import TileLayout
from ParallelGenomeLayout import ParallelLayout
from ChainParser import ChainParser


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
        output_image = args.output_name if args.output_name.lower().endswith(".png") else args.output_name + ".png"

        layout = TileLayout()
        layout.process_file(args.fasta_input, output_dir, output_image)
        create_deepzoom_stack(os.path.join(output_dir, output_image), os.path.join(output_dir, 'GeneratedImages', "dzc_output.xml"))
        sys.exit(0)
    elif args.layout_type == "parallel":  # Parallel genome column layout OR quad comparison columns
        output_image = args.output_name if args.output_name.lower().endswith(".png") else args.output_name + ".png"

        # layout = ParallelLayout(n_arguments - 3)
        # additional_files = argv[4:]
        # layout.process_file(input_file_path, folder, image, additional_files)
        create_deepzoom_stack(os.path.join(output_dir, output_image), os.path.join(output_dir, 'GeneratedImages', "dzc_output.xml"))
        sys.exit(0)
    elif args.layout_type == "original":
        raise NotImplementedError("Original layout is not implemented!")
    else:
        raise NotImplementedError("What you are trying to do is not currently implemented!")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(usage="%(prog)s [options]",
                                     description="Creates visualizations of FASTA formatted DNA nucleotide data.",
                                     add_help=True)

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
                        help="What to name the DeepZoom output DZI file (not a path).",
                        dest="output_name")
    parser.add_argument("-l", "--layout",
                        type=str,
                        help="The type of layout to perform. Will autodetect between Tiled and Parallel. Really only need if you want the Original DDV layout.",
                        choices=["original", "tiled", "parallel"],
                        dest="layout_type")  # Don't set a default so we can do error checking on it later
    parser.add_argument("-x", "--comparisonfasta",
                        nargs='+',
                        type=str,
                        help="Path to secondary FASTA file to process when doing Parallel Comparisons layout.",
                        dest="second_fasta")
    parser.add_argument("-c", "--chainfile",
                        type=str,
                        help="Path to Chain File when doing Parallel Comparisons layout.",
                        dest="chain_file")

    args = parser.parse_args()

    # Errors
    if args.layout_type == "original":
        parser.error("The 'original' layout is not yet implemented in Python!")  # TOOD: Implement the original layout

    if args.image and (args.input_fasta or args.layout_type or args.second_fasta or args.chain_file):
        parser.error("No layout will be performed if an existing image is passed in! Please only define an existing 'image' and the desired 'outfile'.")
    if not args.image and not args.input_fasta:
        parser.error("Please either define a 'fasta' file or an 'image' file!")

    if args.layout_type == "parallel" and not args.second_fasta:
        parser.error("When doing a Parallel layout, you must at least define 'comparisonfasta' if not 'comparisonfasta' and a 'chainfile'!")
    if args.second_fasta and not args.layout_type:
        args.layout_type = "parallel"
    if args.second_fasta and args.layout_type != "parallel":
        parser.error("The 'comparisonfasta' argument is only used when doing a Parallel layout!")
    if args.chain_file and args.layout_type != "parallel":
        parser.error("The 'chainfile' argument is only used when doing a Parallel layout!")

    # Set post error checking defaults
    if not args.image and not args.layout_type:
        args.layout_type = "tiled"

    # Set dependent defaults
    if not args.output_name:
        args.output_name = os.path.splitext(os.path.basename(args.input_fasta or args.image))[0]

    ddv(args)
