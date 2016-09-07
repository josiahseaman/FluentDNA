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


def ddv(argv):
    if getattr(sys, 'frozen', False):
        BASE_DIR = os.path.dirname(sys.executable)
    else:
        BASE_DIR = os.path.dirname(os.path.abspath(__file__))

    input_file_path = argv[1]
    chromosome_name = os.path.basename(input_file_path[:input_file_path.rfind(".")])  # name between /<path>/ and .png
    image = chromosome_name
    n_arguments = len(argv)

    if n_arguments == 2:  # Shortcut for old visualizations to create dz stack from existing large image
        output_file = chromosome_name + '.dzi'
        output_dir = os.path.dirname(input_file_path)
        create_deepzoom_stack(input_file_path, os.path.join(output_dir, output_file))
        sys.exit(0)

    if n_arguments == 3:
        raise ValueError("You need to specify an output folder and an image name")

    if n_arguments >= 4:  # Typical use case
        image = argv[3]
        chromosome_name = image  # based on image name, not fasta
        folder = os.path.join(argv[2], chromosome_name)  # place inside a folder with chromosome_name

    if n_arguments > 4:  # Multiple inputs => Parallel genome column layout
        pass
        layout = ParallelLayout(n_arguments - 3)
        additional_files = argv[4:]
        layout.process_file(input_file_path, folder, image, additional_files)
    else:  # Typical use case
        layout = TileLayout()
        layout.process_file(input_file_path, folder, image)

    create_deepzoom_stack(os.path.join(folder, image + '.png'), os.path.join(folder, 'GeneratedImages', 'dzc_output.xml'))
    print("Done creating Deep Zoom Structure\nCopying Source File:", input_file_path)
    destination = os.path.join(folder, os.path.basename(input_file_path))
    if not os.path.exists(destination):  # could have been created by ChainParser.py
        shutil.copy(input_file_path, destination)  # copy source file

    sys.exit(0)


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
                        help="What to name the DeepZoom output DZI file.",
                        dest="output_name")
    parser.add_argument("-l", "--layout",
                        type=str,
                        help="The type of layout to perform. Will autodetect between Tiled and Parallel. Really only need if you want the Original DDV layout.",
                        choices=["original", "tiled", "parallel"],
                        dest="layout_type")  # Don't set a default so we can do error checking on it later
    parser.add_argument("-x", "--extrafastas",
                        nargs='+',
                        type=str,
                        help="Path to secondary FASTA files to process when doing Parallel layout.",
                        dest="extra_fastas")
    parser.add_argument("-c", "--chainfile",
                        type=str,
                        help="Path to Chain File when doing Parallel layout.",
                        dest="chain_file")

    args = parser.parse_args()

    # Errors
    if args.image and (args.input_fasta or args.layout_type or args.extra_fastas or args.chain_file):
        parser.error("No layout will be performed if an existing image is passed in! Please only define an existing 'image' and the desired 'outfile'.")
    if not args.image and not args.input_fasta:
        parser.error("Please either define a 'fasta' file or an 'image' file!")

    if args.layout_type == "parallel" and not args.extra_fastas:
        parser.error("When doing a Parallel layout, you must at least define 'extrafastas' if not 'extrafastas' and a 'chainfile'!")
    if args.extra_fastas and not args.layout_type:
        args.layout_type = "parallel"
    if args.extra_fastas and args.layout_type != "parallel":
        parser.error("The 'extrafastas' argument is only used when doing a Parallel layout!")
    if args.chain_file and args.layout_type != "parallel":
        parser.error("The 'chainfile' argument is only used when doing a Parallel layout!")

    # Set post error checking defaults
    if not args.image and not args.layout_type:
        args.layout_type = "tiled"

    # Set dependent defaults
    if not args.output_name:
        args.output_name = os.path.basename(args.input_fasta)
    if not args.output_name.lower().endswith(".dzi"):
        args.output_name += ".dzi"

    ddv(sys.argv)  # TODO: Pass in the newly parsed args
