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
import psutil
import subprocess

from http import server

from DDVUtils import create_deepzoom_stack
from TileLayout import TileLayout
from ParallelGenomeLayout import ParallelLayout
from ChainParser import ChainParser


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
            handler_class = server.BaseHTTPRequestHandler
            server.test(HandlerClass=handler_class, port=8000, bind='')

        sys.exit(0)
    elif args.layout_type == "parallel":  # Parallel genome column layout OR quad comparison columns
        n_genomes = len(args.extra_fastas) + 1

        if args.chain_file:
            print("Created Gapped and Unique Fastas from Chain File...")
            chain_parser = ChainParser(args.input_fasta, args.extra_fastas[0], args.chain_file, output_dir)
            chain_parser.parse_chain(args.chromosomes)
            n_genomes = 4
            args.extra_fastas = chain_parser.extra_generated_fastas
            args.fasta = args.extra_fastas.pop()
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
            handler_class = server.BaseHTTPRequestHandler
            server.test(HandlerClass=handler_class, port=8000, bind='')

        sys.exit(0)
    elif args.layout_type == "original":
        raise NotImplementedError("Original layout is not implemented!")
    else:
        raise NotImplementedError("What you are trying to do is not currently implemented!")


def launch_external_program_and_exit(launch, code=0, close_self=True, cmd_args=None, launch_args=None):
    if not launch_args:
        launch_args = {}
    if not cmd_args:
        cmd_args = []
    launch = [launch, ]
    if cmd_args:
        for cmd_arg in cmd_args:
            launch.append(cmd_arg)
    launch = ' '.join(launch)
    if sys.platform == 'win32':  # Yes, this is also x64.
        CREATE_NEW_PROCESS_GROUP = 0x00000200
        DETACHED_PROCESS = 0x00000008
        launch_args.update(creationflags=DETACHED_PROCESS | CREATE_NEW_PROCESS_GROUP)
    else:
        launch_args.update(preexec_fn=os.setsid)
        launch_args.update(start_new_session=True)
    subprocess.Popen(launch, stdin=subprocess.PIPE, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, **launch_args)
    if close_self:
        sys.exit(code)


def check_update():
    new_version = None
    try:
        npu = os.path.join(BASE_DIR, 'npu'+EXTENSION)
        process = subprocess.Popen(npu + " --check_update --silent", shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        try:
            output, error = process.communicate(timeout=60000)
            exit_code = process.returncode
        except:
            exit_code = 1
            output = None
        try:
            process.kill()
        except:
            pass

        if output:
            new_version = output.splitlines()[-1].decode().strip()
    except:
        print("Unable to get DDV Version!")

    return new_version


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
                        help="What to name the DeepZoom output DZI file (not a path).",
                        dest="output_name")
    parser.add_argument("-l", "--layout",
                        type=str,
                        help="The type of layout to perform. Will autodetect between Tiled and Parallel. Really only need if you want the Original DDV layout.",
                        choices=["original", "tiled", "parallel"],
                        dest="layout_type")  # Don't set a default so we can do error checking on it later
    parser.add_argument("-x", "--extrafastas",
                        nargs='+',
                        type=str,
                        help="Path to secondary FASTA file to process when doing Parallel Comparisons layout.",
                        dest="extra_fastas")
    parser.add_argument("-c", "--chainfile",
                        type=str,
                        help="Path to Chain File when doing Parallel Comparisons layout.",
                        dest="chain_file")
    parser.add_argument("-ch", "--chromosomes",
                        nargs='+',
                        type=str,
                        help="Chromosome to parse from Chain File. NOTE: Defaults to 'chr21' for testing.",
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

    # Check that another instance of the program isn't running
    for proc in psutil.process_iter():
        try:
            proc_name = proc.name().lower()
        except psutil.AccessDenied as e:
            continue
        if 'DDV'.lower() in proc_name:  # TODO: Is DDV actually the process name? Hopefully...
            print("\nThere is already an instance of ADSM running!")
            print("\nPress any key to exit...")
            input()
            sys.exit(1)

    if getattr(sys, 'frozen', False):
        print("Checking for updates...")
        version = check_update()
        if version and version != 'False' and version != '0':
            print("New version available:", version)
            if query_yes_no("Do you want to update now?"):
                try:
                    print("Launching Newline Program Updater...")
                    npu = os.path.join(BASE_DIR, 'npu'+EXTENSION)
                    launch_external_program_and_exit(npu)#cmd_args=['--silent'])
                except:
                    print("Failed to run the Program Updater!")
            else:
                print("Skipping update.")
        else:
            print("No updates currently.")

    # Errors
    if args.layout_type == "original":
        parser.error("The 'original' layout is not yet implemented in Python!")  # TOOD: Implement the original layout

    if args.image and (args.input_fasta or args.layout_type or args.extra_fastas or args.chain_file):
        parser.error("No layout will be performed if an existing image is passed in! Please only define an existing 'image' and the desired 'outfile'.")
    if not args.image and not args.input_fasta:
        parser.error("Please either define a 'fasta' file or an 'image' file!")

    if args.layout_type == "parallel" and not args.extra_fastas:
        parser.error("When doing a Parallel layout, you must at least define 'extrafastas' if not 'extrafastas' and a 'chainfile'!")
    if args.extra_fastas and not args.layout_type:
        args.layout_type = "parallel"
    if args.chromosomes and not args.layout_type == "parallel" or not args.chain_file:
        parser.error("Listing 'Chromosomes' is only relevant when parsing Chain Files!")
    if args.extra_fastas and args.layout_type != "parallel":
        parser.error("The 'extrafastas' argument is only used when doing a Parallel layout!")
    if args.chain_file and args.layout_type != "parallel":
        parser.error("The 'chainfile' argument is only used when doing a Parallel layout!")
    if args.chain_file and len(args.extra_fastas) > 1:
        parser.error("Chaining more than two samples is currently not supported! Please only specify one 'extrafastas' when using a Chain input.")

    # Set post error checking defaults
    if not args.image and not args.layout_type:
        args.layout_type = "tiled"

    if args.chain_file and not args.chromosomes:
        args.chromosomes = ['chr21']

    # Set dependent defaults
    if not args.output_name:
        if args.chain_file:
            args.output_name = os.path.splitext(os.path.basename(args.input_fasta))[0] + '_AND_' + os.path.splitext(os.path.basename(args.extra_fastas[0]))[0]
        else:
            args.output_name = os.path.splitext(os.path.basename(args.input_fasta or args.image))[0]

    ddv(args)
