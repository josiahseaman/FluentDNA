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

from DDVUtils import create_deepzoom_stack


def ddv(argv):
    from TileLayout import TileLayout
    from ParallelGenomeLayout import ParallelLayout

    folder = "."
    input_file_path = argv[1]
    chromosome_name = os.path.basename(input_file_path[:input_file_path.rfind(".")])  # name between /<path>/ and .png
    image = chromosome_name
    n_arguments = len(argv)

    if n_arguments == 2:  # Shortcut for old visualizations
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

    # sys.exit(0)


if __name__ == "__main__":
    ddv(sys.argv)
