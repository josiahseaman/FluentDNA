"""This is a convenient script to have for working with very large images
when you want to create high resolution, yet managable files for posters
at 300dpi, you'll want 7,000+ pixel width.  That's not possible through
screenshots and the 1GB image is usually too large to open.

You'll need to manually edit this script with you own paths and parameters."""
import sys
from os.path import dirname, join, splitext

from PIL import Image as PILImage
PILImage.MAX_IMAGE_PIXELS = 10 * 1024 * 1024 * 1024  # 10GB

def simple_resize(source, scale_factor):
    image = PILImage.open(source)
    print("Successfully opened image ", source)
    width, height = image.width, image.height
    width, height = int(width * scale_factor), int(height * scale_factor)
    print("Starting resize:", width, height)
    filter = PILImage.ANTIALIAS
    if scale_factor > 1:
        filter = PILImage.NEAREST  # don't blur pixel lines when upscaling
    resized = image.resize((width, height), filter)
    print("Successfully resized to", width, height)
    return resized

def save_resized(in_path, out_path, scale_factor):
    image = simple_resize(in_path, scale_factor)
    image.save(out_path, "png")

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Example downsizing:  python FluentDNA\image_resize_script.py /genomes/bigImage.png 1/6")
        sys.exit(1)
    # scale_factor = 1/2
    path = sys.argv[1]
    scale_factor = eval(sys.argv[2])
    out_path = splitext(path)[0] + ' - resized.png'
    save_resized(path, out_path, scale_factor)


    """
    path = r"D:\josiah\Projects\FluentDNA old\www-data\dnadata\Fraxinus pennsylvanica June 2017 sorted\Fraxinus pennsylvanica June 2017 sorted.png"
    path = r"D:\josiah\Projects\FluentDNA\FluentDNA\www-data\dnadata\Test Simple\Test Simple.png"
    path = r"D:\Genomes\Ash BATG-0.5-CLCbioSSPACE\Ash Tree BATG-0.5.png"
    path = r"D:\josiah\Documents\Research\Thesis - Genome Symmetry\Poster\Gene Families.png"
    path = r"D:\josiah\Documents\Research\Thesis - Genome Symmetry\Poster\example block.png"
    path = r"D:\josiah\Documents\Research\Thesis - Genome Symmetry\Bioinformatics Application Note\Figure 1 - manual colors.png"
    scale_factor = 4  # 1/2
    base = r'D:\josiah\Projects\FluentDNA\FluentDNA\www-data\dnadata\Oak\Oak - Quercus robur - chromosome '
    fb = r'\Oak - Quercus robur - chromosome '
    save_resized(base + '01' + fb + '01.png', "Chr01 - resized.png", scale_factor)
    save_resized(base + '02' + fb + '02.png', "Chr02 - resized.png", scale_factor)
    save_resized(base + '03' + fb + '03.png', "Chr03 - resized.png", scale_factor)
    save_resized(base + '04' + fb + '04.png', "Chr04 - resized.png", scale_factor)
    save_resized(base + '05' + fb + '05.png', "Chr05 - resized.png", scale_factor)
    save_resized(base + '06' + fb + '06.png', "Chr06 - resized.png", scale_factor)
    save_resized(base + '07' + fb + '07.png', "Chr07 - resized.png", scale_factor)
    save_resized(base + '08' + fb + '08.png', "Chr08 - resized.png", scale_factor)
    save_resized(base + '09' + fb + '09.png', "Chr09 - resized.png", scale_factor)
    save_resized(base + '10' + fb + '10.png', "Chr10 - resized.png", scale_factor)
    save_resized(base + '11' + fb + '11.png', "Chr11 - resized.png", scale_factor)
    save_resized(base + '12' + fb + '12.png', "Chr12 - resized.png", scale_factor)
    """