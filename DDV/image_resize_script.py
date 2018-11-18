"""This is a convenient script to have for working with very large images
when you want to create high resolution, yet managable files for posters
at 300dpi, you'll want 7,000+ pixel width.  That's not possible through
screenshots and the 1GB image is usually too large to open.

You'll need to manually edit this script with you own paths and parameters."""


from PIL import Image as PILImage


def simple_resize(source, reduction_factor):
    image = PILImage.open(source)
    print("Successfully opened image ", source)
    width, height = image.width, image.height
    width, height = int(width // reduction_factor), int(height // reduction_factor)
    print("Starting resize:", width, height)
    filter = PILImage.ANTIALIAS
    if reduction_factor < 1:
        filter = PILImage.NEAREST  # don't blur pixel lines when upscaling
    resized = image.resize((width, height), filter)
    print("Successfully resized to", width, height)
    return resized

def save_resized(in_path, out_path, reduction_factor):
    image = simple_resize(in_path, reduction_factor)
    image.save(out_path, "png")

if __name__ == "__main__":
    path = r"D:\josiah\Projects\DDV old\www-data\dnadata\Fraxinus pennsylvanica June 2017 sorted\Fraxinus pennsylvanica June 2017 sorted.png"
    path = r"D:\josiah\Projects\DDV\DDV\www-data\dnadata\Test Simple\Test Simple.png"
    path = r"D:\Genomes\Ash BATG-0.5-CLCbioSSPACE\Ash Tree BATG-0.5.png"
    path = r"D:\josiah\Documents\Research\Thesis - Genome Symmetry\Poster\Gene Families.png"
    path = r"D:\josiah\Documents\Research\Thesis - Genome Symmetry\Poster\example block.png"
    path = r"D:\josiah\Documents\Research\Thesis - Genome Symmetry\Bioinformatics Application Note\Figure 1 - manual colors.png"
    reduction_factor = 4  # 1/2
    base = r'D:\josiah\Projects\DDV\DDV\www-data\dnadata\Oak\Oak - Quercus robur - chromosome '
    fb = r'\Oak - Quercus robur - chromosome '
    save_resized(base + '01' + fb + '01.png', "Chr01 - resized.png", reduction_factor)
    save_resized(base + '02' + fb + '02.png', "Chr02 - resized.png", reduction_factor)
    save_resized(base + '03' + fb + '03.png', "Chr03 - resized.png", reduction_factor)
    save_resized(base + '04' + fb + '04.png', "Chr04 - resized.png", reduction_factor)
    save_resized(base + '05' + fb + '05.png', "Chr05 - resized.png", reduction_factor)
    save_resized(base + '06' + fb + '06.png', "Chr06 - resized.png", reduction_factor)
    save_resized(base + '07' + fb + '07.png', "Chr07 - resized.png", reduction_factor)
    save_resized(base + '08' + fb + '08.png', "Chr08 - resized.png", reduction_factor)
    save_resized(base + '09' + fb + '09.png', "Chr09 - resized.png", reduction_factor)
    save_resized(base + '10' + fb + '10.png', "Chr10 - resized.png", reduction_factor)
    save_resized(base + '11' + fb + '11.png', "Chr11 - resized.png", reduction_factor)
    save_resized(base + '12' + fb + '12.png', "Chr12 - resized.png", reduction_factor)
