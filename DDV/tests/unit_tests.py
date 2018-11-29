import os
import unittest

from DNASkittleUtils.Contigs import read_contigs

from DDV.AnnotatedTrackLayout import AnnotatedTrackLayout

class AnnotationTrackTest(unittest.TestCase):
    def test_single_annotation_coord(self):
        print(os.getcwd())
        fa = "../example_data/gnetum_sample.fa"
        layout = AnnotatedTrackLayout(fa, "../tests/tiny_annotation.gff", annotation_width=18)
        layout.read_contigs_and_calc_padding(fa,[])
        layout.prepare_annotation_labels()


        self.assertEqual(True, True)


if __name__ == '__main__':
    unittest.main()
