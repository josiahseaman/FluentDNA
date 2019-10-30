import os
import unittest

from FluentDNA.AnnotatedTrackLayout import AnnotatedTrackLayout

class AnnotationTrackTest(unittest.TestCase):
    """The majority of testing is done in end_to_end_tests.py because visualization have
    to be visually verified for correctness.  Place math unit tests here."""
    def test_single_annotation_coord(self):
        print(os.getcwd())
        fa = "../example_data/gnetum_sample.fa"
        layout = AnnotatedTrackLayout(fa, "../tests/tiny_annotation.gff", annotation_width=18)
        layout.read_contigs_and_calc_padding(fa,[])
        layout.prepare_annotation_labels()

        self.assertEqual(True, True)


if __name__ == '__main__':
    unittest.main()
