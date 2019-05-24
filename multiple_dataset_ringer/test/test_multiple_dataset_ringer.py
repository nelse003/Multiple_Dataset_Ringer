import unittest
import shutil
import os

from multiple_dataset_ringer.multiple_dataset_ringer import (
    run as multiple_dataset_ringer,
)
from multiple_dataset_ringer.phil.multiple_dataset_phil import master_phil


class TestMultipleDatasetRinger(unittest.TestCase):
    """
    Test the main loop of multiple dataset ringer
    """

    def setUp(self):
        """
        Provide minimal parameters to test multiple dataset ringer
        """
        self.params = master_phil.extract()

        self.params.input.dir = ["multiple_dataset_ringer/test/resources/*/"]
        self.params.input.pdb_style = "*-pandda-input.pdb"
        self.params.input.mtz_style = "*_mrflagsref_idxs.mtz"
        self.params.output.out_dir = "output"
        self.params.input.column_labels = "2FOFCWT,PHI2FOFCWT"

    def tearDown(self):
        """
        Remove folder
        """
        shutil.rmtree(os.path.realpath("output"))


    def test_multiple_dataset_rigner(self):
        """
        Test multiple dataset ringer
        """
        multiple_dataset_ringer(self.params)

        # TODO Add a more through set of file existence to test

        self.assertTrue(os.path.isfile(os.path.join("output", "ARG A24","all-ARG A24-257-dataset.png")), True)


if __name__ == "__main__":
    unittest.main()
