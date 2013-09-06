"""Unit tests for Dataframe_pd.py"""

import dataFrame_pd as df
import unittest
import numpy as np


class CorrectInitialization(unittest.TestCase):
    def testCorrectInitialization(self):
        """The data in the data frame should match the data in the array_dict"""
        a = np.arange(100)
        b = np.arange(100)
        b = b + 10
        array_dict = {"item1":a,"item2":b}
        frame = df.Dataframe(array_dict = array_dict)
        for key, value in frame.data.items():
            self.assertEqual(frame.data[key].all(),array_dict[key].all())
        



if __name__ == "__main__":
    unittest.main()
    
