# Mantid Repository : https://github.com/mantidproject/mantid
#
# Copyright &copy; 2018 ISIS Rutherford Appleton Laboratory UKRI,
#   NScD Oak Ridge National Laboratory, European Spallation Source,
#   Institut Laue - Langevin & CSNS, Institute of High Energy Physics, CAS
# SPDX - License - Identifier: GPL - 3.0 +
import unittest
import numpy as np

from abins import KpointsData
from abins.constants import ACOUSTIC_PHONON_THRESHOLD


class KpointsDataTest(unittest.TestCase):

    _good_data_1 = {"k_vectors": np.asarray([[0.2, 0.1, 0.2], [0.1, 0.0, 0.2], [0.2, 0.2, 0.2]]),
                    "weights": np.asarray([0.3, 0.2, 0.5]),
                    "frequencies": np.asarray([[1.0, 2.0, 34.0, 4.9, 1.0, 2.0],
                                               [11.0, 12.0, 134.0, 14.9, 11.0, 12.0],
                                               [1.0, 2.0, 34.0, 4.9, 1.0, 2.0]]),  # 6 frequencies for one k-point
                    "atomic_displacements": np.asarray([[[[1.0, 1.0, 1.0], [1.0, 1.0, 1.0], [1.0, 1.0, 1.0],
                                                        [1.0, 1.0, 1.0], [1.0, 1.0, 1.0], [1.0, 1.0, 1.0]],

                                                       [[1.0, 1.0, 1.0], [1.0, 1.0, 111.0], [1.0, 1.0, 1.0],
                                                       [1.0, 1.0, 1.0], [1.0, 1.0, 1.0], [1.0, 1.0, 1.0]]],

                                                     [[[1.0, 1.0, 1.0], [1.0, 1.0, 1.0], [1.0, 1.0, 1.0],
                                                       [1.0, 1.0, 1.0], [1.0, 1.0, 1.0], [1.0, 1.0, 1.0]],

                                                      [[1.0, 1.0, 1.0], [1.0, 1.0, 221.0], [1.0, 1.0, 1.0],
                                                      [1.0, 1.0, 1.0], [1.0, 1.0, 1.0], [1.0, 1.0, 1.0]]],

                                                     [[[1.0, 1.0, 1.0], [1.0, 1.0, 41.0], [1.0, 1.0, 1.0],
                                                      [1.0, 1.0, 1.0], [1.0, 1.0, 31.0], [1.0, 1.0, 1.0]],

                                                      [[1.0, 1.0, 1.0], [1.0, 1.0, 1.0], [1.0, 1.0, 1.0],
                                                      [1.0, 1.0, 1.0], [1.0, 1.0, 41.0], [1.0, 1.0, 1.0]]]
                                                     # 12 atomic displacements for each k-point
                                                      ]).astype(complex),
                    "unit_cell": np.asarray([[ 7.44,  0.  ,  0.  ],
                                             [ 0.  ,  9.55,  0.  ],
                                             [ 0.  ,  0.  ,  6.92]])}

    # data with soft phonons
    _good_data_2 = {"k_vectors": np.asarray([[0.2, 0.1, 0.2], [0.1, 0.0, 0.2], [0.2, 0.2, 0.2]]),
                    "weights": np.asarray([0.3, 0.2, 0.5]),
                    "frequencies": np.asarray([[-10.0, -2.0, -3.0, 4.9, 1.0, 2.0],
                                               [11.0, 12.0, 134.0, 14.9, 11.0, 12.0],
                                               [1.0, 2.0, 34.0, 4.9, 1.0, 2.0]]),  # 6 frequencies for one k-point
                    "atomic_displacements": np.asarray([[[[1.0, 1.0, 1.0], [1.0, 1.0, 1.0], [1.0, 1.0, 1.0],
                                                          [1.0, 1.0, 1.0], [1.0, 121.0, 1.0], [1.0, 1.0, 1.0]],

                                                         [[1.0, 1.0, 1.0], [1.0, 1.0, 1.0], [1.0, 1.0, 1.0],
                                                          [1.0, 1.0, 1.0], [1.0, 1.0, 131.0], [1.0, 1.0, 1.0]]],


                                                        [[[1.0, 1.0, 1.0], [1.0, 1.0, 1.0], [1.0, 1.0, 1.0],
                                                          [1.0, 1.0, 1.0], [1.0, 1.0, 1.0], [1.0, 1.0, 1.0]],

                                                         [[1.0, 1.0, 1.0], [1.0, 1.0, 221.0], [1.0, 1.0, 1.0],
                                                          [1.0, 1.0, 1.0], [1.0, 1.0, 1.0], [1.0, 1.0, 1.0]]],


                                                        [[[1.0, 1.0, 1.0], [1.0, 1.0, 41.0], [1.0, 1.0, 1.0],
                                                          [1.0, 1.0, 1.0], [1.0, 1.0, 31.0], [1.0, 1.0, 1.0]],

                                                         [[1.0, 1.0, 1.0], [1.0, 1.0, 1.0], [1.0, 1.0, 1.0],
                                                          [1.0, 1.0, 1.0], [1.0, 1.0, 41.0], [1.0, 1.0, 1.0]]]
                                                        # 12 atomic displacements for each k-point
                                                        ]).astype(complex),
                    "unit_cell": np.asarray([[7.44, 0., 0.],
                                             [0., 9.55, 0.],
                                             [0., 0., 6.92]])
                    }

    def test_bad_items(self):
        # Dict has wrong contents
        with self.assertRaises(ValueError):
            wrong_dict = {"k_vectors": [], "freq": []}
            KpointsData.from_dict(wrong_dict)

        with self.assertRaises(TypeError):
            wrong_dict= ["k_vectors", 3, "freq"]
            KpointsData.from_dict(wrong_dict)

    def test_missing_key(self):
        # missing atomic_displacements
        items = self._good_data_1.copy()
        del items["atomic_displacements"]

        with self.assertRaises(ValueError):
            KpointsData.from_dict(items)

    def test_wrong_value(self):
        # All values should be numpy arrays
        items = self._good_data_1.copy()
        items["k_vectors"] = "wrong_value"

        with self.assertRaises(TypeError):
            KpointsData(**items)

    def test_wrong_weight(self):
        # negative weight (weight should be represented as a positive real number)
        items = self._good_data_1.copy()
        items["weights"] = np.asarray([-0.1, 0.3, 0.2])

        with self.assertRaises(ValueError):
            KpointsData(**items)

    def test_wrong_freq(self):
        # frequencies as a string
        wrong_items = self._good_data_1.copy()
        wrong_items["frequencies"] =  "Wrong_freq"

        with self.assertRaises(TypeError):
            KpointsData(**wrong_items)

        # complex frequencies
        wrong_items = self._good_data_1.copy()
        wrong_items["frequencies"] = wrong_items["frequencies"].astype(complex)

        with self.assertRaises(ValueError):
            KpointsData(**wrong_items)

        # frequencies as 2D arrays but with a bad shape
        wrong_items = self._good_data_1.copy()
        wrong_items["frequencies"] = np.asarray([[1.0, 2.0, 34.0], [4.9, 1.0, 1.0]])

        with self.assertRaises(ValueError):
            KpointsData(**wrong_items)

    def test_wrong_displacements(self):
        # displacements as a number
        wrong_items = self._good_data_1.copy()
        wrong_items["atomic_displacements"] = 1

        with self.assertRaises(TypeError):
            KpointsData(**wrong_items)

        # wrong size of the second dimension
        wrong_items = self._good_data_1.copy()
        wrong_items["atomic_displacements"] = np.asarray(
            [[[[1., 1., 11.],  [1.,  1., 1., 1.0], [1.0, 1.0, 1.0],
               [1., 1.0, 1.0], [1., 1., 11.],     [1., 1.,  11.]],
              [[1., 1.0, 1.0], [1., 1., 11.],     [1., 1.,  11.],
               [1., 1.0, 1.0],  [1., 1., 11.],     [1., 1.,  11.]]],
             wrong_items["atomic_displacements"][0, 0],
             wrong_items["atomic_displacements"][0, 1]])

        with self.assertRaises(ValueError):
            KpointsData(**wrong_items)

        # displacements as numpy arrays with integers
        wrong_items = self._good_data_1.copy()
        wrong_items["atomic_displacements"] = wrong_items["atomic_displacements"].astype(int)

        with self.assertRaises(ValueError):
            KpointsData(**wrong_items)

        # displacements as a 1D array
        wrong_items = self._good_data_1.copy()
        wrong_items["atomic_displacements"] = np.ravel(wrong_items["atomic_displacements"])

        with self.assertRaises(ValueError):
            KpointsData(**wrong_items)

    def test_set_good_case(self):

        self._set_good_case_core(data=self._good_data_1)
        self._set_good_case_core(data=self._good_data_2)

    def _set_good_case_core(self, data):
        kpd = KpointsData(**data)
        collected_data = kpd.extract()

        for k in range(data["frequencies"].shape[0]):
            indices = data["frequencies"][k] > ACOUSTIC_PHONON_THRESHOLD
            temp_f = data["frequencies"][k]
            self.assertEqual(True, np.allclose(temp_f[indices],
                                               collected_data["frequencies"][str(k)]))
            temp_a = data["atomic_displacements"][k]
            self.assertEqual(True, np.allclose(temp_a[:, indices],
                                               collected_data["atomic_displacements"][str(k)]))
            self.assertEqual(True, np.allclose(data["k_vectors"][k], collected_data["k_vectors"][str(k)]))
            self.assertEqual(data["weights"][k], collected_data["weights"][str(k)])


if __name__ == "__main__":
    unittest.main()
