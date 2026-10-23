import unittest

import numpy as np

from hexwatershed_utility.postprocess.river.find_contributing_cells import (
    build_contributing_cells_all,
    find_contributing_cells,
)


class TestContributingCells(unittest.TestCase):
    def setUp(self):
        # Directed acyclic graph (tree):
        # Cell 1 drains to Cell 3
        # Cell 2 drains to Cell 3
        # Cell 3 drains to Cell 4 (outlet)
        # Cell 4 drains to -1 (ocean/outlet)
        self.aCellID = np.array([1, 2, 3, 4])
        self.aCellID_downslope = np.array([3, 3, 4, -1])

    def test_find_contributing_cells_single(self):
        # Target = Cell 3
        cell_ids, cell_indices = find_contributing_cells(
            self.aCellID, self.aCellID_downslope, lCellID=3
        )
        self.assertIn(3, cell_ids)
        self.assertIn(1, cell_ids)
        self.assertIn(2, cell_ids)
        self.assertNotIn(4, cell_ids)
        self.assertEqual(len(cell_ids), 3)

    def test_find_contributing_cells_outlet(self):
        # Target = Cell 4 (outlet should have all cells contributing)
        cell_ids, cell_indices = find_contributing_cells(
            self.aCellID, self.aCellID_downslope, lCellID=4
        )
        self.assertEqual(set(cell_ids), {1, 2, 3, 4})

    def test_build_contributing_cells_all(self):
        contrib_all = build_contributing_cells_all(self.aCellID, self.aCellID_downslope)
        self.assertEqual(len(contrib_all), 4)

        # Index 0 (Cell 1): only itself [0]
        self.assertEqual(contrib_all[0], [0])

        # Index 1 (Cell 2): only itself [1]
        self.assertEqual(contrib_all[1], [1])

        # Index 2 (Cell 3): 3 itself + 1 and 2
        self.assertEqual(set(contrib_all[2]), {0, 1, 2})

        # Index 3 (Cell 4): all indices {0, 1, 2, 3}
        self.assertEqual(set(contrib_all[3]), {0, 1, 2, 3})


if __name__ == "__main__":
    unittest.main()
