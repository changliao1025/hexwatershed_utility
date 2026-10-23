import os
import tempfile
import unittest

import numpy as np

from hexwatershed_utility.core.graph import (
    build_downstream_map,
    find_downstream_cells,
    meshjson_to_connectivity,
)
from tests.fixtures.synthetic_data import create_synthetic_mesh_json


class TestGraph(unittest.TestCase):
    def test_build_downstream_map(self):
        cell_ids = [10, 20, 30]
        downslope_ids = [20, 30, -1]
        dn_map = build_downstream_map(cell_ids, downslope_ids)

        self.assertEqual(dn_map, {10: 20, 20: 30, 30: -1})
        self.assertEqual(build_downstream_map([], []), {})

    def test_meshjson_to_connectivity(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            json_path = os.path.join(tmpdir, "test_mesh.json")
            create_synthetic_mesh_json(json_path, n_cells=7)

            id_arr, dn_arr = meshjson_to_connectivity(json_path)

            self.assertEqual(len(id_arr), 7)
            self.assertEqual(len(dn_arr), 7)
            np.testing.assert_array_equal(id_arr, np.arange(1, 8))

            # Cell 1 (index 0) drains to cell 2 (1-based index 2)
            self.assertEqual(dn_arr[0], 2)
            # Cell 2 (index 1) drains to cell 3 (1-based index 3)
            self.assertEqual(dn_arr[1], 3)
            # Cell 7 (index 6) is outlet -> -9999
            self.assertEqual(dn_arr[6], -9999)

    def test_find_downstream_cells_mosart(self):
        # 1 -> 2 -> 3 -> -9999
        cell_ids = [1, 2, 3]
        downslope_ids = [2, 3, -9999]

        indices, ids = find_downstream_cells(
            cell_ids=cell_ids,
            downslope_ids=downslope_ids,
            query_indices=[0],
            id_type="mosart",
        )

        self.assertEqual(indices[0], [1, 2])
        self.assertEqual(ids[0], [2, 3])

    def test_find_downstream_cells_hexwatershed(self):
        # 10 -> 20 -> -1
        cell_ids = [10, 20]
        downslope_ids = [20, -1]

        indices, ids = find_downstream_cells(
            cell_ids=cell_ids,
            downslope_ids=downslope_ids,
            query_indices=[0],
            id_type="hexwatershed",
        )

        self.assertEqual(indices[0], [1])
        self.assertEqual(ids[0], [20])

    def test_find_downstream_cells_cycle_protection(self):
        # 1 -> 2 -> 1 (cycle)
        cell_ids = [1, 2]
        downslope_ids = [2, 1]

        indices, ids = find_downstream_cells(
            cell_ids=cell_ids,
            downslope_ids=downslope_ids,
            query_indices=[0],
            id_type="mosart",
            max_hops=10,
        )

        # Loop detected and broke cleanly without hanging
        self.assertTrue(len(indices[0]) <= 2)


if __name__ == "__main__":
    unittest.main()
