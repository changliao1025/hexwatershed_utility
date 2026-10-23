import unittest

from hexwatershed_utility.core.exceptions import (
    DataInconsistencyError,
    HexWatershedError,
    MeshFormatError,
    TopologicalLoopError,
)


class TestExceptions(unittest.TestCase):
    def test_inheritance(self):
        self.assertTrue(issubclass(MeshFormatError, HexWatershedError))
        self.assertTrue(issubclass(TopologicalLoopError, HexWatershedError))
        self.assertTrue(issubclass(DataInconsistencyError, HexWatershedError))

    def test_raise_and_catch_base(self):
        with self.assertRaises(HexWatershedError):
            raise MeshFormatError("Invalid mesh format")

        with self.assertRaises(HexWatershedError):
            raise TopologicalLoopError("Detected topological cycle")

        with self.assertRaises(HexWatershedError):
            raise DataInconsistencyError("Array lengths do not match")

    def test_exception_messages(self):
        err = MeshFormatError("Corrupted file")
        self.assertEqual(str(err), "Corrupted file")


if __name__ == "__main__":
    unittest.main()
