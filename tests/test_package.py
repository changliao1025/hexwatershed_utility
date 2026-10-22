import unittest
import hexwatershed_utility


class TestPackage(unittest.TestCase):
    def test_version(self):
        self.assertTrue(hasattr(hexwatershed_utility, "__version__"))
        self.assertIsInstance(hexwatershed_utility.__version__, str)

    def test_author(self):
        self.assertTrue(hasattr(hexwatershed_utility, "__author__"))
        self.assertEqual(hexwatershed_utility.__author__, "Chang Liao")

    def test_subpackages(self):
        import hexwatershed_utility.preprocess
        import hexwatershed_utility.postprocess
        import hexwatershed_utility.postprocess.land
        import hexwatershed_utility.postprocess.river
        self.assertIsNotNone(hexwatershed_utility.preprocess)
        self.assertIsNotNone(hexwatershed_utility.postprocess)
        self.assertIsNotNone(hexwatershed_utility.postprocess.land)
        self.assertIsNotNone(hexwatershed_utility.postprocess.river)


if __name__ == "__main__":
    unittest.main()
