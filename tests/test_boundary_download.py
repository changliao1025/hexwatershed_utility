import unittest
from unittest.mock import MagicMock, patch

from hexwatershed_utility.preprocess.feature.boundary import (
    download_boundary,
    download_conus_boundary,
)


class TestBoundaryDownload(unittest.TestCase):
    @patch("hexwatershed_utility.preprocess.feature.boundary.download_boundary.requests.get")
    def test_download_boundary_mock(self, mock_get):
        mock_response = MagicMock()
        mock_response.status_code = 200
        mock_response.json.return_value = {"type": "FeatureCollection", "features": []}
        mock_get.return_value = mock_response

        out_path = "/tmp/test_boundary_mock.geojson"
        res = download_boundary(sUrl="https://example.com/query", sFilename_geojson_out=out_path)
        self.assertEqual(res, out_path)


if __name__ == "__main__":
    unittest.main()
