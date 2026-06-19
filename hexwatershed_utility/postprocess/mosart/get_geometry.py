
import sys
import numpy as np
from pyearth.gis.gdal.read.raster.gdal_read_envi_file import gdal_read_envi_file_multiple_band

from hexwatershed_utility.postprocess.mosart.find_contributing_cells import build_contributing_cells_all


def get_geometry(aLongitude_in, aLatitude_in, aCellID, aCellID_downslope, aArea,
                 pWidth_in=None, pDepth_in=None, iFlag_read_runoff=0, nChunk_size=1000):
    """
    Retrieve river width, depth and 2-year flood discharge for each mesh cell.

    Optimisations over the original implementation
    -----------------------------------------------
    1. Contributing-cell search: replaced O(N²) repeated BFS with a single
       O(N) topological sort via build_contributing_cells_all().
    2. Nearest-neighbour mapping: replaced Python loop + per-cell distance
       calculation with fully vectorised NumPy broadcasting (batched to
       control peak memory).
    3. Runoff extraction: pre-extract runoff for ALL cells at once so that
       the "outside-chunk re-read" anti-pattern is eliminated entirely.

    Args:
        aLongitude_in  : Array of cell longitudes, shape (N,)
        aLatitude_in   : Array of cell latitudes,  shape (N,)
        aCellID        : Array of cell IDs,         shape (N,)
        aCellID_downslope: Array of downstream cell IDs, shape (N,)
        aArea          : Array of cell areas (m²),  shape (N,)
        pWidth_in      : Width scaling coefficient (default 7.2).
        pDepth_in      : Depth scaling coefficient (default 0.27).
        iFlag_read_runoff: 1 = use runoff data; 0 = empirical area-based estimate.
        nChunk_size    : Batch size for vectorised nearest-neighbour mapping
                         (controls peak memory; default 1000).

    Returns:
        tuple: (aWidth_out, aDepth_out, aFlood_2yr_out)
    """
    nCell = len(aLongitude_in)

    pWidth = 7.2  if pWidth_in is None else pWidth_in
    pDepth = 0.27 if pDepth_in is None else pDepth_in

    # ------------------------------------------------------------------
    # 1. Build contributing-cell lists for ALL cells in O(N)
    # ------------------------------------------------------------------
    print('Building contributing-cell lists (O(N) topological sort)...')
    sys.stdout.flush()
    aCellIndex_contribution_all = build_contributing_cells_all(aCellID, aCellID_downslope)
    print(f'  Done. {nCell} cells processed.\n')

    # ------------------------------------------------------------------
    # 2. Runoff path
    # ------------------------------------------------------------------
    if iFlag_read_runoff == 1:
        nyear = 2009 - 1978
        nday  = nyear * 365

        # --- 2a. Runoff grid definition ---
        dResolution_runoff = 0.5
        aLongitude0 = np.arange(-179.75, 180,  dResolution_runoff)
        aLatitude0  = np.arange( -59.75,  90,  dResolution_runoff)
        aLongitude_grid, aLatitude_grid = np.meshgrid(aLongitude0, aLatitude0)
        nrow, ncol = aLongitude_grid.shape

        # --- 2b. Vectorised nearest-neighbour mapping (batched) ---
        # For each mesh cell find the (row, col) of the closest runoff grid point.
        # Processing in batches of nChunk_size avoids allocating an
        # (N × nrow × ncol) array all at once.
        print('Building nearest-neighbour mapping (vectorised, batched)...')
        sys.stdout.flush()

        row_idx = np.empty(nCell, dtype=np.int32)
        col_idx = np.empty(nCell, dtype=np.int32)

        nBatches = int(np.ceil(nCell / nChunk_size))
        for b in range(nBatches):
            i0 = b * nChunk_size
            i1 = min(i0 + nChunk_size, nCell)
            if b % 10 == 0:
                print(f'  Batch {b+1}/{nBatches}  (cells {i0}–{i1})...')

            lon_b = aLongitude_in[i0:i1]   # (B,)
            lat_b = aLatitude_in[i0:i1]    # (B,)

            # Squared angular distance — sufficient for argmin
            # shapes: (B,1,1) - (1,nrow,ncol) → (B,nrow,ncol)
            dist2 = ((lon_b[:, None, None] - aLongitude_grid[None, :, :]) ** 2 +
                     (lat_b[:, None, None] - aLatitude_grid[None, :, :])  ** 2)

            flat = np.argmin(dist2.reshape(i1 - i0, -1), axis=1)  # (B,)
            row_idx[i0:i1] = flat // ncol
            col_idx[i0:i1] = flat %  ncol

        print('  Nearest-neighbour mapping complete.\n')

        # --- 2c. Load runoff file once ---
        sWorkspace_runoff = '/compyfs/liao313/00raw'
        sFilename = sWorkspace_runoff + '/runoff.dat'
        print('Loading runoff data file...')
        sys.stdout.flush()
        dummy0 = gdal_read_envi_file_multiple_band(sFilename)
        aRunoff_source = dummy0[0]   # shape: (nday, nrow, ncol)

        # --- 2d. Extract per-cell runoff in one shot ---
        # aRunoff_all[t, i] = runoff at the grid point nearest to cell i on day t
        print('Extracting per-cell runoff (vectorised)...')
        sys.stdout.flush()
        aRunoff_all = aRunoff_source[:, row_idx, col_idx]   # (nday, nCell)
        aRunoff_all[aRunoff_all == -9999.0] = 0.0
        del aRunoff_source   # free the large source array

        # --- 2e. Compute discharge for every cell ---
        # Q[t, i] = sum over contributing cells k of (runoff[t,k] * area[k])
        # Unit conversion: mm/day → m³/s  (÷1000 for mm→m, ÷(3*3600) for day→s)
        print('Computing discharge for all cells...')
        sys.stdout.flush()
        aDischarge = np.zeros((nday, nCell), dtype=np.float64)
        for i in range(nCell):
            if i % 10000 == 0 and i > 0:
                print(f'  {i}/{nCell} cells...')
            contrib = aCellIndex_contribution_all[i]
            if len(contrib) == 0:
                continue
            contrib_arr = np.asarray(contrib, dtype=np.int64)
            # (nday, K) * (K,) → sum over K → (nday,)
            aDischarge[:, i] = np.dot(aRunoff_all[:, contrib_arr],
                                      aArea[contrib_arr])

        aDischarge /= (1000.0 * 3.0 * 3600.0)

        # --- 2f. Annual maximum flood (AMF) and 2-year return period ---
        print('Computing AMF and 2-year flood discharge...')
        sys.stdout.flush()
        AMF = np.zeros((nyear, nCell), dtype=np.float64)
        for y in range(nyear):
            AMF[y] = np.max(aDischarge[y * 365:(y + 1) * 365], axis=0)

        aFlood_2yr_out = np.percentile(AMF, 50, axis=0)

    else:
        # ------------------------------------------------------------------
        # Simplified path: empirical drainage-area relationship
        # Q_2yr (m³/s) ≈ 0.01 * Area(m²)^0.7
        # ------------------------------------------------------------------
        print('Skipping runoff data — using simplified area-based geometry...')
        aFlood_2yr_out = 0.01 * np.power(aArea, 0.7)

    # ------------------------------------------------------------------
    # 3. Width and depth from power-law scaling
    # ------------------------------------------------------------------
    print('Computing river width and depth...')
    aWidth_out = pWidth * np.power(aFlood_2yr_out, 0.52)
    aDepth_out = pDepth * np.power(aFlood_2yr_out, 0.31)

    print('Geometry calculation complete!\n')
    return aWidth_out, aDepth_out, aFlood_2yr_out
