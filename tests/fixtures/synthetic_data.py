"""
Synthetic test fixtures generator for HexWatershed Utility.
"""

import csv
import json

import netCDF4 as nc
import numpy as np


def create_synthetic_mesh_json(file_path: str, n_cells: int = 7) -> str:
    """
    Create a synthetic HexWatershed JSON file representing a simple river mesh.
    Default: 7 cells with dendritic drainage towards cell 7 (outlet).
    """
    cells = []
    # Base centers
    lons = [-100.0, -99.9, -99.8, -99.7, -99.9, -99.8, -99.6]
    lats = [40.0, 40.0, 40.0, 40.0, 40.1, 40.1, 40.0]
    downslopes = [2, 3, 4, 7, 3, 4, -1]
    elevations = [150.0, 140.0, 130.0, 120.0, 145.0, 135.0, 110.0]
    drainage_areas = [1.0e7, 2.5e7, 5.0e7, 8.5e7, 1.2e7, 1.8e7, 1.2e8]

    # Hexagon offsets (~0.04 degrees)
    r = 0.04
    angles = np.linspace(0, 2 * np.pi, 7)[:-1]

    for i in range(min(n_cells, len(lons))):
        cid = i + 1
        cx, cy = lons[i], lats[i]

        vertices = []
        for angle in angles:
            vx = cx + r * np.cos(angle)
            vy = cy + r * np.sin(angle)
            vertices.append({"dLongitude_degree": float(vx), "dLatitude_degree": float(vy)})

        cell_data = {
            "lCellID": cid,
            "lCellID_downslope": downslopes[i],
            "dLongitude_center_degree": cx,
            "dLatitude_center_degree": cy,
            "dElevation": elevations[i],
            "dDrainage_area": drainage_areas[i],
            "dSlope_between": 0.002,
            "dLength_flowline": 5000.0,
            "dArea": 1.0e7,
            "vVertex": vertices,
        }
        cells.append(cell_data)

    with open(file_path, "w") as fp:
        json.dump(cells, fp, indent=2)

    return file_path


def create_synthetic_mosart_files(domain_path: str, parameter_path: str, n_cells: int = 7) -> None:
    """
    Create minimal valid domain and parameter NetCDF files for MOSART testing.
    """
    lons = np.array([-100.0, -99.9, -99.8, -99.7, -99.9, -99.8, -99.6][:n_cells])
    lats = np.array([40.0, 40.0, 40.0, 40.0, 40.1, 40.1, 40.0][:n_cells])
    cell_ids = np.arange(1, n_cells + 1, dtype=np.int32)
    ids = np.arange(1, n_cells + 1, dtype=np.int32)
    # Downstream IDs (-9999 for outlet)
    dn_ids = np.array([2, 3, 4, 7, 3, 4, -9999][:n_cells], dtype=np.int32)
    area_accum = np.array([10.0, 25.0, 50.0, 85.0, 12.0, 18.0, 120.0][:n_cells]) * 1.0e6

    # 1. Domain NetCDF
    nv = 6
    r = 0.04
    angles = np.linspace(0, 2 * np.pi, nv, endpoint=False)

    xv = np.zeros((n_cells, 1, nv), dtype=np.float64)
    yv = np.zeros((n_cells, 1, nv), dtype=np.float64)

    for i in range(n_cells):
        for j, ang in enumerate(angles):
            xv[i, 0, j] = lons[i] + r * np.cos(ang)
            yv[i, 0, j] = lats[i] + r * np.sin(ang)

    with nc.Dataset(domain_path, "w", format="NETCDF4") as ds_dom:
        ds_dom.createDimension("gridcell", n_cells)
        ds_dom.createDimension("ni", 1)
        ds_dom.createDimension("nv", nv)

        v_xc = ds_dom.createVariable("xc", "f8", ("gridcell",))
        v_yc = ds_dom.createVariable("yc", "f8", ("gridcell",))
        v_xv = ds_dom.createVariable("xv", "f8", ("gridcell", "ni", "nv"))
        v_yv = ds_dom.createVariable("yv", "f8", ("gridcell", "ni", "nv"))
        v_area = ds_dom.createVariable("area", "f8", ("gridcell",))
        v_mask = ds_dom.createVariable("mask", "i4", ("gridcell",))

        v_xc[:] = lons
        v_yc[:] = lats
        v_xv[:] = xv
        v_yv[:] = yv
        v_area[:] = np.full(n_cells, 1.0e7)
        v_mask[:] = np.ones(n_cells, dtype=np.int32)

    # 2. Parameter NetCDF
    with nc.Dataset(parameter_path, "w", format="NETCDF4") as ds_par:
        ds_par.createDimension("gridcell", n_cells)

        v_cellid = ds_par.createVariable("CellID", "i4", ("gridcell",))
        v_id = ds_par.createVariable("ID", "i4", ("gridcell",))
        v_dnid = ds_par.createVariable("dnID", "i4", ("gridcell",))
        v_fdir = ds_par.createVariable("fdir", "i4", ("gridcell",))
        v_lat = ds_par.createVariable("latixy", "f8", ("gridcell",))
        v_lon = ds_par.createVariable("longxy", "f8", ("gridcell",))
        v_area = ds_par.createVariable("area", "f8", ("gridcell",))
        v_area2 = ds_par.createVariable("areaTotal2", "f8", ("gridcell",))
        v_rlen = ds_par.createVariable("rlen", "f8", ("gridcell",))
        v_rwid = ds_par.createVariable("rwid", "f8", ("gridcell",))
        v_rdep = ds_par.createVariable("rdep", "f8", ("gridcell",))
        v_rslp = ds_par.createVariable("rslp", "f8", ("gridcell",))

        v_cellid[:] = cell_ids
        v_id[:] = ids
        v_dnid[:] = dn_ids
        v_fdir[:] = np.ones(n_cells, dtype=np.int32)
        v_lat[:] = lats
        v_lon[:] = lons
        v_area[:] = np.full(n_cells, 1.0e7)
        v_area2[:] = area_accum
        v_rlen[:] = np.full(n_cells, 5000.0)
        v_rwid[:] = np.full(n_cells, 20.0)
        v_rdep[:] = np.full(n_cells, 2.0)
        v_rslp[:] = np.full(n_cells, 0.001)


def create_synthetic_bifurcation_csv(csv_path: str) -> str:
    """
    Create a synthetic bifurcation/IBT CSV file for testing without requiring pandas.
    """
    header = [
        "cell_lat",
        "cell_long",
        "ibt",
        "note",
        "dncell1_lat",
        "dncell1_long",
        "dncell2_lat",
        "dncell2_long",
        "ratio1",
        "ratio2",
        "demand1",
        "demand2",
    ]
    row = [
        40.0,
        -99.8,
        0,
        "Synthetic bifurcation point",
        40.0,
        -99.7,
        40.1,
        -99.8,
        0.6,
        0.4,
        0.0,
        0.0,
    ]
    with open(csv_path, "w", newline="") as fp:
        writer = csv.writer(fp)
        writer.writerow(header)
        writer.writerow(row)
    return csv_path
