# HexWatershed Utility Package: Codebase Review and Improvement Plan

**Date:** October 2026  
**Package:** `hexwatershed_utility` (v0.1.0)  
**Author / Maintainer:** Chang Liao  
**Review Target:** Core library (`hexwatershed_utility/`), workflow scripts (`examples/`), build & testing configurations  

---

## 1. Executive Summary & Codebase Health Assessment

`hexwatershed_utility` is a specialized scientific Python package designed to bridge geospatial workflows, hydrologic preprocessing (watershed boundaries, river networks, coastlines, dams, lakes), computational meshes (MPAS polygonal meshes), and Earth System Models (ESMs, particularly E3SM MOSART and ELM).

While the library implements sophisticated hydrologic and mesh algorithms, the codebase exhibits critical architectural, packaging, and performance challenges resulting from rapid research development. 

### Summary Scorecard

| Area | Current Rating | Core Challenge |
| :--- | :--- | :--- |
| **Packaging & Dependencies** | ⚠️ Critical | `pyproject.toml` omits essential dependencies (GDAL, pandas, rtree, scipy, pyearth, pyesm, etc.). |
| **Architectural Separation** | ⚠️ Critical | Reusable core domain logic is stranded in `examples/`; unparameterized one-off HPC scripts reside in `hexwatershed_utility/`. |
| **Algorithmic Performance** | ⚠️ Moderate–High | Multiple $O(N^2)$ linear search bottlenecks in cell graph traversals; deep recursion requiring `sys.setrecursionlimit(100000)`. |
| **Resource & Error Handling** | ⚠️ Moderate | Unclosed NetCDF file descriptors; module-level side effects (GDAL error handler overrides); silent failures via `print()` + `return None`. |
| **Code Consistency & Typing** | ⚠️ Moderate | Variable naming combines Hungarian notation, camelCase, and snake_case; limited type hinting; dead `pass` statements. |
| **Testing & CI/CD** | ⚠️ Critical | Only 4 test files covering ~10% of functionality; no mock fixtures; zero GitHub Actions / automated CI pipelines. |
| **Documentation** | ⚠️ Critical | `README.md` is minimal (5 lines); `.readthedocs.yaml` points to non-existent `docs/source/conf.py`. |

---

## 2. In-Depth Codebase Audit & Key Findings

### 2.1 Packaging & Dependency Management
1. **Undeclared Dependencies in `pyproject.toml`:**
   - Declared dependencies: `numpy>=1.20.0`, `netCDF4>=1.5.0`, `shapely>=1.8.0`, `requests`.
   - **Actually imported across library modules:**
     - `osgeo` (`gdal`, `ogr`, `osr`)
     - `pandas` (used in `create_bifurc_netcdf.py`)
     - `rtree` (used in `simplify_hydrorivers_network.py` and mesh searches)
     - `scipy` (used for `KDTree` in auxiliary utilities)
     - Proprietary/ecosystem packages: `pyearth`, `pyearthbuffer`, `pyearthviz`, `pyearthviz3d`, `pyesm`, `pyflowline`, `geovista`, `tinyr`.
   - *Impact:* `pip install hexwatershed_utility` fails immediately on import due to `ModuleNotFoundError`.
2. **Environment Inconsistency:**
   - `environment.yml` defines some packages (`gdal`, `rtree`, `cartopy`, `matplotlib`), but omits `pyearth`, `pyesm`, `pyflowline`, and `pandas`.
3. **Pre-commit vs Linter Configuration Mismatch:**
   - `pyproject.toml` sets Black `line-length = 100`.
   - `.pre-commit-config.yaml` configures Black and isort with `--line-length=88`. This triggers continuous formatting churn on every commit.

### 2.2 Architectural Inversion Between `hexwatershed_utility/` and `examples/`
There is a major role reversal between the library package and the examples directory:
- **Core Domain Logic Stranded in `examples/`:**
  - `examples/auxiliary/find_cells_on_vector_flowline.py`: Spatial mapping between vector lines and mesh cells.
  - `examples/auxiliary/find_downstream_cells.py`: Downstream cell graph traversal.
  - `examples/auxiliary/meshjson_dnID.py`: Conversion of HexWatershed JSON connectivity to downstream index arrays.
  - `examples/auxiliary/make_dam_dependency.py`: Upstream/downstream dam topological dependency analysis.
  - `examples/postprocess/river/calculate_nse.py`: Reusable hydrological validation metric (Nash–Sutcliffe Efficiency).
  - `examples/postprocess/river/find_gage_mesh_cell_id.py` (261 lines): Production-grade R-tree spatial snapping for streamflow gauges to mesh cells.
  - `examples/general/convert_hexwatershed_output_to_raster.py` (266 lines): Flow direction raster re-gridding.
  - `examples/postprocess/mpas/netcdf/inject_dem_into_mpas_mesh_netcdf.py` (381 lines): DEM elevation interpolation into MPAS unstructured mesh.
- **Procedural Scripts Placed Inside the Library:**
  - `hexwatershed_utility/preprocess/feature/city/create_city_point.py`: Bare executable script reading `/compyfs/liao313/00raw/city/...` at import time with no functions defined.
  - `hexwatershed_utility/preprocess/feature/city/create_city_mask.py`: Bare executable script.
  - `hexwatershed_utility/preprocess/feature/lake/extract_lake_polygon.py` and `create_lake_polygon.py`: Contain no function definitions, but their subpackage `__init__.py` attempts to import `from ... import extract_lake_polygon`, masking the failure behind `try ... except ImportError: pass`.
  - `hexwatershed_utility/preprocess/feature/boundary/download_conus_boundary.py`: Executes an ArcGIS REST HTTP download on import.
  - `hexwatershed_utility/postprocess/river/create_bifurc_netcdf.m`: A 487-line MATLAB script placed in the Python source directory.

### 2.3 Algorithmic Performance & Scaling Bottlenecks
1. **$O(N^2)$ Graph Traversal via Linear Scanning:**
   - In `hexwatershed_utility/postprocess/river/convert_hexwatershed_outputs_to_river_inputs.py`:
     ```python
     for i in range(nCell):
         lCellID_downslope = aCellID_downslope[i]
         dummy_index = np.where(aCellID == lCellID_downslope)
         dnID = aID[dummy_index[0]]
         aDnID.append(dnID[0])
     ```
     For $N = 1.5 \times 10^6$ cells, this runs $1.5 \times 10^6 \times 1.5 \times 10^6 \approx 2.25 \times 10^{12}$ operations. A hash map `dict(zip(aCellID, aID))` reduces this to $O(N)$ (fractions of a second).
   - Identical linear scanning patterns exist in `examples/auxiliary/meshjson_dnID.py` and `examples/auxiliary/find_downstream_cells.py`.
2. **Monolithic Functions & Excessive Recursion:**
   - `hexwatershed_utility/preprocess/feature/river/simplify_hydrorivers_network.py` is 1,428 lines long.
   - It sets `sys.setrecursionlimit(100000)` at module import time due to deep recursive calls in `tag_upstream()` and `is_downstream()`. This risks C-level stack overflow crashes on large continental river networks.
3. **Memory Optimization Status:**
   - The memory patch in `get_geometry.py` (`MEMORY_FIX_README.md`) resolved the 264 GiB allocation error using chunking (`nChunk_size`). However, similar streaming/chunking patterns must be audited across the remaining NetCDF and raster utilities.

### 2.4 Resource Safety, Modularity & Error Handling
1. **Unclosed File Descriptors:**
   - Multiple functions instantiate NetCDF datasets (`pDatasets_out = nc.Dataset(..., "w")`) without using `with` context managers or calling `.close()`. This leaves files susceptible to truncation or filesystem locks.
2. **Global Side-Effects on Import:**
   - `convert_mpas_mesh_to_vector.py` invokes `gdal.UseExceptions()` and `gdal.PushErrorHandler("CPLQuietErrorHandler")` at the top level, silently mutating global GDAL configuration for any application importing the package.
3. **Error Reporting:**
   - Functions often catch errors or check preconditions using:
     ```python
     if not os.path.exists(mesh_file):
         print("The mesh file does not exist, please check the file path")
         return
     ```
     This should raise `FileNotFoundError` or custom domain exceptions (`MeshValidationError`), allowing callers to handle or recover from failure states cleanly.

### 2.5 Testing & CI/CD
- **Current Tests:** 4 unit tests covering `TestPackage`, `reorganize_boundary`, `find_contributing_cells`, and empirical `get_geometry`.
- **Untested Systems:**
  - River and land domain NetCDF generation
  - Bifurcation NetCDF creation
  - MPAS to vector conversion
  - Dam snapping and HydroRivers simplification
  - Gauge-to-mesh matching algorithm
- **No CI Workflow:** There is no `.github/workflows/` directory. Commits and pull requests are not automatically validated against linting, type checks, or test suites.

---

## 3. Target Architecture & Package Layout

```
hexwatershed_utility/
├── __init__.py                     # Package metadata and public API exports
├── py.typed
├── core/                           # Shared domain models and foundational algorithms
│   ├── __init__.py
│   ├── graph.py                    # O(N) downstream/upstream traversal, cell connectivity
│   ├── geometry.py                 # Spherical coordinates, polygon areas, IDL wrapping
│   ├── spatial_index.py            # R-Tree and KDTree helper wrappers
│   └── exceptions.py               # Custom package exceptions
├── io/                             # Robust, managed file I/O
│   ├── __init__.py
│   ├── netcdf.py                   # Context-managed NetCDF read/write helpers
│   ├── geojson.py                  # Standardized GeoJSON reading/writing
│   └── hex_json.py                 # HexWatershed JSON reader/validator
├── preprocess/                     # Pre-processing modules (Feature & Mesh)
│   ├── __init__.py
│   ├── boundary/                   # Watershed and administrative boundary operations
│   │   ├── __init__.py
│   │   ├── reorganize.py
│   │   └── download.py             # Parameterized download functions (not run on import)
│   ├── coastline/                  # Land-ocean masks and HydroBASIN/NaturalEarth fixes
│   │   ├── __init__.py
│   │   ├── land_ocean_mask.py
│   │   └── naturalearth_fix.py
│   ├── river/                      # River network simplification and tagging
│   │   ├── __init__.py
│   │   ├── simplify.py             # Iterative (non-recursive) river simplifier
│   │   └── outlet.py
│   ├── dam/                        # Dam snapping and dependency graphs
│   │   ├── __init__.py
│   │   └── snap.py
│   ├── lake/                       # Parameterized lake polygon and mask generation
│   │   ├── __init__.py
│   │   └── lake_features.py
│   └── mesh/                       # MPAS and computational mesh pre-processing
│       ├── __init__.py
│       ├── mpas_to_vector.py
│       └── dem_injection.py        # Promoted from examples/
├── postprocess/                    # Model post-processing modules
│   ├── __init__.py
│   ├── river/                      # MOSART river parameters, geometry & bifurcations
│   │   ├── __init__.py
│   │   ├── parameter_generator.py
│   │   ├── geometry.py             # Chunked hydraulic geometry
│   │   ├── bifurcation.py          # Bifurcation NetCDF generation
│   │   └── gauge_snapping.py       # Promoted from examples/
│   ├── land/                       # ELM / land domain and SCRIP grid generation
│   │   ├── __init__.py
│   │   └── domain_generator.py
│   └── metrics/                    # Hydrologic validation metrics
│       ├── __init__.py
│       └── efficiency.py           # NSE, KGE, RMSE (promoted from examples/)
└── visual/                         # Visualization routines (optional [viz] dependency)
    ├── __init__.py
    ├── flow_direction.py
    └── mesh_visualization.py
```

---

## 4. Phased Improvement Plan

### Phase 1: Dependency Harmonization, Packaging & CI/CD Pipeline
*Objective: Ensure the package can be installed, formatted, and tested reproducibly on clean environments.*

1. **Synchronize Dependencies (`pyproject.toml`):**
   - Core dependencies:
     ```toml
     dependencies = [
         "numpy>=1.22.0",
         "netCDF4>=1.5.8",
         "shapely>=2.0.0",
         "pandas>=1.5.0",
         "scipy>=1.9.0",
         "rtree>=1.0.0",
         "requests>=2.28.0",
     ]
     ```
   - Optional dependency groups:
     ```toml
     [project.optional-dependencies]
     geo = ["gdal>=3.4.0"]          # Documented conda/system requirements
     viz = ["matplotlib>=3.5.0", "cartopy>=0.21.0", "pyearthviz>=0.1.1"]
     viz3d = ["pyearthviz3d>=0.1.1", "geovista>=0.4.0"]
     dev = ["pytest>=7.4", "pytest-cov>=4.1", "black==24.2.0", "flake8>=6.1", "mypy>=1.5", "isort>=5.12"]
     docs = ["sphinx>=7.0", "sphinx-rtd-theme>=2.0", "myst-parser>=2.0"]
     all = ["hexwatershed_utility[geo,viz,viz3d,dev,docs]"]
     ```
2. **Align Tool Formatting Configurations:**
   - Standardize line length to **88** or **100** consistently across `pyproject.toml`, `.flake8`, and `.pre-commit-config.yaml`.
3. **Establish GitHub Actions CI (`.github/workflows/ci.yml`):**
   - Matrix builds testing Python 3.10, 3.11, 3.12 across Ubuntu and macOS.
   - Run `pre-commit` hooks (Black, Flake8, isort).
   - Install Conda environment with GDAL via `conda-incubator/setup-miniconda`.
   - Run `pytest --cov=hexwatershed_utility`.

### Phase 2: Architectural Reorganization & Logic Relocation
*Objective: Clean separation of library code vs. user examples.*

1. **Promote Reusable Algorithms from `examples/` into `hexwatershed_utility/`:**
   - Move `examples/postprocess/river/find_gage_mesh_cell_id.py` → `hexwatershed_utility/postprocess/river/gauge_snapping.py`.
   - Move `examples/postprocess/mpas/netcdf/inject_dem_into_mpas_mesh_netcdf.py` → `hexwatershed_utility/preprocess/mesh/dem_injection.py`.
   - Move `examples/postprocess/river/calculate_nse.py` → `hexwatershed_utility/postprocess/metrics/efficiency.py` (add KGE and RMSE).
   - Move `examples/auxiliary/find_cells_on_vector_flowline.py` and `find_downstream_cells.py` → `hexwatershed_utility/core/graph.py`.
2. **Transform Library Scripts into Parameterized Callable Functions:**
   - Refactor `preprocess/feature/city/create_city_point.py` and `create_city_mask.py` into functions taking input/output paths and parameters (e.g., `def create_city_points_from_csv(csv_path, output_geojson, min_population=1e6): ...`).
   - Refactor `preprocess/feature/lake/*.py` into clean library functions.
   - Remove bare execution on module import (wrap script runners with `if __name__ == "__main__":`).
3. **Archive / Relocate Non-Python Assets:**
   - Move `hexwatershed_utility/postprocess/river/create_bifurc_netcdf.m` out of package sources to an archive directory (e.g. `legacy/matlab/`).

### Phase 3: Algorithmic Optimization & Robust Resource Management
*Objective: Eliminate scaling bottlenecks and ensure fault-tolerant I/O.*

1. **Convert $O(N^2)$ Cell Scans to $O(N)$ Hash-Lookups:**
   - In `convert_hexwatershed_outputs_to_river_inputs.py`:
     ```python
     # Replace:
     # dummy_index = np.where(aCellID == lCellID_downslope)
     # With:
     cell_id_to_idx = {cid: idx for idx, cid in enumerate(aCellID)}
     downstream_indices = [
         cell_id_to_idx[dn_id] if dn_id in cell_id_to_idx else -9999
         for dn_id in aCellID_downslope
     ]
     ```
2. **Refactor River Network Simplification (`simplify_hydrorivers_network.py`):**
   - Replace recursive upstream traversal with iterative stack/queue traversal (e.g., Breadth-First or Depth-First Search with an explicit Python `collections.deque`).
   - Remove `sys.setrecursionlimit(100000)`.
   - Decompose the 1,428-line function into focused helper modules:
     - `_segment_indexer.py`
     - `_rtree_manager.py`
     - `_geometry_cache.py`
3. **Resource Safety via Context Managers:**
   - Replace raw `nc.Dataset(path, 'w')` calls with `with nc.Dataset(path, 'w') as ds:` to ensure datasets are reliably flushed and closed on exceptions.
   - Enforce localized GDAL error handling without altering global state.

### Phase 4: API Modernization, Error Handling & Type Annotations
*Objective: Deliver an idiomatic, professional developer experience.*

1. **Structured Exception Hierarchy:**
   - Create `hexwatershed_utility.core.exceptions`:
     - `HexWatershedError` (base)
     - `MeshFormatError`
     - `TopologicalLoopError`
     - `DataInconsistencyError`
2. **Replace `print()` with Standard Python `logging`:**
   - Configure a package logger `logger = logging.getLogger(__name__)`.
   - Provide clean log levels (`debug`, `info`, `warning`, `error`) so users can silence or redirect output.
3. **Add Complete Type Annotations:**
   - Annotate all public functions with typed signatures (`PathLike`, `NDArray`, `Optional[Dict]`).
   - Enforce `mypy` type checking in CI.

### Phase 5: Comprehensive Testing Suite & Fixtures
*Objective: Elevate test coverage to >85% without requiring external HPC cluster datasets.*

1. **Build Lightweight Test Fixtures (`tests/fixtures/`):**
   - Provide minimal synthetic data fixtures:
     - Synthetic 7-cell hexagonal mesh GeoJSON
     - Synthetic 20-reach river network shapefile / GeoJSON
     - Minimal mock MOSART NetCDF file generated via `pytest.fixture`
2. **Expand Test Suites:**
   - `test_graph.py`: Upstream/downstream tracing, cycle detection, outlet identification.
   - `test_parameter_conversion.py`: HexWatershed JSON to MOSART NetCDF conversion with synthetic data.
   - `test_bifurcation.py`: Structured and unstructured bifurcation NetCDF generation.
   - `test_gauge_snapping.py`: Snapping precision against synthetic stream gauges.
   - `test_metrics.py`: Edge cases for NSE/KGE calculations (all NaN, constant series, perfect correlation).

### Phase 6: Documentation & Example Overhaul
*Objective: Enable new users and researchers to learn and apply the tools easily.*

1. **Complete ReadTheDocs Setup:**
   - Create `docs/source/conf.py`, `docs/source/index.rst`, and `docs/requirements.txt`.
   - Configure Sphinx with MyST parser to render markdown guides.
   - Auto-generate API documentation using `sphinx.ext.autodoc` and `sphinx.ext.napoleon`.
2. **Revamp `examples/` as Reproducible Tutorials:**
   - Rewrite examples to import from `hexwatershed_utility`.
   - Include self-contained sample data in `data/sample/` so examples can run locally out-of-the-box without requiring access to `/compyfs/`.
   - Provide a Quickstart Jupyter Notebook demonstrating the end-to-end preprocessing to MOSART input workflow.

---

## 5. Implementation Roadmap & Priority Matrix

| Task | Priority | Effort | Target Milestone | Impact |
| :--- | :--- | :--- | :--- | :--- |
| Fix `pyproject.toml` dependencies & pre-commit config | **P0 (Immediate)** | Low (0.5 day) | v0.1.1 | Makes package cleanly installable |
| Replace $O(N^2)$ linear scans with hash lookups | **P0 (Immediate)** | Low (0.5 day) | v0.1.1 | 100x–1000x speedup on large meshes |
| Wrap NetCDF/file operations in context managers | **P0 (Immediate)** | Low (0.5 day) | v0.1.1 | Prevents file corruption and leaks |
| Move MATLAB script (`.m`) to legacy archive | **P1 (High)** | Minimal (1 hr) | v0.1.1 | Cleans up package source tree |
| Setup GitHub Actions CI workflow | **P1 (High)** | Medium (1 day) | v0.1.2 | Automated testing on every PR |
| Promote gauge snapping & DEM injection from `examples/` | **P1 (High)** | Medium (2 days) | v0.2.0 | Consolidates core algorithms |
| Convert script-like modules in `preprocess/` to functions | **P1 (High)** | Medium (2 days) | v0.2.0 | Eliminates hardcoded compyfs paths |
| Deconstruct `simplify_hydrorivers_network.py` & fix recursion | **P2 (Medium)** | High (3–4 days) | v0.2.0 | Prevents recursion crashes on large basins |
| Build synthetic test fixtures and expand test suite | **P2 (Medium)** | Medium (3 days) | v0.2.1 | Guarantees test stability across platforms |
| Create Sphinx / ReadTheDocs documentation suite | **P2 (Medium)** | Medium (2 days) | v0.3.0 | Public documentation & user guides |

---

## 6. Conclusion

By addressing these core structural and algorithmic items, `hexwatershed_utility` will transition from a collection of interconnected research scripts into a production-grade, highly scalable scientific modeling library ready for broader community adoption within E3SM and the hydrologic modeling community.
