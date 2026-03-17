# p3d2foam

Convert NASA TMR Plot3D multi-block structured grids to OpenFOAM format.

## Background

NASA's [Turbulence Modeling Resource (TMR)](https://tmbwg.github.io/turbmodels/) publishes structured multi-block grids in Plot3D format for CFD validation cases (ONERA M6, flat plate, bump, etc.). These grids come with Neutral Map Files (`.nmf`) that describe block boundaries and connectivity.

This tool converts those grids into OpenFOAM polyMesh format through a pipeline:

```
Plot3D (.p3dfmt / .p3d) + NMF (.nmf)
    → Gmsh MSH 2.2 (.msh)
        → gmshToFoam -keepOrientation
            → stitchMesh (if multi-block)
                → renumberMesh + checkMesh
```

The Gmsh intermediate format is used because `gmshToFoam` handles the complex conversion from hex/quad element connectivity to OpenFOAM's face-based polyMesh format (owner/neighbour/faces/points/boundary). NMF boundary names pass through directly as OpenFOAM patch names.

## Installation

```bash
pip install -e /path/to/p3d2foam
```

This installs the dependencies automatically:
- [plot3d](https://pypi.org/project/plot3d/) (NASA's Plot3D_utilities) -- for `Block` data structure and binary file reading
- [foamlib](https://pypi.org/project/foamlib/) -- for running OpenFOAM utilities and writing OpenFOAM dictionaries
- numpy, pyyaml

## Quick Start

### Generate Gmsh mesh only (no OpenFOAM required)

```bash
p3d2foam grid.p3dfmt --msh-only
```

This reads `grid.p3dfmt` + `grid.nmf` (auto-detected) and writes `grid.msh`.

### Full conversion to OpenFOAM

```bash
# Single-block grid (no stitching needed)
p3d2foam flatplate.p3dfmt -d myCase

# Multi-block grid with a YAML config for interface pairs
p3d2foam m6wing.p3dfmt -c m6wing.yaml -d meshCase

# With scaling (e.g., ONERA M6 wing span normalization)
p3d2foam m6wing.p3dfmt -c m6wing.yaml -d meshCase --scale 1.1963 1.1963 1.1963
```

### Using from Python

```python
from p3d2foam import NeutralMapFile, GmshWriter, read_plot3d

blocks = read_plot3d("m6wing.p3dfmt")
nmf = NeutralMapFile("m6wing.nmf")

writer = GmshWriter()
writer.build(blocks, nmf)
writer.write("m6wing.msh")
```

Or run the full OpenFOAM pipeline:

```python
from p3d2foam.types import ConversionConfig
from p3d2foam.pipeline import FoamPipeline

config = ConversionConfig(
    p3d_file="m6wing.p3dfmt",
    nmf_file="m6wing.nmf",
    case_dir="./meshCase",
    interfaces=[("B1KM", "B4KM"), ("B1J1", "B4J1"), ...],
    scale=(1.1963, 1.1963, 1.1963),
)

pipeline = FoamPipeline(config)
pipeline.run()
```

## CLI Reference

```
p3d2foam [-h] [-m NMF_FILE] [-c CONFIG] [-d CASE_DIR] [-o OUTPUT]
         [--msh-only] [--binary] [--big-endian] [--scale SX SY SZ]
         [-v] p3d_file
```

| Argument | Description |
|----------|-------------|
| `p3d_file` | Plot3D grid file (ASCII or binary) |
| `-m, --nmf-file` | Neutral Map File (default: `<p3d_stem>.nmf`) |
| `-c, --config` | YAML config file for interfaces, scale, etc. |
| `-d, --case-dir` | OpenFOAM case directory (default: `.`) |
| `-o, --output` | Output `.msh` filename (for `--msh-only`) |
| `--msh-only` | Only generate `.msh` file, skip OpenFOAM pipeline |
| `--binary` | Read Plot3D file as binary (default: ASCII) |
| `--big-endian` | Binary file is big-endian |
| `--scale SX SY SZ` | Scale factor for `transformPoints` |
| `-v, --verbose` | Enable debug logging |

## YAML Config File

For multi-block grids where the NMF does not use `ONE-TO-ONE` entries, you need a config file to specify which boundary pairs should be stitched:

```yaml
p3d_file: m6wing.p3dfmt
nmf_file: m6wing.nmf
case_dir: ./mesh

scale: [1.1963, 1.1963, 1.1963]

interfaces:
  - [B1KM, B4KM]
  - [B1J1, B4J1]
  - [B1IM, B2I1]
  - [B2KM, B3KM]
  - [B2IM, B3I1]
  - [B3IM, B4I1]
```

Interface pairs are the NMF boundary names that share coincident nodes between blocks. These get passed to `stitchMesh -perfect` to merge them.

If the NMF file uses `ONE-TO-ONE` or `Patched` entries, interfaces are detected automatically and no config file is needed.

## Pipeline Details

The full OpenFOAM pipeline performs these steps:

1. **Read inputs** -- Parse the Plot3D grid and NMF boundary file
2. **Validate** -- Check that block counts and dimensions match between the two files
3. **Bootstrap case** -- Create minimal `system/controlDict`, `fvSchemes`, `fvSolution` (required by OpenFOAM utilities)
4. **Generate Gmsh mesh** -- Write `.msh` file with hex cells and named quad boundary faces
5. **gmshToFoam** -- Convert to OpenFOAM polyMesh (`-keepOrientation` flag is critical)
6. **stitchMesh** -- Merge multi-block interfaces (skipped for single-block grids)
7. **transformPoints** -- Apply scaling if specified
8. **renumberMesh** -- Optimize cell ordering for solver performance
9. **checkMesh** -- Validate the final mesh

NMF boundary names become OpenFOAM patch names directly. All patches are created as type `patch`. To reclassify them (e.g., set walls, symmetry planes), use OpenFOAM's `changeDictionary` or `createPatch` utilities afterward.

### Interface Resolution

For multi-block grids, stitch pairs are resolved in this order:

1. **YAML config** -- Explicit `interfaces` list (highest priority)
2. **NMF file** -- `ONE-TO-ONE` or `Patched` entries parsed from the NMF
3. **Geometric detection** -- `plot3d.connectivity_fast()` auto-detects coincident block faces (best-effort fallback)

### The -keepOrientation Flag

The `gmshToFoam -keepOrientation` flag is essential. Without it, `gmshToFoam` attempts to "fix" hex orientation by inverting cells it deems incorrectly ordered. For Plot3D-derived meshes, this check incorrectly selects valid cells for reordering, producing negative cell volumes and wrong-oriented faces. The Plot3D vertex ordering is already correct -- `keepOrientation` tells `gmshToFoam` to trust it.

## NMF Format Reference

The Neutral Map File ([official spec](https://tmbwg.github.io/turbmodels/nmf_documentation.html)) describes block boundaries and connectivity for multi-block Plot3D grids.

### Face ID Convention

| Face | Constant Index | Primary (S1-E1) | Secondary (S2-E2) |
|------|----------------|------------------|--------------------|
| 1    | kmin           | i                | j                  |
| 2    | kmax           | i                | j                  |
| 3    | imin           | j                | k                  |
| 4    | imax           | j                | k                  |
| 5    | jmin           | k                | i                  |
| 6    | jmax           | k                | i                  |

### Supported Boundary Types

- **Named boundaries** (quoted): `'wall'`, `'symm1'`, `'farfield'`, etc. -- single-sided, 7 fields
- **Standard keywords** (unquoted): `WALL`, `Symmetry-X`, `Symmetry-Y`, `Symmetry-Z`, `Inflow`, `Outflow` -- single-sided, 7 fields
- **ONE-TO-ONE / ONE_TO_ONE**: Point-matched block interface -- two-sided, 13-14 fields
- **Patched**: Non-conforming block interface -- two-sided, 13-14 fields

## Package Structure

```
p3d2foam/
    pyproject.toml
    README.md
    p3d2foam/
        __init__.py          # Public API exports
        __main__.py          # CLI entry point (python -m p3d2foam)
        types.py             # Boundary, InterfacePair, ConversionConfig
        nmf.py               # NeutralMapFile parser
        p3d_reader.py        # ASCII Plot3D reader with repeat-notation expansion
        gmsh_writer.py       # Gmsh MSH 2.2 writer (hex + quad elements)
        pipeline.py          # FoamPipeline: full conversion via foamlib
        config.py            # YAML config loader
```

## Provenance

The Gmsh writer logic (hex vertex ordering, boundary quad winding, node ID scheme) is ported from [p3d2gmsh.py](https://github.com/furstj/myFoam) by Alexey Matveichev (MIT License, 2015-2016), originally bundled with the HiSA solver's ONERA M6 example. The NMF parser is a rewrite based on the [official NMF specification](https://tmbwg.github.io/turbmodels/nmf_documentation.html). The Plot3D reader handles both standard ASCII format and Fortran repeat notation (e.g., `3578*0.0`).
