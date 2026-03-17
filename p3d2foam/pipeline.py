"""OpenFOAM mesh conversion pipeline using foamlib.

Orchestrates the full Plot3D -> Gmsh -> OpenFOAM conversion:
  1. Read Plot3D grid (via plot3d library)
  2. Read NMF boundaries
  3. Generate Gmsh .msh
  4. gmshToFoam -keepOrientation
  5. stitchMesh (if multi-block)
  6. transformPoints (if scale specified)
  7. renumberMesh + checkMesh
"""

from __future__ import annotations

import logging
from pathlib import Path

from foamlib import FoamCase, FoamFile

from .gmsh_writer import GmshWriter
from .nmf import NeutralMapFile
from .types import ConversionConfig

logger = logging.getLogger(__name__)


class FoamPipeline:
    """Convert a Plot3D grid to an OpenFOAM polyMesh via Gmsh intermediate."""

    def __init__(self, config: ConversionConfig) -> None:
        self.config = config

    def run(self) -> None:
        """Execute the full conversion pipeline."""
        from .p3d_reader import read_plot3d

        cfg = self.config
        case_dir = Path(cfg.case_dir).resolve()

        # Step 1: Read inputs
        logger.info("Reading Plot3D file: %s", cfg.p3d_file)
        blocks = read_plot3d(cfg.p3d_file, binary=cfg.binary, big_endian=cfg.big_endian)
        logger.info("Read %d blocks", len(blocks))

        logger.info("Reading NMF file: %s", cfg.nmf_file)
        nmf = NeutralMapFile(cfg.nmf_file)
        logger.info("%s", nmf)

        # Step 2: Validate
        if nmf.nblocks != len(blocks):
            raise ValueError(
                f"NMF declares {nmf.nblocks} blocks but Plot3D file has {len(blocks)}"
            )
        for i, (blk, dims) in enumerate(zip(blocks, nmf.block_dims)):
            blk_shape = blk.X.shape
            if blk_shape != dims:
                raise ValueError(
                    f"Block {i+1} dims mismatch: Plot3D {blk_shape} vs NMF {dims}"
                )

        # Step 3: Bootstrap case directory
        self._bootstrap_case(case_dir)
        case = FoamCase(case_dir)

        # Step 4: Generate Gmsh .msh
        msh_path = case_dir / "mesh.msh"
        writer = GmshWriter()
        writer.build(blocks, nmf)
        writer.write(msh_path)

        # Step 5: gmshToFoam
        logger.info("Running gmshToFoam")
        case.run(["gmshToFoam", "-keepOrientation", str(msh_path)])

        # Step 6: stitchMesh (if multi-block)
        if len(blocks) > 1:
            pairs = self._resolve_interfaces(nmf, blocks)
            if pairs:
                self._write_stitch_tolerances(case_dir)
                for name_a, name_b in pairs:
                    logger.info("Stitching %s <-> %s", name_a, name_b)
                    case.run([
                        "stitchMesh", "-perfect",
                        "-toleranceDict", "stitchMeshToleranceDict",
                        name_a, name_b, "-overwrite",
                    ])
                    # Clean up meshPhi if created
                    mesh_phi = case_dir / "0" / "meshPhi"
                    if mesh_phi.exists():
                        mesh_phi.unlink()

        # Step 7: transformPoints (if scale specified)
        if cfg.scale:
            sx, sy, sz = cfg.scale
            logger.info("Scaling mesh by (%s, %s, %s)", sx, sy, sz)
            case.run(["transformPoints", "-scale", f"({sx} {sy} {sz})"])

        # Step 8: renumberMesh + checkMesh
        logger.info("Running renumberMesh")
        case.run(["renumberMesh", "-overwrite"])

        logger.info("Running checkMesh")
        case.run(["checkMesh"])

        logger.info("Conversion complete: %s", case_dir / "constant" / "polyMesh")

    def _resolve_interfaces(
        self, nmf: NeutralMapFile, blocks: list
    ) -> list[tuple[str, str]]:
        """Determine stitch pairs from config, NMF, or geometric detection."""
        # Priority 1: Explicit config
        if self.config.interfaces:
            logger.info("Using interface pairs from config")
            return list(self.config.interfaces)

        # Priority 2: NMF ONE-TO-ONE / Patched entries
        if nmf.interfaces:
            logger.info("Using interface pairs from NMF")
            return [
                (ip.boundary_a.name, ip.boundary_b.name)
                for ip in nmf.interfaces
            ]

        # Priority 3: Geometric auto-detection via Plot3D_utilities
        logger.info("Auto-detecting interfaces via connectivity_fast()")
        try:
            from plot3d import connectivity_fast
            face_matches, _ = connectivity_fast(blocks)
            # face_matches contains dicts with block indices — map to NMF names
            # This is a best-effort fallback; for complex cases, use config
            pairs = self._face_matches_to_pairs(face_matches, nmf)
            if pairs:
                return pairs
        except ImportError:
            logger.warning("plot3d.connectivity_fast not available")
        except Exception:
            logger.warning("Geometric interface detection failed", exc_info=True)

        logger.warning("No interface pairs found — skipping stitchMesh")
        return []

    @staticmethod
    def _face_matches_to_pairs(
        face_matches: list[dict], nmf: NeutralMapFile
    ) -> list[tuple[str, str]]:
        """Convert plot3d connectivity face_matches to NMF boundary name pairs.

        This is approximate — it matches by block index and face type.
        For reliable results, use explicit config or NMF ONE-TO-ONE entries.
        """
        # Build lookup from (block_0based, face_type) to boundary name
        # face_type in plot3d: 'imin','imax','jmin','jmax','kmin','kmax'
        face_type_to_id = {
            "kmin": 1, "kmax": 2,
            "imin": 3, "imax": 4,
            "jmin": 5, "jmax": 6,
        }
        name_lookup: dict[tuple[int, int], str] = {}
        for bdry in nmf.boundaries:
            name_lookup[(bdry.block - 1, bdry.face_id)] = bdry.name

        pairs: list[tuple[str, str]] = []
        for match in face_matches:
            b1 = match.get("block1", match.get("block_index_1"))
            b2 = match.get("block2", match.get("block_index_2"))
            f1_type = match.get("face1", match.get("face_type_1", ""))
            f2_type = match.get("face2", match.get("face_type_2", ""))

            fid1 = face_type_to_id.get(f1_type)
            fid2 = face_type_to_id.get(f2_type)

            if fid1 and fid2:
                name1 = name_lookup.get((b1, fid1))
                name2 = name_lookup.get((b2, fid2))
                if name1 and name2:
                    pairs.append((name1, name2))

        return pairs

    def _write_stitch_tolerances(self, case_dir: Path) -> None:
        """Write constant/stitchMeshToleranceDict."""
        tol_path = case_dir / "constant" / "stitchMeshToleranceDict"

        # Use config tolerances or sensible defaults
        tols = self.config.stitch_tolerances or {
            "pointMergeTol": 0.3,
            "edgeMergeTol": 0.05,
            "integralAdjTol": 0.15,
            "edgeMasterCatchFraction": 0.4,
            "edgeEndCutoffTol": 0.0001,
            "edgeCoPlanarTol": 0.8,
        }

        f = FoamFile(tol_path)
        for key, val in tols.items():
            f[key] = val

    @staticmethod
    def _bootstrap_case(case_dir: Path) -> None:
        """Ensure minimal OpenFOAM case structure exists.

        OpenFOAM utilities require system/controlDict, system/fvSchemes,
        and system/fvSolution to exist even for mesh-only operations.
        """
        (case_dir / "system").mkdir(parents=True, exist_ok=True)
        (case_dir / "constant").mkdir(exist_ok=True)

        # Minimal controlDict
        cd_path = case_dir / "system" / "controlDict"
        if not cd_path.exists():
            cd = FoamFile(cd_path)
            cd["application"] = "simpleFoam"
            cd["startFrom"] = "startTime"
            cd["startTime"] = 0
            cd["endTime"] = 1
            cd["deltaT"] = 1
            cd["writeControl"] = "timeStep"
            cd["writeInterval"] = 1
            cd["writeFormat"] = "ascii"
            cd["writePrecision"] = 6
            cd["writeCompression"] = "off"
            cd["timeFormat"] = "general"
            cd["timePrecision"] = 6
            cd["runTimeModifiable"] = True

        # Minimal fvSchemes (empty is sufficient for mesh utilities)
        fvs_path = case_dir / "system" / "fvSchemes"
        if not fvs_path.exists():
            fvs = FoamFile(fvs_path)
            fvs["ddtSchemes"] = {}
            fvs["gradSchemes"] = {}
            fvs["divSchemes"] = {}
            fvs["laplacianSchemes"] = {}
            fvs["interpolationSchemes"] = {}
            fvs["snGradSchemes"] = {}

        # Minimal fvSolution
        fvsol_path = case_dir / "system" / "fvSolution"
        if not fvsol_path.exists():
            fvsol = FoamFile(fvsol_path)
            fvsol["solvers"] = {}
