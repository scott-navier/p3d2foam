"""Write Gmsh MSH 2.2 format from Plot3D blocks and NMF boundaries.

The vertex ordering for hexahedra and boundary quads is preserved exactly
from the original p3d2gmsh.py (Alexey Matveichev, MIT License, 2015-2016).
This ordering is critical for ``gmshToFoam -keepOrientation`` to produce
valid OpenFOAM meshes without inverted cells.
"""

from __future__ import annotations

import logging
from pathlib import Path

from .nmf import NeutralMapFile
from .types import Boundary

logger = logging.getLogger(__name__)

# Try importing plot3d Block type for type hints.
# At runtime we just duck-type on .X, .Y, .Z attributes.
try:
    from plot3d import Block
except ImportError:
    Block = object  # type: ignore[misc,assignment]

# Hex vertex shifts relative to grid point (i, j, k).
# Each hex cell at (i, j, k) uses the 8 nodes from (i-1, j-1, k-1) to (i, j, k).
# Order matches Gmsh hexahedron node numbering (type 5).
_HEX_SHIFTS = [
    (-1, -1, -1),
    (-1, 0, -1),
    (-1, 0, 0),
    (-1, -1, 0),
    (0, -1, -1),
    (0, 0, -1),
    (0, 0, 0),
    (0, -1, 0),
]


class GmshWriter:
    """Build and write a Gmsh MSH 2.2 file from Plot3D blocks and NMF boundaries."""

    def __init__(self) -> None:
        self._nodes: list[tuple[int, float, float, float]] = []
        self._elements: list[list[int]] = []
        self._groups: list[tuple[int, int, str]] = []
        self._element_id: int = 0

    def build(self, blocks: list[Block], nmf: NeutralMapFile) -> None:
        """Consume Plot3D blocks and NMF to populate nodes, elements, groups."""
        # Physical group 1 is always the 3D volume ("internal")
        self._groups.append((3, 1, "internal"))

        # Add nodes and hex cells for each block
        for blk_idx in range(len(blocks)):
            logger.info(
                "Processing block %d/%d (%s)",
                blk_idx + 1,
                len(blocks),
                blocks[blk_idx].X.shape,
            )
            self._add_block_nodes(blocks, blk_idx)
            self._add_block_hexes(blocks, blk_idx)

        # Add boundary quad faces from NMF
        for bdry in nmf.boundaries:
            self._add_boundary_quads(blocks, bdry)

        # Also add interface faces (so they become patches for stitchMesh)
        for iface in nmf.interfaces:
            self._add_boundary_quads(blocks, iface.boundary_a)
            self._add_boundary_quads(blocks, iface.boundary_b)

        logger.info(
            "Built mesh: %d nodes, %d elements, %d groups",
            len(self._nodes),
            len(self._elements),
            len(self._groups),
        )

    def write(self, filename: str | Path) -> None:
        """Write the mesh to a Gmsh MSH 2.2 ASCII file."""
        path = Path(filename)
        logger.info("Writing %s", path)

        with open(path, "w") as fp:
            # Header
            fp.write("$MeshFormat\n2.2 0 8\n$EndMeshFormat\n")

            # Physical names
            fp.write("$PhysicalNames\n")
            fp.write(f"{len(self._groups)}\n")
            for dim, gid, name in self._groups:
                fp.write(f'{dim} {gid} "{name}"\n')
            fp.write("$EndPhysicalNames\n")

            # Nodes
            fp.write("$Nodes\n")
            fp.write(f"{len(self._nodes)}\n")
            for nid, x, y, z in self._nodes:
                fp.write(f"{nid} {x:15.13e} {y:15.13e} {z:15.13e}\n")
            fp.write("$EndNodes\n")

            # Elements
            fp.write("$Elements\n")
            fp.write(f"{len(self._elements)}\n")
            for el in self._elements:
                fp.write(" ".join(str(v) for v in el) + "\n")
            fp.write("$EndElements\n")

    # ---- Node ID calculation ----

    @staticmethod
    def _node_id(
        blocks: list[Block], block_idx: int, i: int, j: int, k: int
    ) -> int:
        """Compute global 1-based node ID from block-local (i, j, k).

        Nodes are numbered sequentially across blocks. Within each block,
        the ordering is column-major: k varies fastest, then j, then i.
        """
        base = 1
        for idx in range(block_idx):
            di, dj, dk = blocks[idx].X.shape
            base += di * dj * dk

        _, dj, dk = blocks[block_idx].X.shape
        return base + k + dk * j + dk * dj * i

    # ---- Block processing ----

    def _add_block_nodes(self, blocks: list[Block], blk_idx: int) -> None:
        """Add all grid points of a block as Gmsh nodes."""
        blk = blocks[blk_idx]
        idim, jdim, kdim = blk.X.shape

        for i in range(idim):
            for j in range(jdim):
                for k in range(kdim):
                    nid = self._node_id(blocks, blk_idx, i, j, k)
                    self._nodes.append(
                        (nid, float(blk.X[i, j, k]), float(blk.Y[i, j, k]), float(blk.Z[i, j, k]))
                    )

    def _add_block_hexes(self, blocks: list[Block], blk_idx: int) -> None:
        """Generate hexahedral elements for all cells in a block."""
        idim, jdim, kdim = blocks[blk_idx].X.shape

        for i in range(1, idim):
            for j in range(1, jdim):
                for k in range(1, kdim):
                    self._element_id += 1
                    # Element format: id type ntags tag1 tag2 n0..n7
                    # type=5 (hex), ntags=2, tag1=physical_group(1=internal), tag2=-1
                    el: list[int] = [self._element_id, 5, 2, 1, -1]
                    for di, dj, dk in _HEX_SHIFTS:
                        el.append(
                            self._node_id(blocks, blk_idx, i + di, j + dj, k + dk)
                        )
                    self._elements.append(el)

    # ---- Boundary face generation ----

    def _next_group_id(self) -> int:
        return max(g[1] for g in self._groups) + 1

    def _add_boundary_quads(self, blocks: list[Block], bdry: Boundary) -> None:
        """Generate quad face elements for a named boundary.

        The quad node winding for each face ID follows the NMF spec:
          Face 1 (kmin): loop i(s1..e1), j(s2..e2) at k=0
          Face 2 (kmax): loop i(s1..e1), j(s2..e2) at k=kmax
          Face 3 (imin): loop j(s1..e1), k(s2..e2) at i=0
          Face 4 (imax): loop j(s1..e1), k(s2..e2) at i=imax
          Face 5 (jmin): loop k(s1..e1), i(s2..e2) at j=0
          Face 6 (jmax): loop k(s1..e1), i(s2..e2) at j=jmax
        """
        gid = self._next_group_id()
        self._groups.append((2, gid, bdry.name))

        blk_idx = bdry.block - 1  # convert to 0-based
        x_shape = blocks[blk_idx].X.shape
        imax = x_shape[0] - 1
        jmax = x_shape[1] - 1
        kmax = x_shape[2] - 1

        s1, e1, s2, e2 = bdry.s1, bdry.e1, bdry.s2, bdry.e2
        fid = bdry.face_id
        nid = self._node_id

        if fid == 1:
            # kmin face: primary=i, secondary=j, at k=0
            for j in range(s2 - 1, e2 - 1):
                for i in range(s1 - 1, e1 - 1):
                    self._element_id += 1
                    n1 = nid(blocks, blk_idx, i, j, 0)
                    n2 = nid(blocks, blk_idx, i + 1, j, 0)
                    n3 = nid(blocks, blk_idx, i + 1, j + 1, 0)
                    n4 = nid(blocks, blk_idx, i, j + 1, 0)
                    self._elements.append(
                        [self._element_id, 3, 2, gid, -1, n1, n2, n3, n4]
                    )

        elif fid == 2:
            # kmax face: primary=i, secondary=j, at k=kmax
            for j in range(s2 - 1, e2 - 1):
                for i in range(s1 - 1, e1 - 1):
                    self._element_id += 1
                    n1 = nid(blocks, blk_idx, i, j, kmax)
                    n2 = nid(blocks, blk_idx, i + 1, j, kmax)
                    n3 = nid(blocks, blk_idx, i + 1, j + 1, kmax)
                    n4 = nid(blocks, blk_idx, i, j + 1, kmax)
                    self._elements.append(
                        [self._element_id, 3, 2, gid, -1, n1, n2, n3, n4]
                    )

        elif fid == 3:
            # imin face: primary=j, secondary=k, at i=0
            for k in range(s2 - 1, e2 - 1):
                for j in range(s1 - 1, e1 - 1):
                    self._element_id += 1
                    n1 = nid(blocks, blk_idx, 0, j, k)
                    n2 = nid(blocks, blk_idx, 0, j + 1, k)
                    n3 = nid(blocks, blk_idx, 0, j + 1, k + 1)
                    n4 = nid(blocks, blk_idx, 0, j, k + 1)
                    self._elements.append(
                        [self._element_id, 3, 2, gid, -1, n1, n2, n3, n4]
                    )

        elif fid == 4:
            # imax face: primary=j, secondary=k, at i=imax
            for k in range(s2 - 1, e2 - 1):
                for j in range(s1 - 1, e1 - 1):
                    self._element_id += 1
                    n1 = nid(blocks, blk_idx, imax, j, k)
                    n2 = nid(blocks, blk_idx, imax, j + 1, k)
                    n3 = nid(blocks, blk_idx, imax, j + 1, k + 1)
                    n4 = nid(blocks, blk_idx, imax, j, k + 1)
                    self._elements.append(
                        [self._element_id, 3, 2, gid, -1, n1, n2, n3, n4]
                    )

        elif fid == 5:
            # jmin face: primary=k, secondary=i, at j=0
            for i in range(s2 - 1, e2 - 1):
                for k in range(s1 - 1, e1 - 1):
                    self._element_id += 1
                    n1 = nid(blocks, blk_idx, i, 0, k)
                    n2 = nid(blocks, blk_idx, i + 1, 0, k)
                    n3 = nid(blocks, blk_idx, i + 1, 0, k + 1)
                    n4 = nid(blocks, blk_idx, i, 0, k + 1)
                    self._elements.append(
                        [self._element_id, 3, 2, gid, -1, n1, n2, n3, n4]
                    )

        elif fid == 6:
            # jmax face: primary=k, secondary=i, at j=jmax
            for i in range(s2 - 1, e2 - 1):
                for k in range(s1 - 1, e1 - 1):
                    self._element_id += 1
                    n1 = nid(blocks, blk_idx, i, jmax, k)
                    n2 = nid(blocks, blk_idx, i + 1, jmax, k)
                    n3 = nid(blocks, blk_idx, i + 1, jmax, k + 1)
                    n4 = nid(blocks, blk_idx, i, jmax, k + 1)
                    self._elements.append(
                        [self._element_id, 3, 2, gid, -1, n1, n2, n3, n4]
                    )

        else:
            raise ValueError(f"Unknown face ID {fid} for boundary '{bdry.name}'")
