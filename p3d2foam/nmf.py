"""Parser for NASA's Neutral Map File (NMF) format.

Spec: https://tmbwg.github.io/turbmodels/nmf_documentation.html

The NMF describes boundary conditions and block-to-block connectivity
for multi-block structured grids in Plot3D format.
"""

from __future__ import annotations

import logging
from pathlib import Path

from .types import Boundary, InterfacePair

logger = logging.getLogger(__name__)

# NMF boundary types that define block-to-block interfaces (two-sided).
# These have 14 fields: Type B1 F1 S1 E1 S2 E2 B2 F2 S1 E1 S2 E2 Swap
_INTERFACE_TYPES = frozenset({"ONE-TO-ONE", "ONE_TO_ONE", "PATCHED"})


class NeutralMapFile:
    """Parse and represent a NASA Neutral Map File.

    Attributes:
        nblocks: Number of blocks in the grid.
        block_dims: List of (IDIM, JDIM, KDIM) tuples for each block.
        boundaries: Named physical boundaries (single-sided).
        interfaces: Block-to-block interface pairs (ONE-TO-ONE or Patched).
    """

    def __init__(self, filename: str | Path) -> None:
        self._nblocks: int = 0
        self._block_dims: list[tuple[int, int, int]] = []
        self._boundaries: list[Boundary] = []
        self._interfaces: list[InterfacePair] = []
        self._parse(Path(filename))

    @property
    def nblocks(self) -> int:
        return self._nblocks

    @property
    def block_dims(self) -> list[tuple[int, int, int]]:
        return self._block_dims

    @property
    def boundaries(self) -> list[Boundary]:
        return self._boundaries

    @property
    def interfaces(self) -> list[InterfacePair]:
        return self._interfaces

    def _parse(self, path: Path) -> None:
        # Read all lines, strip trailing backslashes (common NMF artifact),
        # then process as a flat list of cleaned lines.
        with open(path) as fp:
            raw_lines = fp.readlines()

        lines: list[str] = []
        for raw in raw_lines:
            cleaned = raw.rstrip("\n\r")
            if cleaned.endswith("\\"):
                cleaned = cleaned[:-1]
            lines.append(cleaned)

        idx = self._parse_blocks(lines, 0)
        self._parse_boundaries(lines, idx)

    def _parse_blocks(self, lines: list[str], start: int) -> int:
        """Parse block dimensions section. Returns index of next unparsed line."""
        idx = start

        # Skip comments and blanks to find the block count
        while idx < len(lines):
            s = lines[idx].strip()
            if s and not s.startswith("#"):
                break
            idx += 1

        self._nblocks = int(lines[idx].strip())
        idx += 1

        # Skip blanks between count and dimension rows
        while idx < len(lines) and not lines[idx].strip():
            idx += 1

        # Each block: Block# IDIM JDIM KDIM
        for _ in range(self._nblocks):
            while idx < len(lines) and not lines[idx].strip():
                idx += 1
            parts = lines[idx].split()
            if len(parts) >= 4:
                idim, jdim, kdim = int(parts[1]), int(parts[2]), int(parts[3])
            elif len(parts) == 3:
                idim, jdim, kdim = int(parts[0]), int(parts[1]), int(parts[2])
            else:
                raise ValueError(f"Cannot parse block dimensions from: {lines[idx]!r}")
            self._block_dims.append((idim, jdim, kdim))
            idx += 1

        logger.debug("Parsed %d blocks: %s", self._nblocks, self._block_dims)
        return idx

    def _parse_boundaries(self, lines: list[str], start: int) -> None:
        """Parse the topology & boundary conditions section."""
        interface_counter = 0

        for idx in range(start, len(lines)):
            line = lines[idx].strip()
            if not line or line.startswith("#"):
                continue

            fields = line.split()
            if not fields:
                continue

            # First field is the type/name — strip quotes if present
            type_or_name = fields[0].strip("'\"")

            if type_or_name.upper() in _INTERFACE_TYPES:
                self._parse_interface(fields, type_or_name, interface_counter)
                interface_counter += 1
            else:
                self._parse_single_boundary(fields, type_or_name)

        logger.debug(
            "Parsed %d boundaries, %d interfaces",
            len(self._boundaries),
            len(self._interfaces),
        )

    def _parse_interface(
        self, fields: list[str], itype: str, counter: int
    ) -> None:
        """Parse a ONE-TO-ONE or Patched interface line.

        Format: Type B1 F1 S1 E1 S2 E2 B2 F2 S1 E1 S2 E2 [Swap]
        """
        ints = [int(f) for f in fields[1:13]]
        b1, f1, s1a, e1a, s2a, e2a = ints[0:6]
        b2, f2, s1b, e1b, s2b, e2b = ints[6:12]

        swap = False
        if len(fields) > 13:
            swap = fields[13].upper() == "TRUE"

        name_a = f"interface{counter}a_B{b1}F{f1}"
        name_b = f"interface{counter}b_B{b2}F{f2}"

        boundary_a = Boundary(
            name=name_a, block=b1, face_id=f1, s1=s1a, e1=e1a, s2=s2a, e2=e2a
        )
        boundary_b = Boundary(
            name=name_b, block=b2, face_id=f2, s1=s1b, e1=e1b, s2=s2b, e2=e2b
        )

        self._interfaces.append(InterfacePair(boundary_a, boundary_b, swap))

        logger.debug(
            "Interface (%s): %s <-> %s (swap=%s)", itype, name_a, name_b, swap
        )

    def _parse_single_boundary(
        self, fields: list[str], name: str
    ) -> None:
        """Parse a single-sided boundary line.

        Format: Name/Type B1 F1 S1 E1 S2 E2
        """
        ints = [int(f) for f in fields[1:7]]
        block, face_id, s1, e1, s2, e2 = ints

        self._boundaries.append(
            Boundary(name=name, block=block, face_id=face_id, s1=s1, e1=e1, s2=s2, e2=e2)
        )

    def __str__(self) -> str:
        return (
            f"NeutralMapFile: {self._nblocks} blocks, "
            f"{len(self._boundaries)} boundaries, "
            f"{len(self._interfaces)} interfaces"
        )
