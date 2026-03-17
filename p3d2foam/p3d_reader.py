"""Read Plot3D formatted (ASCII) grid files.

Handles both standard multi-line dimension headers (one block per line)
and single-line dimension headers (all blocks on one line), as well as
Fortran repeat notation (e.g., ``3578*0.0``).

Uses ``plot3d.Block`` from Plot3D_utilities for the output data structure.
Falls back to ``plot3d.read_plot3D`` for binary files.
"""

from __future__ import annotations

import logging
import re
from pathlib import Path

import numpy as np

logger = logging.getLogger(__name__)


def read_plot3d(
    filename: str | Path,
    binary: bool = False,
    big_endian: bool = False,
) -> list:
    """Read a multi-block Plot3D grid file.

    For binary files, delegates to ``plot3d.read_plot3D``.
    For ASCII files, handles both dimension header formats and
    Fortran repeat notation.

    Returns a list of ``plot3d.Block`` objects.
    """
    if binary:
        from plot3d import read_plot3D as _read_binary
        return _read_binary(str(filename), binary=True, big_endian=big_endian)

    return _read_ascii(Path(filename))


def _read_ascii(path: Path) -> list:
    """Read an ASCII Plot3D file, returning list of plot3d.Block."""
    from plot3d import Block

    logger.info("Reading ASCII Plot3D: %s", path)

    with open(path) as fp:
        raw = fp.read()

    # Expand Fortran repeat notation (e.g., "3578*0.0" -> "0.0 0.0 0.0 ...")
    # and replace commas with spaces
    text = _expand_repeats(raw)
    tokens = text.split()

    idx = 0

    # Number of blocks
    nblocks = int(tokens[idx])
    idx += 1

    # Block dimensions: may be all on one line or one per block
    dims: list[tuple[int, int, int]] = []
    for _ in range(nblocks):
        idim = int(tokens[idx])
        jdim = int(tokens[idx + 1])
        kdim = int(tokens[idx + 2])
        dims.append((idim, jdim, kdim))
        idx += 3

    # Read coordinate arrays
    blocks = []
    for b, (idim, jdim, kdim) in enumerate(dims):
        n = idim * jdim * kdim
        coords = []
        for coord_name in ("X", "Y", "Z"):
            vals = np.array(tokens[idx : idx + n], dtype=np.float64)
            idx += n
            # File order: i varies fastest (innermost), k slowest
            # (standard Fortran column-major for Plot3D)
            # Reshape to (KMAX, JMAX, IMAX) then transpose to (IMAX, JMAX, KMAX)
            arr = vals.reshape((kdim, jdim, idim)).transpose(2, 1, 0)
            coords.append(arr)

        blocks.append(Block(coords[0], coords[1], coords[2]))
        logger.debug(
            "Block %d: (%d, %d, %d) = %d nodes",
            b + 1, idim, jdim, kdim, n,
        )

    logger.info("Read %d blocks, %d total tokens consumed", len(blocks), idx)
    return blocks


# Pattern for Fortran repeat notation: "N*value" where N is an integer
_REPEAT_RE = re.compile(r"(\d+)\*([^\s,]+)")


def _expand_repeats(text: str) -> str:
    """Expand Fortran repeat notation and normalize separators.

    Converts patterns like ``3578*0.0`` to ``0.0 0.0 0.0 ...`` (3578 times).
    Also replaces commas with spaces.
    """
    # Replace commas with spaces first
    text = text.replace(",", " ")

    # Find all repeat patterns and expand them
    def _replacer(match: re.Match) -> str:
        count = int(match.group(1))
        value = match.group(2)
        return " ".join([value] * count)

    return _REPEAT_RE.sub(_replacer, text)
