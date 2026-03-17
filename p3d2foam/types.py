"""Data model for p3d2foam conversions."""

from __future__ import annotations

from dataclasses import dataclass, field

# NMF face ID constants (official spec)
# Face 1: kmin  primary=i  secondary=j
# Face 2: kmax  primary=i  secondary=j
# Face 3: imin  primary=j  secondary=k
# Face 4: imax  primary=j  secondary=k
# Face 5: jmin  primary=k  secondary=i
# Face 6: jmax  primary=k  secondary=i
FACE_K_MIN = 1
FACE_K_MAX = 2
FACE_I_MIN = 3
FACE_I_MAX = 4
FACE_J_MIN = 5
FACE_J_MAX = 6


@dataclass
class Boundary:
    """A named boundary face region from the NMF file.

    All indices are 1-based per the NMF convention.
    """

    name: str
    block: int  # 1-based block number
    face_id: int  # 1-6 per NMF face convention
    s1: int  # start index in primary coordinate direction
    e1: int  # end index in primary coordinate direction
    s2: int  # start index in secondary coordinate direction
    e2: int  # end index in secondary coordinate direction


@dataclass
class InterfacePair:
    """A pair of boundaries forming a block-block interface.

    Created from ONE-TO-ONE or Patched entries in the NMF file.
    """

    boundary_a: Boundary
    boundary_b: Boundary
    swap: bool = False  # TRUE if primary directions are crossed between blocks


@dataclass
class ConversionConfig:
    """Full configuration for a Plot3D to OpenFOAM conversion."""

    p3d_file: str
    nmf_file: str
    case_dir: str = "."
    interfaces: list[tuple[str, str]] | None = None
    scale: tuple[float, float, float] | None = None
    stitch_tolerances: dict | None = None
    binary: bool = False
    big_endian: bool = False
    fortran: bool = False  # Fortran unformatted binary (.ufmt)
