"""CLI entry point for p3d2foam.

Usage:
    python -m p3d2foam grid.p3dfmt                     # auto-find .nmf
    python -m p3d2foam grid.p3dfmt -c config.yaml      # use YAML config
    python -m p3d2foam grid.p3dfmt --msh-only           # just produce .msh
    python -m p3d2foam grid.p3d --binary --scale 1.2 1.2 1.2
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(
        prog="p3d2foam",
        description="Convert NASA TMR Plot3D grids to OpenFOAM format",
    )
    parser.add_argument("p3d_file", help="Plot3D grid file (ASCII or binary)")
    parser.add_argument(
        "-m", "--nmf-file",
        help="Neutral Map File (default: <p3d_file_stem>.nmf)",
    )
    parser.add_argument(
        "-c", "--config",
        help="YAML config file for interfaces, scale, etc.",
    )
    parser.add_argument(
        "-d", "--case-dir",
        default=".",
        help="OpenFOAM case directory (default: current dir)",
    )
    parser.add_argument(
        "-o", "--output",
        help="Output .msh filename (for --msh-only mode)",
    )
    parser.add_argument(
        "--msh-only",
        action="store_true",
        help="Only generate .msh file, skip OpenFOAM pipeline",
    )
    parser.add_argument(
        "--binary",
        action="store_true",
        help="Read Plot3D file as binary (default: ASCII)",
    )
    parser.add_argument(
        "--big-endian",
        action="store_true",
        help="Binary file is big-endian",
    )
    parser.add_argument(
        "--fortran",
        action="store_true",
        help="Read Fortran unformatted binary (.ufmt) with record markers",
    )
    parser.add_argument(
        "--scale",
        type=float,
        nargs=3,
        metavar=("SX", "SY", "SZ"),
        help="Scale factor for transformPoints",
    )
    parser.add_argument(
        "-v", "--verbose",
        action="store_true",
        help="Enable verbose/debug logging",
    )

    args = parser.parse_args(argv)

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(levelname)s: %(message)s",
    )

    # Resolve NMF file
    nmf_file = args.nmf_file
    if nmf_file is None:
        nmf_file = str(Path(args.p3d_file).with_suffix(".nmf"))

    if args.config:
        from .config import load_config

        config = load_config(args.config)
        # CLI args override config file
        config.p3d_file = args.p3d_file
        config.nmf_file = nmf_file
        if args.case_dir != ".":
            config.case_dir = args.case_dir
        if args.binary:
            config.binary = True
        if args.big_endian:
            config.big_endian = True
        if args.fortran:
            config.fortran = True
        if args.scale:
            config.scale = tuple(args.scale)
    else:
        from .types import ConversionConfig

        config = ConversionConfig(
            p3d_file=args.p3d_file,
            nmf_file=nmf_file,
            case_dir=args.case_dir,
            binary=args.binary,
            big_endian=args.big_endian,
            fortran=args.fortran,
            scale=tuple(args.scale) if args.scale else None,
        )

    if args.msh_only:
        from .gmsh_writer import GmshWriter
        from .nmf import NeutralMapFile
        from .p3d_reader import read_plot3d

        blocks = read_plot3d(config.p3d_file, binary=config.binary, big_endian=config.big_endian, fortran=config.fortran)
        nmf = NeutralMapFile(config.nmf_file)

        writer = GmshWriter()
        writer.build(blocks, nmf)

        out_path = args.output or str(Path(config.p3d_file).with_suffix(".msh"))
        writer.write(out_path)
        logging.info("Wrote %s", out_path)
    else:
        from .pipeline import FoamPipeline

        pipeline = FoamPipeline(config)
        pipeline.run()


if __name__ == "__main__":
    main()
