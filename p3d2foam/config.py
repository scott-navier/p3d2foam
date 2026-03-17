"""YAML configuration loader for p3d2foam."""

from __future__ import annotations

from pathlib import Path

import yaml

from .types import ConversionConfig


def load_config(config_path: str | Path) -> ConversionConfig:
    """Load a conversion configuration from a YAML file.

    Example YAML::

        p3d_file: m6wing.p3dfmt
        nmf_file: m6wing.nmf
        case_dir: ./mesh
        binary: false
        scale: [1.1963, 1.1963, 1.1963]
        interfaces:
          - [B1KM, B4KM]
          - [B1J1, B4J1]
    """
    with open(config_path) as f:
        data = yaml.safe_load(f)

    interfaces = None
    if "interfaces" in data:
        interfaces = [tuple(pair) for pair in data["interfaces"]]

    scale = None
    if "scale" in data:
        scale = tuple(data["scale"])

    return ConversionConfig(
        p3d_file=data["p3d_file"],
        nmf_file=data["nmf_file"],
        case_dir=data.get("case_dir", "."),
        interfaces=interfaces,
        scale=scale,
        stitch_tolerances=data.get("stitch_tolerances"),
        binary=data.get("binary", False),
        big_endian=data.get("big_endian", False),
    )
