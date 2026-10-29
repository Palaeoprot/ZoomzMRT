"""
Unit tests for 2-column text and mMass MSD spectrum ingestion.
"""

import base64
import zlib
from pathlib import Path

import numpy as np
import pytest

from zoomzmrt.msd_reader import parse_mmass_msd
from zoomzmrt.txt_reader import parse_summed_text_spectrum


def test_parse_summed_text_spectrum(tmp_path: Path):
    txt_file = tmp_path / "sample_spot_A1_LMC.txt"
    # Create test 2-column data
    lines = [
        "m/z\tIntensity\n",
        "1105.5754\t50000.0\n",
        "1106.5594\t25000.0\n",
        "1106.5841\t25000.0\n",
        "2019.9736\t80000.0\n",
    ]
    txt_file.write_text("".join(lines), encoding="utf-8")

    rec, qc = parse_summed_text_spectrum(
        file_path=txt_file,
        dataset_id="test_txt_ds",
    )

    assert rec["n_peaks"] == 4
    assert qc.is_precalibrated is True  # Detected from _LMC
    assert abs(qc.base_peak_intensity - 80000.0) < 1e-3
    assert abs(qc.base_peak_mz - 2019.9736) < 1e-3
    assert len(rec["mz"]) == 4


def test_parse_mmass_msd(tmp_path: Path):
    msd_file = tmp_path / "test_sample.msd"

    # Encode float32 arrays
    mzs = np.array([1105.5754, 1106.5594, 2019.9736], dtype="<f4")
    ints = np.array([50000.0, 25000.0, 80000.0], dtype="<f4")

    mz_b64 = base64.b64encode(zlib.compress(mzs.tobytes())).decode("ascii")
    int_b64 = base64.b64encode(zlib.compress(ints.tobytes())).decode("ascii")

    xml_content = f"""<?xml version="1.0" encoding="utf-8" ?>
<mSD version="2.2">
  <description>
    <title>ROM_1555</title>
  </description>
  <spectrum points="3" msLevel="1" precision="32">
    <mzArray precision="32" compression="zlib" endian="little">{mz_b64}</mzArray>
    <intArray precision="32" compression="zlib" endian="little">{int_b64}</intArray>
  </spectrum>
  <peaklist>
    <peak mz="1105.5754" intensity="50000.0" sn="25.0" />
    <peak mz="2019.9736" intensity="80000.0" sn="40.0" />
  </peaklist>
</mSD>
"""
    msd_file.write_text(xml_content, encoding="utf-8")

    rec, qc = parse_mmass_msd(
        file_path=msd_file,
        dataset_id="test_msd_ds",
    )

    assert rec["sample_id"] == "ROM_1555"
    assert rec["n_peaks"] == 3
    assert qc.profile_points == 3
    assert qc.peaklist_peaks == 2
    assert abs(qc.base_peak_intensity - 80000.0) < 1e-1
