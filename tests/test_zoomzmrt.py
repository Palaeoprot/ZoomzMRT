"""
Tests for ZoomzMRT package.
"""

from pathlib import Path
import numpy as np
import pyarrow.parquet as pq
import pytest

from zoomzmrt.fticr_reader import parse_fticr_mzml, _fticr_centroid
from zoomzmrt.mrt_reader import parse_waters_mrt_mzml, _high_res_centroid
from zoomzmrt.cli import run_pipeline


def test_high_res_centroid_synthetic():
    mz1 = np.linspace(1570.670, 1570.684, 15)
    int1 = np.exp(-0.5 * ((mz1 - 1570.6774) / 0.002) ** 2) * 10000.0

    mz2 = np.linspace(1570.686, 1570.700, 15)
    int2 = np.exp(-0.5 * ((mz2 - 1570.6930) / 0.002) ** 2) * 8000.0

    mzs = np.concatenate([mz1, mz2])
    ints = np.concatenate([int1, int2])

    c_mz, c_int = _high_res_centroid(mzs, ints, bin_size_da=0.003, sn_threshold=2.0)
    assert len(c_mz) == 2
    assert abs(c_mz[0] - 1570.6774) < 0.0005
    assert abs(c_mz[1] - 1570.6930) < 0.0005


def test_fticr_fine_isotopic_preservation():
    mz_arr = np.array([1200.500, 1200.501, 1200.502, 1200.508, 1200.509, 1200.510])
    int_arr = np.array([5000.0, 10000.0, 5000.0, 4000.0, 8000.0, 4000.0])

    c_mz, c_int = _fticr_centroid(mz_arr, int_arr, sn_threshold=2.0)
    assert len(c_mz) == 2
    assert abs(c_mz[0] - 1200.501) < 0.0005
    assert abs(c_mz[1] - 1200.509) < 0.0005


def test_direct_pqi_estimation():
    from zoomzmrt.deamidation import compute_high_res_deamidation

    # Simulate P1105 undeamidated (1105.5807) and deamidated (1106.5647) and 13C (1106.5841)
    mzs = [1105.5807, 1106.5647, 1106.5841]
    ints = [10000.0, 5000.0, 500.0]  # 5000/(10000+5000) = 33.3% deamidation, PQI = 66.7%

    summary = compute_high_res_deamidation(mzs, ints, sample_id="test_bone", dataset_id="test_ds")
    assert summary.n_markers_detected >= 1

    p1105 = [r for r in summary.results if r.marker_name == "COL1A1_P1105"][0]
    assert p1105.pqi_fraction is not None
    assert abs(p1105.pqi_fraction - (10000.0 / 15000.0)) < 0.001
    assert abs(p1105.deamidation_fraction - (5000.0 / 15000.0)) < 0.001


def test_zoomzmrt_pipeline(tmp_path):
    sample_mzml = Path(r"C:\Users\matth\Downloads\NickMitchell\extracted\ExampleInput-20260927T070858Z-1-001\ExampleInput\110924_h12_glufib.mzML")
    if not sample_mzml.exists():
        pytest.skip("Sample mzML file not found")

    res = run_pipeline(
        input_path=sample_mzml,
        dataset_id="Test_ZoomzMRT_Dataset",
        output_root=tmp_path,
        instrument_type="mrt",
        lockmass_name="glufib",
    )

    assert res["n_rows"] == 1
    p_path = Path(res["parquet_path"])
    assert p_path.exists()
    table = pq.read_table(p_path)
    assert len(table) == 1
    assert "mz" in table.column_names
    assert "intensity" in table.column_names


