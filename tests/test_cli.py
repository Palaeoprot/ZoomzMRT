"""
Unit tests for ZoomzMRT CLI pipeline execution.
"""

from pathlib import Path
import pyarrow.parquet as pq
import pytest

from zoomzmrt.cli import detect_instrument_type, run_pipeline


def test_detect_instrument_type():
    mrt_files = [Path("sample_glufib_mrt_01.mzML"), Path("sample_glufib_mrt_02.mzML")]
    assert detect_instrument_type(mrt_files) == "mrt"

    fticr_files = [Path("solarix_fticr_bone_01.mzML")]
    assert detect_instrument_type(fticr_files) == "fticr"

    with pytest.raises(ValueError, match="Could not automatically detect"):
        detect_instrument_type([Path("unknown_spectrum.mzML")])
