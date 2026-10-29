"""
Unit tests for interactive HTML dashboard generator.
"""

from pathlib import Path

import pandas as pd
import pytest

from zoomzmrt.dashboard import generate_interactive_dashboard_html, parse_spot_coordinate


def test_parse_spot_coordinate():
    assert parse_spot_coordinate("20260327MRT_KR71_0_K10_1") == ("K", 10)
    assert parse_spot_coordinate("sample_B19") == ("B", 19)
    assert parse_spot_coordinate("E5") == ("E", 5)
    assert parse_spot_coordinate("ROM_1555") is None


def test_generate_interactive_dashboard_html(tmp_path: Path):
    pqi_csv = tmp_path / "pqi_report.csv"
    plate_map_csv = tmp_path / "plate_map.csv"
    out_html = tmp_path / "dashboard.html"

    # Create synthetic pqi_report.csv
    df_pqi = pd.DataFrame([
        {
            "sample_id": "test_spot_A1_1",
            "n_markers_detected": 6,
            "pqi_mean": 0.75,
            "deamidation_mean": 0.25,
            "COL1A1_P1105_pqi": 0.70,
            "COL1A1_P1105_deam": 0.30,
            "COL1A1_P1105_flag": "OK",
        },
        {
            "sample_id": "test_spot_A3_1",
            "n_markers_detected": 6,
            "pqi_mean": 0.71,
            "deamidation_mean": 0.29,
            "COL1A1_P1105_pqi": 0.68,
            "COL1A1_P1105_deam": 0.32,
            "COL1A1_P1105_flag": "OK",
        },
    ])
    df_pqi.to_csv(pqi_csv, index=False)

    # Create synthetic plate map
    df_pm = pd.DataFrame([
        {
            "SampleID": "ROM_1001",
            "Spot1": "A1",
            "Spot2": "A3",
        }
    ])
    df_pm.to_csv(plate_map_csv, index=False)

    html_file = generate_interactive_dashboard_html(
        pqi_report_path=pqi_csv,
        output_html_path=out_html,
        plate_map_path=plate_map_csv,
        dataset_title="Test Dashboard",
    )

    assert html_file.exists()
    content = html_file.read_text(encoding="utf-8")
    assert "<!DOCTYPE html>" in content
    assert "Test Dashboard" in content
    assert "ROM_1001" in content
    assert "plotly" in content.lower()
