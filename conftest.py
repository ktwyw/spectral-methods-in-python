"""Lets `pytest tests/` import spectral_tools from the repository root."""
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).parent))
