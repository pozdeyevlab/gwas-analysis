"""Basic tests of the bed builder module."""

from pathlib import Path
from tempfile import NamedTemporaryFile as NamedTemp
import pandas as pd

import pytest

from modules.align_non_palindromic_variants import harmonize


@pytest.mark.parametrize(
    "study_id_exact, study_id_inverse, gnoamd_id_one, gnomad_id_two, AF_one, AF_two, study_AF, alignment",
    [
        (
            "1:794332:G:A",
            "1:794332:G:A",
            "1:794332:G:A",
            "1:794332:G:A",
            "0.001",
            "0.85",
            "0.01",
            "exact_match",
        ),
    ],
)
def test_build_bed(
    snpid: str,
    expected_chr: str,
    expected_start: str,
    expected_end: str,
    tmp_path: Path,
) -> None:
    """
    Creates temp input file, creates bed and tests output
    """
    data_dict = {"SNPID": [snpid]}
    data = pd.DataFrame(data_dict)
    with NamedTemp(suffix=".out", dir=tmp_path, mode="w", delete=True) as out:
        data.to_csv(out.name, sep=" ")

        with NamedTemp(suffix=".bed", dir=tmp_path, mode="w", delete=True) as bed:
            make_bed(input=out.name, output=bed.name)

            with open(bed.name, "r") as read_bed:
                for line in read_bed.readlines():
                    print(line)
                    assert line.split(" ")[0] == expected_chr
                    assert line.split(" ")[1] == expected_start
                    assert line.split(" ")[2] == expected_end
                    assert line.split(" ")[3].rstrip() == snpid
