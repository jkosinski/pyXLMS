#!/usr/bin/env python3

# pyXLMS - TESTS

import pytest


def test1():
    # residue pair notation
    from pyXLMS.parser import parse_site_from_scout

    site = parse_site_from_scout("IDSLSAQLSQLQ [K] QLAAK", 13, "sp|P02545|LMNA", "311")
    assert site["peptide"] == "IDSLSAQLSQLQKQLAAK"
    assert site["peptide_position"] == 13
    assert site["proteins"] == ["sp|P02545|LMNA"]
    assert site["protein_positions"] == [311]


def test2():
    # filtered CSM notation with several proteins
    from pyXLMS.parser import parse_site_from_scout

    site = parse_site_from_scout(
        "VTKGSLEEEK", "3", "sp|P0DMV9|HSPA1B; sp|P0DMV8|HSPA1A", "159; 159"
    )
    assert site["peptide"] == "VTKGSLEEEK"
    assert site["peptide_position"] == 3
    assert site["proteins"] == ["sp|P0DMV9|HSPA1B", "sp|P0DMV8|HSPA1A"]
    assert site["protein_positions"] == [159, 159]


def test3():
    from pyXLMS.parser import parse_site_from_scout

    with pytest.raises(ValueError):
        _ = parse_site_from_scout("KLGGTIDDCELVEGLVLTQK", 1, "sp|P50991|CCT4", "---")
