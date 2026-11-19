#!/usr/bin/env python3

# pyXLMS - TESTS

import pytest


def test1():
    from pyXLMS.parser import parse_site_from_merox

    site = parse_site_from_merox("[VKYVTEGmR]", ">Cas9", 531, "K2")
    assert site["peptide"] == "VKYVTEGMR"
    assert site["peptide_position"] == 2
    assert site["proteins"] == ["Cas9"]
    assert site["protein_positions"] == [532]


def test2():
    # additional proteins are listed without position and are not parsed
    from pyXLMS.parser import parse_site_from_merox

    site = parse_site_from_merox("[GGYDSPR]", ">npl3 and >npl3@345-351", "345", "S5")
    assert site["peptide"] == "GGYDSPR"
    assert site["peptide_position"] == 5
    assert site["proteins"] == ["npl3 and >npl3@345-351"]
    assert site["protein_positions"] == [349]


def test3():
    from pyXLMS.parser import parse_site_from_merox

    site = parse_site_from_merox("[GKSDNVPSEEVVK]", ">REV__Cas9", "869", "K2")
    assert site["proteins"] == ["REV__Cas9"]
    assert site["protein_positions"] == [870]


def test4():
    from pyXLMS.parser import parse_site_from_merox

    with pytest.raises(KeyError, match="Key x not found in parameter 'modifications'"):
        _ = parse_site_from_merox("[LTTNLExTLPESVEILGPQNK]", ">Ssz1", 346, "x7")

    site = parse_site_from_merox(
        "[LTTNLExTLPESVEILGPQNK]", ">Ssz1", 346, "x7", parse_modifications=False
    )
    assert site["peptide"] == "LTTNLExTLPESVEILGPQNK"
    assert site["protein_positions"] == [352]


def test5():
    from pyXLMS.parser import parse_site_from_merox

    with pytest.raises(RuntimeError, match="Could not parse position from K"):
        _ = parse_site_from_merox("[GKSDNVPSEEVVK]", ">Cas9", 869, "K")
