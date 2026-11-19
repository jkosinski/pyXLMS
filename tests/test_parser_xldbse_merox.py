#!/usr/bin/env python3

# pyXLMS - TESTS
# 2025 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

import pytest

MEROX1 = "data/merox/XLpeplib_Beveridge_QEx-HFX_DSS_R1.csv"
MEROX2 = "data/merox/XLpeplib_Beveridge_QEx-HFX_DSS_R1.zhrm"
MEROX3 = "data/merox/cas11.csv"
# excerpt of public PRIDE deposit PXD025064 (rows copied verbatim)
MEROX4 = "data/merox/PXD025064_AGAL_300DSBU_Fr7_excerpt.zhrm"


def test1():
    from pyXLMS.parser import read_merox
    from pyXLMS.transform import modifications_to_str as mts

    pr = read_merox(MEROX1, crosslinker="DSS")
    assert pr["data_type"] == "parser_result"
    assert pr["completeness"] == "partial"
    assert pr["search_engine"] == "MeroX"
    assert pr["crosslink-spectrum-matches"] is not None
    assert pr["crosslinks"] is None

    csms = pr["crosslink-spectrum-matches"]
    assert len(csms) == 93

    csm = csms[0]
    assert csm["data_type"] == "crosslink-spectrum-match"
    assert csm["completeness"] == "partial"
    assert csm["alpha_peptide"] == "GKSDNVPSEEVVK"
    assert mts(csm["alpha_modifications"]) == "(2:[DSS|138.06808])"
    assert csm["alpha_peptide_crosslink_position"] == 2
    assert csm["alpha_proteins"] == ["Cas9"]
    assert csm["alpha_proteins_crosslink_positions"] == [870]
    assert csm["alpha_proteins_peptide_positions"] == [869]
    assert csm["alpha_score"] == pytest.approx(13.14747)
    assert not csm["alpha_decoy"]
    assert csm["beta_peptide"] == "VKYVTEGMR"
    assert (
        mts(csm["beta_modifications"])
        == "(2:[DSS|138.06808]);(8:[Oxidation|15.994915])"
    )
    assert csm["beta_peptide_crosslink_position"] == 2
    assert csm["beta_proteins"] == ["Cas9"]
    assert csm["beta_proteins_crosslink_positions"] == [532]
    assert csm["beta_proteins_peptide_positions"] == [531]
    assert csm["beta_score"] == pytest.approx(12.143336)
    assert not csm["beta_decoy"]
    assert csm["crosslink_type"] == "intra"
    assert csm["score"] == pytest.approx(97.0)
    assert csm["spectrum_file"] == "XLpeplib_Beveridge_QEx-HFX_DSS_R1"
    assert csm["scan_nr"] == 10061
    assert csm["charge"] == 4
    assert csm["retention_time"] == pytest.approx(3099.0)
    assert csm["ion_mobility"] is None

    csm = csms[-1]
    assert csm["data_type"] == "crosslink-spectrum-match"
    assert csm["completeness"] == "partial"
    assert csm["alpha_peptide"] == "KLVDSTDK"
    assert mts(csm["alpha_modifications"]) == "(1:[DSS|138.06808])"
    assert csm["alpha_peptide_crosslink_position"] == 1
    assert csm["alpha_proteins"] == ["Cas9"]
    assert csm["alpha_proteins_crosslink_positions"] == [145]
    assert csm["alpha_proteins_peptide_positions"] == [145]
    assert csm["alpha_score"] == pytest.approx(77.979659)
    assert not csm["alpha_decoy"]
    assert csm["beta_peptide"] == "KLVDSTDK"
    assert mts(csm["beta_modifications"]) == "(1:[DSS|138.06808])"
    assert csm["beta_peptide_crosslink_position"] == 1
    assert csm["beta_proteins"] == ["Cas9"]
    assert csm["beta_proteins_crosslink_positions"] == [145]
    assert csm["beta_proteins_peptide_positions"] == [145]
    assert csm["beta_score"] == pytest.approx(77.979659)
    assert not csm["beta_decoy"]
    assert csm["crosslink_type"] == "intra"
    assert csm["score"] == pytest.approx(139)
    assert csm["spectrum_file"] == "XLpeplib_Beveridge_QEx-HFX_DSS_R1"
    assert csm["scan_nr"] == 8258
    assert csm["charge"] == 4
    assert csm["retention_time"] == pytest.approx(2605.0)
    assert csm["ion_mobility"] is None


def test2():
    from pyXLMS.parser import read_merox
    from pyXLMS.transform import modifications_to_str as mts

    pr = read_merox(MEROX2, crosslinker="DSS")
    assert pr["data_type"] == "parser_result"
    assert pr["completeness"] == "partial"
    assert pr["search_engine"] == "MeroX"
    assert pr["crosslink-spectrum-matches"] is not None
    assert pr["crosslinks"] is None

    csms = pr["crosslink-spectrum-matches"]
    assert len(csms) == 93

    csm = csms[0]
    assert csm["data_type"] == "crosslink-spectrum-match"
    assert csm["completeness"] == "partial"
    assert csm["alpha_peptide"] == "GKSDNVPSEEVVK"
    assert mts(csm["alpha_modifications"]) == "(2:[DSS|138.06808])"
    assert csm["alpha_peptide_crosslink_position"] == 2
    assert csm["alpha_proteins"] == ["Cas9"]
    assert csm["alpha_proteins_crosslink_positions"] == [870]
    assert csm["alpha_proteins_peptide_positions"] == [869]
    assert csm["alpha_score"] == pytest.approx(13.14747)
    assert not csm["alpha_decoy"]
    assert csm["beta_peptide"] == "VKYVTEGMR"
    assert (
        mts(csm["beta_modifications"])
        == "(2:[DSS|138.06808]);(8:[Oxidation|15.994915])"
    )
    assert csm["beta_peptide_crosslink_position"] == 2
    assert csm["beta_proteins"] == ["Cas9"]
    assert csm["beta_proteins_crosslink_positions"] == [532]
    assert csm["beta_proteins_peptide_positions"] == [531]
    assert csm["beta_score"] == pytest.approx(12.143336)
    assert not csm["beta_decoy"]
    assert csm["crosslink_type"] == "intra"
    assert csm["score"] == pytest.approx(97.0)
    assert csm["spectrum_file"] == "XLpeplib_Beveridge_QEx-HFX_DSS_R1"
    assert csm["scan_nr"] == 10061
    assert csm["charge"] == 4
    assert csm["retention_time"] == pytest.approx(3099.0)
    assert csm["ion_mobility"] is None

    csm = csms[-1]
    assert csm["data_type"] == "crosslink-spectrum-match"
    assert csm["completeness"] == "partial"
    assert csm["alpha_peptide"] == "KLVDSTDK"
    assert mts(csm["alpha_modifications"]) == "(1:[DSS|138.06808])"
    assert csm["alpha_peptide_crosslink_position"] == 1
    assert csm["alpha_proteins"] == ["Cas9"]
    assert csm["alpha_proteins_crosslink_positions"] == [145]
    assert csm["alpha_proteins_peptide_positions"] == [145]
    assert csm["alpha_score"] == pytest.approx(77.979659)
    assert not csm["alpha_decoy"]
    assert csm["beta_peptide"] == "KLVDSTDK"
    assert mts(csm["beta_modifications"]) == "(1:[DSS|138.06808])"
    assert csm["beta_peptide_crosslink_position"] == 1
    assert csm["beta_proteins"] == ["Cas9"]
    assert csm["beta_proteins_crosslink_positions"] == [145]
    assert csm["beta_proteins_peptide_positions"] == [145]
    assert csm["beta_score"] == pytest.approx(77.979659)
    assert not csm["beta_decoy"]
    assert csm["crosslink_type"] == "intra"
    assert csm["score"] == pytest.approx(139)
    assert csm["spectrum_file"] == "XLpeplib_Beveridge_QEx-HFX_DSS_R1"
    assert csm["scan_nr"] == 8258
    assert csm["charge"] == 4
    assert csm["retention_time"] == pytest.approx(2605.0)
    assert csm["ion_mobility"] is None


def test3():
    from pyXLMS.parser import read_merox
    from pyXLMS.transform import modifications_to_str as mts

    pr = read_merox(MEROX3, crosslinker="DSS")
    assert pr["data_type"] == "parser_result"
    assert pr["completeness"] == "partial"
    assert pr["search_engine"] == "MeroX"
    assert pr["crosslink-spectrum-matches"] is not None
    assert pr["crosslinks"] is None

    csms = pr["crosslink-spectrum-matches"]
    assert len(csms) == 93

    csm = csms[0]
    assert csm["data_type"] == "crosslink-spectrum-match"
    assert csm["completeness"] == "partial"
    assert csm["alpha_peptide"] == "GKSDNVPSEEVVK"
    assert mts(csm["alpha_modifications"]) == "(2:[DSS|138.06808])"
    assert csm["alpha_peptide_crosslink_position"] == 2
    assert csm["alpha_proteins"] == ["Cas9"]
    assert csm["alpha_proteins_crosslink_positions"] == [870]
    assert csm["alpha_proteins_peptide_positions"] == [869]
    assert csm["alpha_score"] == pytest.approx(13.14747)
    assert not csm["alpha_decoy"]
    assert csm["beta_peptide"] == "VKYVTEGMR"
    assert (
        mts(csm["beta_modifications"])
        == "(2:[DSS|138.06808]);(8:[Oxidation|15.994915])"
    )
    assert csm["beta_peptide_crosslink_position"] == 2
    assert csm["beta_proteins"] == ["Cas9"]
    assert csm["beta_proteins_crosslink_positions"] == [532]
    assert csm["beta_proteins_peptide_positions"] == [531]
    assert csm["beta_score"] == pytest.approx(12.143336)
    assert not csm["beta_decoy"]
    assert csm["crosslink_type"] == "intra"
    assert csm["score"] == pytest.approx(97.0)
    assert csm["spectrum_file"] == "XLpeplib_Beveridge_QEx-HFX_DSS_R1"
    assert csm["scan_nr"] == 10061
    assert csm["charge"] == 4
    assert csm["retention_time"] == pytest.approx(3099.0)
    assert csm["ion_mobility"] is None

    csm = csms[-1]
    assert csm["data_type"] == "crosslink-spectrum-match"
    assert csm["completeness"] == "partial"
    assert csm["alpha_peptide"] == "KLVDSTDK"
    assert mts(csm["alpha_modifications"]) == "(1:[DSS|138.06808])"
    assert csm["alpha_peptide_crosslink_position"] == 1
    assert csm["alpha_proteins"] == ["Cas9"]
    assert csm["alpha_proteins_crosslink_positions"] == [145]
    assert csm["alpha_proteins_peptide_positions"] == [145]
    assert csm["alpha_score"] == pytest.approx(77.979659)
    assert not csm["alpha_decoy"]
    assert csm["beta_peptide"] == "KLVDSTDK"
    assert mts(csm["beta_modifications"]) == "(1:[DSS|138.06808])"
    assert csm["beta_peptide_crosslink_position"] == 1
    assert csm["beta_proteins"] == ["Cas9"]
    assert csm["beta_proteins_crosslink_positions"] == [145]
    assert csm["beta_proteins_peptide_positions"] == [145]
    assert csm["beta_score"] == pytest.approx(77.979659)
    assert not csm["beta_decoy"]
    assert csm["crosslink_type"] == "intra"
    assert csm["score"] == pytest.approx(139)
    assert csm["spectrum_file"] == "XLpeplib_Beveridge_QEx-HFX_DSS_R1"
    assert csm["scan_nr"] == 8258
    assert csm["charge"] == 4
    assert csm["retention_time"] == pytest.approx(2605.0)
    assert csm["ion_mobility"] is None


def test4():
    from pyXLMS.parser import read_merox
    from pyXLMS.transform import modifications_to_str as mts

    with open(MEROX1, "r", encoding="utf-8") as f:
        pr = read_merox(f, crosslinker="DSS")
        f.close()

    assert pr["data_type"] == "parser_result"
    assert pr["completeness"] == "partial"
    assert pr["search_engine"] == "MeroX"
    assert pr["crosslink-spectrum-matches"] is not None
    assert pr["crosslinks"] is None

    csms = pr["crosslink-spectrum-matches"]
    assert len(csms) == 93

    csm = csms[0]
    assert csm["data_type"] == "crosslink-spectrum-match"
    assert csm["completeness"] == "partial"
    assert csm["alpha_peptide"] == "GKSDNVPSEEVVK"
    assert mts(csm["alpha_modifications"]) == "(2:[DSS|138.06808])"
    assert csm["alpha_peptide_crosslink_position"] == 2
    assert csm["alpha_proteins"] == ["Cas9"]
    assert csm["alpha_proteins_crosslink_positions"] == [870]
    assert csm["alpha_proteins_peptide_positions"] == [869]
    assert csm["alpha_score"] == pytest.approx(13.14747)
    assert not csm["alpha_decoy"]
    assert csm["beta_peptide"] == "VKYVTEGMR"
    assert (
        mts(csm["beta_modifications"])
        == "(2:[DSS|138.06808]);(8:[Oxidation|15.994915])"
    )
    assert csm["beta_peptide_crosslink_position"] == 2
    assert csm["beta_proteins"] == ["Cas9"]
    assert csm["beta_proteins_crosslink_positions"] == [532]
    assert csm["beta_proteins_peptide_positions"] == [531]
    assert csm["beta_score"] == pytest.approx(12.143336)
    assert not csm["beta_decoy"]
    assert csm["crosslink_type"] == "intra"
    assert csm["score"] == pytest.approx(97.0)
    assert csm["spectrum_file"] == "XLpeplib_Beveridge_QEx-HFX_DSS_R1"
    assert csm["scan_nr"] == 10061
    assert csm["charge"] == 4
    assert csm["retention_time"] == pytest.approx(3099.0)
    assert csm["ion_mobility"] is None

    csm = csms[-1]
    assert csm["data_type"] == "crosslink-spectrum-match"
    assert csm["completeness"] == "partial"
    assert csm["alpha_peptide"] == "KLVDSTDK"
    assert mts(csm["alpha_modifications"]) == "(1:[DSS|138.06808])"
    assert csm["alpha_peptide_crosslink_position"] == 1
    assert csm["alpha_proteins"] == ["Cas9"]
    assert csm["alpha_proteins_crosslink_positions"] == [145]
    assert csm["alpha_proteins_peptide_positions"] == [145]
    assert csm["alpha_score"] == pytest.approx(77.979659)
    assert not csm["alpha_decoy"]
    assert csm["beta_peptide"] == "KLVDSTDK"
    assert mts(csm["beta_modifications"]) == "(1:[DSS|138.06808])"
    assert csm["beta_peptide_crosslink_position"] == 1
    assert csm["beta_proteins"] == ["Cas9"]
    assert csm["beta_proteins_crosslink_positions"] == [145]
    assert csm["beta_proteins_peptide_positions"] == [145]
    assert csm["beta_score"] == pytest.approx(77.979659)
    assert not csm["beta_decoy"]
    assert csm["crosslink_type"] == "intra"
    assert csm["score"] == pytest.approx(139)
    assert csm["spectrum_file"] == "XLpeplib_Beveridge_QEx-HFX_DSS_R1"
    assert csm["scan_nr"] == 8258
    assert csm["charge"] == 4
    assert csm["retention_time"] == pytest.approx(2605.0)
    assert csm["ion_mobility"] is None


def test5():
    from pyXLMS.parser import read_merox
    from pyXLMS.transform import modifications_to_str as mts

    with open(MEROX2, "rb") as f:
        pr = read_merox(f, crosslinker="DSS")
        f.close()

    assert pr["data_type"] == "parser_result"
    assert pr["completeness"] == "partial"
    assert pr["search_engine"] == "MeroX"
    assert pr["crosslink-spectrum-matches"] is not None
    assert pr["crosslinks"] is None

    csms = pr["crosslink-spectrum-matches"]
    assert len(csms) == 93

    csm = csms[0]
    assert csm["data_type"] == "crosslink-spectrum-match"
    assert csm["completeness"] == "partial"
    assert csm["alpha_peptide"] == "GKSDNVPSEEVVK"
    assert mts(csm["alpha_modifications"]) == "(2:[DSS|138.06808])"
    assert csm["alpha_peptide_crosslink_position"] == 2
    assert csm["alpha_proteins"] == ["Cas9"]
    assert csm["alpha_proteins_crosslink_positions"] == [870]
    assert csm["alpha_proteins_peptide_positions"] == [869]
    assert csm["alpha_score"] == pytest.approx(13.14747)
    assert not csm["alpha_decoy"]
    assert csm["beta_peptide"] == "VKYVTEGMR"
    assert (
        mts(csm["beta_modifications"])
        == "(2:[DSS|138.06808]);(8:[Oxidation|15.994915])"
    )
    assert csm["beta_peptide_crosslink_position"] == 2
    assert csm["beta_proteins"] == ["Cas9"]
    assert csm["beta_proteins_crosslink_positions"] == [532]
    assert csm["beta_proteins_peptide_positions"] == [531]
    assert csm["beta_score"] == pytest.approx(12.143336)
    assert not csm["beta_decoy"]
    assert csm["crosslink_type"] == "intra"
    assert csm["score"] == pytest.approx(97.0)
    assert csm["spectrum_file"] == "XLpeplib_Beveridge_QEx-HFX_DSS_R1"
    assert csm["scan_nr"] == 10061
    assert csm["charge"] == 4
    assert csm["retention_time"] == pytest.approx(3099.0)
    assert csm["ion_mobility"] is None

    csm = csms[-1]
    assert csm["data_type"] == "crosslink-spectrum-match"
    assert csm["completeness"] == "partial"
    assert csm["alpha_peptide"] == "KLVDSTDK"
    assert mts(csm["alpha_modifications"]) == "(1:[DSS|138.06808])"
    assert csm["alpha_peptide_crosslink_position"] == 1
    assert csm["alpha_proteins"] == ["Cas9"]
    assert csm["alpha_proteins_crosslink_positions"] == [145]
    assert csm["alpha_proteins_peptide_positions"] == [145]
    assert csm["alpha_score"] == pytest.approx(77.979659)
    assert not csm["alpha_decoy"]
    assert csm["beta_peptide"] == "KLVDSTDK"
    assert mts(csm["beta_modifications"]) == "(1:[DSS|138.06808])"
    assert csm["beta_peptide_crosslink_position"] == 1
    assert csm["beta_proteins"] == ["Cas9"]
    assert csm["beta_proteins_crosslink_positions"] == [145]
    assert csm["beta_proteins_peptide_positions"] == [145]
    assert csm["beta_score"] == pytest.approx(77.979659)
    assert not csm["beta_decoy"]
    assert csm["crosslink_type"] == "intra"
    assert csm["score"] == pytest.approx(139)
    assert csm["spectrum_file"] == "XLpeplib_Beveridge_QEx-HFX_DSS_R1"
    assert csm["scan_nr"] == 8258
    assert csm["charge"] == 4
    assert csm["retention_time"] == pytest.approx(2605.0)
    assert csm["ion_mobility"] is None


def test6():
    from pyXLMS.parser import read_merox
    from pyXLMS.transform import modifications_to_str as mts

    with open(MEROX3, "r", encoding="utf-8") as f:
        pr = read_merox(f, crosslinker="DSS")
        f.close()

    assert pr["data_type"] == "parser_result"
    assert pr["completeness"] == "partial"
    assert pr["search_engine"] == "MeroX"
    assert pr["crosslink-spectrum-matches"] is not None
    assert pr["crosslinks"] is None

    csms = pr["crosslink-spectrum-matches"]
    assert len(csms) == 93

    csm = csms[0]
    assert csm["data_type"] == "crosslink-spectrum-match"
    assert csm["completeness"] == "partial"
    assert csm["alpha_peptide"] == "GKSDNVPSEEVVK"
    assert mts(csm["alpha_modifications"]) == "(2:[DSS|138.06808])"
    assert csm["alpha_peptide_crosslink_position"] == 2
    assert csm["alpha_proteins"] == ["Cas9"]
    assert csm["alpha_proteins_crosslink_positions"] == [870]
    assert csm["alpha_proteins_peptide_positions"] == [869]
    assert csm["alpha_score"] == pytest.approx(13.14747)
    assert not csm["alpha_decoy"]
    assert csm["beta_peptide"] == "VKYVTEGMR"
    assert (
        mts(csm["beta_modifications"])
        == "(2:[DSS|138.06808]);(8:[Oxidation|15.994915])"
    )
    assert csm["beta_peptide_crosslink_position"] == 2
    assert csm["beta_proteins"] == ["Cas9"]
    assert csm["beta_proteins_crosslink_positions"] == [532]
    assert csm["beta_proteins_peptide_positions"] == [531]
    assert csm["beta_score"] == pytest.approx(12.143336)
    assert not csm["beta_decoy"]
    assert csm["crosslink_type"] == "intra"
    assert csm["score"] == pytest.approx(97.0)
    assert csm["spectrum_file"] == "XLpeplib_Beveridge_QEx-HFX_DSS_R1"
    assert csm["scan_nr"] == 10061
    assert csm["charge"] == 4
    assert csm["retention_time"] == pytest.approx(3099.0)
    assert csm["ion_mobility"] is None

    csm = csms[-1]
    assert csm["data_type"] == "crosslink-spectrum-match"
    assert csm["completeness"] == "partial"
    assert csm["alpha_peptide"] == "KLVDSTDK"
    assert mts(csm["alpha_modifications"]) == "(1:[DSS|138.06808])"
    assert csm["alpha_peptide_crosslink_position"] == 1
    assert csm["alpha_proteins"] == ["Cas9"]
    assert csm["alpha_proteins_crosslink_positions"] == [145]
    assert csm["alpha_proteins_peptide_positions"] == [145]
    assert csm["alpha_score"] == pytest.approx(77.979659)
    assert not csm["alpha_decoy"]
    assert csm["beta_peptide"] == "KLVDSTDK"
    assert mts(csm["beta_modifications"]) == "(1:[DSS|138.06808])"
    assert csm["beta_peptide_crosslink_position"] == 1
    assert csm["beta_proteins"] == ["Cas9"]
    assert csm["beta_proteins_crosslink_positions"] == [145]
    assert csm["beta_proteins_peptide_positions"] == [145]
    assert csm["beta_score"] == pytest.approx(77.979659)
    assert not csm["beta_decoy"]
    assert csm["crosslink_type"] == "intra"
    assert csm["score"] == pytest.approx(139)
    assert csm["spectrum_file"] == "XLpeplib_Beveridge_QEx-HFX_DSS_R1"
    assert csm["scan_nr"] == 8258
    assert csm["charge"] == 4
    assert csm["retention_time"] == pytest.approx(2605.0)
    assert csm["ion_mobility"] is None


def test7():
    from pyXLMS.parser import read_merox

    err_str = (
        r"Cannot infer crosslinker mass because crosslinker is unknown\. "
        r"Please specify crosslinker mass manually!"
    )
    with pytest.raises(KeyError, match=err_str):
        _ = read_merox(MEROX1, crosslinker="DiSPASO")


def test8():
    # protein-terminal peptides: "{" protein N-terminus, "}" protein C-terminus
    from pyXLMS.parser import read_merox
    from pyXLMS.transform import modifications_to_str as mts

    pr = read_merox(MEROX4, crosslinker="DSBU")
    csms = pr["crosslink-spectrum-matches"]
    assert len(csms) == 3

    # [SPMFEKK] K6 - [KYLLS} K1
    csm = csms[0]
    assert csm["alpha_peptide"] == "KYLLS"
    assert csm["alpha_peptide_crosslink_position"] == 1
    assert csm["alpha_proteins_crosslink_positions"] == [135]
    assert csm["alpha_proteins_peptide_positions"] == [135]
    assert csm["beta_peptide"] == "SPMFEKK"
    assert csm["beta_proteins_crosslink_positions"] == [127]

    # {MIIYGVYIVSK] {0 - [TGISNIDK} K8
    csm = csms[1]
    assert csm["alpha_peptide"] == "MIIYGVYIVSK"
    assert csm["alpha_peptide_crosslink_position"] == 1
    assert mts(csm["alpha_modifications"]) == "(0:[DSBU|196.08479231])"
    assert csm["alpha_proteins"] == ["TRAPPC4"]
    assert csm["alpha_proteins_crosslink_positions"] == [1]
    assert csm["alpha_proteins_peptide_positions"] == [1]
    assert csm["beta_peptide"] == "TGISNIDK"
    assert csm["beta_peptide_crosslink_position"] == 8
    assert csm["beta_proteins_crosslink_positions"] == [219]
    assert csm["beta_proteins_peptide_positions"] == [212]
    assert csm["scan_nr"] == 21256

    # [LQELLAQVEKTGISNIDK} K10 - {MIIYGVYIVSK] {0
    csm = csms[2]
    assert csm["alpha_peptide"] == "LQELLAQVEKTGISNIDK"
    assert csm["alpha_proteins_crosslink_positions"] == [211]
    assert csm["beta_peptide"] == "MIIYGVYIVSK"
    assert csm["beta_peptide_crosslink_position"] == 1
    assert csm["beta_proteins_crosslink_positions"] == [1]
    assert csm["beta_proteins_peptide_positions"] == [1]


@pytest.mark.parametrize(
    "position_str,from_position,sequence,xl_position,pep_position",
    [
        ("K2", 869, "[GKSDNVPSEEVVK]", 2, 869),
        ("{0", 0, "{MIIYGVYIVSK]", 1, 1),
        # protein N-terminal peptide crosslinked at a residue, e.g. in PRIDE deposit PXD008348
        ("K3", 0, "{MSKTAVK]", 3, 1),
        ("K1", 135, "[KYLLS}", 1, 135),
    ],
)
def test9(position_str, from_position, sequence, xl_position, pep_position):
    from pyXLMS.parser._parser_xldbse_merox import __get_merox_crosslink_position
    from pyXLMS.parser._parser_xldbse_merox import __get_merox_peptide_position

    assert __get_merox_crosslink_position(position_str) == xl_position
    assert __get_merox_peptide_position(from_position, sequence) == pep_position
