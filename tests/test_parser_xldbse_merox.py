#!/usr/bin/env python3

# pyXLMS - TESTS
# 2025 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

import pytest

MEROX1 = "data/merox/XLpeplib_Beveridge_QEx-HFX_DSS_R1.csv"
MEROX2 = "data/merox/XLpeplib_Beveridge_QEx-HFX_DSS_R1.zhrm"
MEROX3 = "data/merox/cas11.csv"


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


@pytest.mark.parametrize(
    "retention_time,expected",
    [
        (3099, 3099.0),
        ("1748", 1748.0),
        ("1904.59542", 1904.59542),
        # range of merged scans, e.g. in PRIDE deposits PXD008348 and PXD019215
        ("1854-1886", 1854.0),
        ("2241-2241", 2241.0),
        # empty in exported tables, e.g. in PRIDE deposits PXD012127 and PXD017299
        (" ", None),
    ],
)
def test8(retention_time, expected):
    from pyXLMS.parser._parser_xldbse_merox import __get_merox_retention_time

    if expected is None:
        assert __get_merox_retention_time(retention_time) is None
    else:
        assert __get_merox_retention_time(retention_time) == pytest.approx(expected)


def test9():
    from pyXLMS.parser._parser_xldbse_merox import __get_merox_retention_time

    with pytest.raises(RuntimeError, match="Could not parse retention time from"):
        _ = __get_merox_retention_time("1854-")


@pytest.mark.parametrize(
    "retention_time,expected",
    [
        ("1854,5-1886,5", 1854.5),
        ("1904,5", 1904.5),
        (1904.5, 1904.5),
    ],
)
def test10(retention_time, expected):
    # retention times with decimal comma
    from pyXLMS.parser._parser_xldbse_merox import __get_merox_retention_time

    assert __get_merox_retention_time(retention_time, ",") == pytest.approx(expected)


def test11(tmp_path):
    # pandas does not apply the decimal comma to a column with retention time ranges
    import pandas as pd
    from pyXLMS.parser import read_merox

    data = pd.read_csv(MEROX1, sep=";").head(2)
    data["Retention time in sec"] = ["1854,5-1886,5", "1904,5"]
    file = tmp_path / "decimal_comma.csv"
    data.to_csv(file, sep=";", index=False, decimal=",")
    pr = read_merox(str(file), crosslinker="DSS", decimal=",")
    csms = pr["crosslink-spectrum-matches"]
    assert [csm["retention_time"] for csm in csms] == pytest.approx([1854.5, 1904.5])
