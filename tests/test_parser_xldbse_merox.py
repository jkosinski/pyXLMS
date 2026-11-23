#!/usr/bin/env python3

# pyXLMS - TESTS
# 2025 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

import pytest

MEROX1 = "data/merox/XLpeplib_Beveridge_QEx-HFX_DSS_R1.csv"
MEROX2 = "data/merox/XLpeplib_Beveridge_QEx-HFX_DSS_R1.zhrm"
MEROX3 = "data/merox/cas11.csv"
# excerpts of public PRIDE deposits PXD008348 and PXD015080 (rows copied verbatim)
MEROX4 = "data/merox/PXD008348_Npl3HMT1_DSSO_R1_excerpt.zhrm"
MEROX5 = "data/merox/PXD015080_N29_excerpt.zhrm"
# excerpt of public PRIDE deposit PXD025064 (rows copied verbatim)
MEROX7 = "data/merox/PXD025064_AGAL_300DSBU_Fr7_excerpt.zhrm"
# excerpt of public PRIDE deposit PXD033205 (rows copied verbatim)
MEROX9 = "data/merox/PXD033205_Day0_DSBU_1_excerpt.zhrm"


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
    # scan number without spectrum file, e.g. "8464"
    from pyXLMS.parser import read_merox

    pr = read_merox(MEROX4, crosslinker="DSSO")
    csms = pr["crosslink-spectrum-matches"]
    assert len(csms) == 3
    assert [csm["spectrum_file"] for csm in csms] == ["", "", ""]
    assert [csm["scan_nr"] for csm in csms] == [8464, 9519, 10115]
    assert [csm["retention_time"] for csm in csms] == pytest.approx(
        [1875.0, 2081.0, 2192.0]
    )
    assert csms[0]["alpha_peptide"] == "DSATEKTK"
    assert csms[0]["beta_peptide"] == "LTAKR"


def test9():
    # spectrum title with raw file and native id, with and without leading "(-1)"
    from pyXLMS.parser import read_merox

    pr = read_merox(MEROX5, crosslinker="HSAKp", parse_modifications=False)
    csms = pr["crosslink-spectrum-matches"]
    assert len(csms) == 4
    assert [csm["spectrum_file"] for csm in csms] == [
        "20190726_stamch_BSA_HSAKp_MODDE_N29"
    ] * 4
    assert [csm["scan_nr"] for csm in csms] == [10031, 10529, 10498, 10539]
    assert [csm["charge"] for csm in csms] == [4, 5, 4, 4]


@pytest.mark.parametrize(
    "scan_nr_and_file,scan_nr,spectrum_file",
    [
        (
            "10061~XLpeplib_Beveridge_QEx-HFX_DSS_R1",
            10061,
            "XLpeplib_Beveridge_QEx-HFX_DSS_R1",
        ),
        ("7755", 7755, ""),
        ("9512-9606", 9512, ""),
        (
            "10005-10006~f3-1828-MS2MS2-rhogrk-mc4-01",
            10005,
            "f3-1828-MS2MS2-rhogrk-mc4-01",
        ),
        (
            "(-1) 10263-10264~f3-1828-MS2MS2-rhogrk-mc4-01",
            10263,
            "f3-1828-MS2MS2-rhogrk-mc4-01",
        ),
        (
            '(-1) 20190726_stamch_BSA_HSAKp_MODDE_N29.10031.10031.4 File:"20190726_stamch_BSA_HSAKp_MODDE_N29.raw", '
            'NativeID:"controllerType=0 controllerNumber=1 scan=10031"',
            10031,
            "20190726_stamch_BSA_HSAKp_MODDE_N29",
        ),
        (
            'X-link_trypsin_Craig.3721.3721.2 File:"X-link_trypsin_Craig.raw", '
            'NativeID:"controllerType=0 controllerNumber=1 scan=3721"',
            3721,
            "X-link_trypsin_Craig",
        ),
        (
            "NS=2247,2252,2254-2257~APX-BA5_100DSBU_S1-C1_1_4266.d_unsorted",
            2247,
            "APX-BA5_100DSBU_S1-C1_1_4266.d_unsorted",
        ),
        (
            "NS=sn10318~6b705f9e9cdedd27___Xpl1_001268.raw.-1",
            10318,
            "6b705f9e9cdedd27___Xpl1_001268.raw.-1",
        ),
        (
            'Locus:1.1.1.2950.17 File:"TB05_M1.wiff"',
            None,
            'Locus:1.1.1.2950.17 File:"TB05_M1.wiff"',
        ),
        # spectrum title with annotations (PRIDE deposit PXD031339)
        (
            'QE35501LM.1633.1633.4 MSn:"707.854395877474" RTmin:"16.776091999999998" '
            'File:"QE35501LM.raw", NativeID:"controllerType=0 controllerNumber=1 scan=1633"',
            1633,
            "QE35501LM",
        ),
        # the same scan number written with leading zeros
        (
            'run.001633.001633.4 File:"run.raw", NativeID:"controllerType=0 scan=1633"',
            1633,
            "run",
        ),
        # StavroX (PRIDE deposit PXD028039)
        ("3324_AJ_Dimer-Tryp|3953", 3953, "3324_AJ_Dimer-Tryp"),
        # Bruker timsTOF (PRIDE deposit PXD027234), no scan number
        (
            "(-1) 11238 +MS2(650.0308), 32.6eV, 52.2min, 1/K0=0.850 #17538-17540|"
            "count(pasefframemsmsinfo.precursor)=3|delta=25|Id=11238|AverageMz=650.4381726071|"
            "MonoisotopicMz=650.030798142347|Charge=3|Intensity=7751|"
            "ScanNumber=1161.06182295069~BSA_DSAU125_1",
            None,
            "11238 +MS2(650.0308), 32.6eV, 52.2min, 1/K0=0.850 #17538-17540|"
            "count(pasefframemsmsinfo.precursor)=3|delta=25|Id=11238|AverageMz=650.4381726071|"
            "MonoisotopicMz=650.030798142347|Charge=3|Intensity=7751|"
            "ScanNumber=1161.06182295069~BSA_DSAU125_1",
        ),
        (
            "10067, 10078 +MS2(647.7008), 35.6eV, 53.6min, 1/K0=0.858 #17593-17607|"
            "Id=10067|Charge=3~Ribosome_DSBU_1",
            None,
            "10067, 10078 +MS2(647.7008), 35.6eV, 53.6min, 1/K0=0.858 #17593-17607|"
            "Id=10067|Charge=3~Ribosome_DSBU_1",
        ),
        # no spectrum title (PRIDE deposit PXD021872)
        ("N/A", None, ""),
    ],
)
def test10(scan_nr_and_file, scan_nr, spectrum_file):
    from pyXLMS.parser._parser_xldbse_merox import __get_merox_scan_number
    from pyXLMS.parser._parser_xldbse_merox import __get_merox_spectrum_file

    assert __get_merox_scan_number(scan_nr_and_file) == scan_nr
    assert __get_merox_spectrum_file(scan_nr_and_file) == spectrum_file


def test11():
    from pyXLMS.parser._parser_xldbse_merox import __get_merox_scan_number

    with pytest.raises(RuntimeError, match="Could not scan number from"):
        _ = __get_merox_scan_number('scan File:"run.raw"')


def test_merox_sciex_spectrum_titles():
    # SCIEX spectrum titles (PRIDE deposit PXD024519) contain no scan number, the whole
    # title identifies the spectrum
    import io
    from pyXLMS.parser import read_merox
    from pyXLMS.transform import unique

    titles = [
        'Locus:1.1.1.2950.17 File:"TB05_M1.wiff"',
        'Locus:1.1.1.2950.17 File:"TB05_M1.wiff"',
        'Locus:1.1.1.2951.18 File:"TB05_M1.wiff"',
    ]
    with open(MEROX1, "r", encoding="utf-8") as f:
        lines = f.read().splitlines()
    rows = [lines[0]]
    for line, title in zip(lines[1:4], titles):
        cells = line.split(";")
        cells[14] = title
        rows.append(";".join(cells))
    pr = read_merox(io.BytesIO("\n".join(rows).encode("utf-8")), crosslinker="DSS")
    csms = pr["crosslink-spectrum-matches"]
    assert len(csms) == 3
    assert [csm["spectrum_file"] for csm in csms] == titles
    assert [csm["scan_nr"] for csm in csms] == [None, None, None]
    assert [csm["completeness"] for csm in csms] == ["partial"] * 3
    best = unique(csms)
    assert len(best) == 2
    assert [csm["spectrum_file"] for csm in best] == [titles[0], titles[2]]
    assert best[0]["score"] == pytest.approx(151)


def test_merox_scan_number_of_title_and_native_id_differ():
    from pyXLMS.parser._parser_xldbse_merox import __get_merox_scan_number

    with pytest.raises(RuntimeError, match="Could not scan number from"):
        _ = __get_merox_scan_number(
            'run.1633.1633.4 File:"run.raw", NativeID:"controllerType=0 scan=1634"'
        )


def test14():
    # protein-terminal peptides: "{" protein N-terminus, "}" protein C-terminus
    from pyXLMS.parser import read_merox
    from pyXLMS.transform import modifications_to_str as mts

    pr = read_merox(MEROX7, crosslinker="DSBU")
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
def test15(position_str, from_position, sequence, xl_position, pep_position):
    from pyXLMS.parser._parser_xldbse_merox import __get_merox_crosslink_position
    from pyXLMS.parser._parser_xldbse_merox import __get_merox_peptide_position

    assert __get_merox_crosslink_position(position_str) == xl_position
    assert __get_merox_peptide_position(from_position, sequence) == pep_position


# the letters of PXD033205's properties.mxf: "Acetylation N-Term;X;C2OH2" with the static
# modification "{;X", and "Methionine oxidized;J;C5H9NO2S" with the variable modification "M;J;2"
MEROX9_MODIFICATIONS = {
    "J": {"Amino Acid": "M", "Modification": ("Methionine oxidized", 15.994915)},
    "X": {"Amino Acid": "{", "Modification": ("Acetylation N-Term", 42.010565)},
}


def test_merox_protein_n_terminal_group():
    # "[XMDVFMKGLSK]" From 0: X is the acetylated protein N-terminus, counted as position 1
    from pyXLMS.parser import read_merox
    from pyXLMS.constants import MEROX_MODIFICATION_MAPPING
    from pyXLMS.transform import modifications_to_str as mts

    pr = read_merox(
        MEROX9,
        crosslinker="DSBU",
        modifications={**MEROX_MODIFICATION_MAPPING, **MEROX9_MODIFICATIONS},
    )
    csms = pr["crosslink-spectrum-matches"]
    assert len(csms) == 3

    # [AKEGVVAAAEKTK] K11 - [XMDVFMKGLSK] K7
    csm = csms[0]
    assert csm["alpha_peptide"] == "AKEGVVAAAEKTK"
    assert csm["alpha_proteins_crosslink_positions"] == [21]
    assert csm["beta_peptide"] == "MDVFMKGLSK"
    assert csm["beta_peptide_crosslink_position"] == 6
    assert (
        mts(csm["beta_modifications"])
        == "(0:[Acetylation N-Term|42.010565]);(6:[DSBU|196.08479231])"
    )
    assert csm["beta_proteins"] == ["sp|P37840|SYUA_N14"]
    assert csm["beta_proteins_crosslink_positions"] == [6]
    assert csm["beta_proteins_peptide_positions"] == [1]

    # [XMDVFMKGLSK] K7 - [AAGKTKE] K4
    csm = csms[1]
    assert csm["alpha_peptide"] == "AAGKTKE"
    assert csm["alpha_proteins_crosslink_positions"] == [32]
    assert csm["beta_peptide"] == "MDVFMKGLSK"
    assert csm["beta_peptide_crosslink_position"] == 6
    assert csm["beta_proteins_crosslink_positions"] == [6]
    assert csm["beta_proteins_peptide_positions"] == [1]

    # [XMDVFJKGLSK] K7 - [KTKQGVAE] K3
    csm = csms[2]
    assert csm["beta_peptide"] == "MDVFMKGLSK"
    assert csm["beta_peptide_crosslink_position"] == 6
    assert mts(csm["beta_modifications"]) == (
        "(0:[Acetylation N-Term|42.010565]);(5:[Methionine oxidized|15.994915]);"
        "(6:[DSBU|196.08479231])"
    )
    assert csm["beta_proteins_crosslink_positions"] == [6]
    assert csm["beta_proteins_peptide_positions"] == [1]


def test_merox_protein_n_terminal_group_without_modifications():
    from pyXLMS.parser import read_merox
    from pyXLMS.constants import MEROX_MODIFICATION_MAPPING

    pr = read_merox(
        MEROX9,
        crosslinker="DSBU",
        parse_modifications=False,
        modifications={**MEROX_MODIFICATION_MAPPING, **MEROX9_MODIFICATIONS},
    )
    csm = pr["crosslink-spectrum-matches"][2]
    assert csm["beta_peptide"] == "MDVFMKGLSK"
    assert csm["beta_modifications"] is None
    assert csm["beta_peptide_crosslink_position"] == 6
    assert csm["beta_proteins_crosslink_positions"] == [6]
    assert csm["beta_proteins_peptide_positions"] == [1]


def test_merox_protein_n_terminal_group_not_mapped():
    # without the letter of the N-terminus, e.g. in an exported table without the settings
    from pyXLMS.parser import read_merox
    from pyXLMS.constants import MEROX_MODIFICATION_MAPPING

    with pytest.raises(ValueError, match=r"MeroX peptide \[XMDVFMKGLSK\] starts at"):
        _ = read_merox(
            MEROX9,
            crosslinker="DSBU",
            modifications={
                **MEROX_MODIFICATION_MAPPING,
                "J": MEROX9_MODIFICATIONS["J"],
            },
        )


@pytest.mark.parametrize(
    "sequence,position_str,expected_sequence,group,xl_position",
    [
        pytest.param(
            "[XMDVFMKGLSK]",
            "K7",
            "{MDVFMKGLSK]",
            ("Acetylation N-Term", 42.010565),
            6,
            id="group",
        ),
        pytest.param("[KTKQGVAE]", "K3", "[KTKQGVAE]", None, 3, id="no-group"),
        pytest.param("{MSKTAVK]", "{0", "{MSKTAVK]", None, 1, id="protein-n-term"),
        pytest.param("[AXK]", "K3", "[AXK]", None, 3, id="letter-not-first"),
    ],
)
def test_merox_protein_n_terminal_group_helpers(
    sequence, position_str, expected_sequence, group, xl_position
):
    from pyXLMS.parser._parser_xldbse_merox import (
        __get_merox_protein_n_terminal_group,
        __get_merox_group_crosslink_position,
    )

    parsed_sequence, parsed_group = __get_merox_protein_n_terminal_group(
        sequence, MEROX9_MODIFICATIONS
    )
    assert parsed_sequence == expected_sequence
    assert parsed_group == group
    assert (
        __get_merox_group_crosslink_position(position_str, parsed_group) == xl_position
    )


@pytest.mark.parametrize("position_str", ["X1", "{0"])
def test_merox_crosslink_on_protein_n_terminal_group(position_str):
    from pyXLMS.parser._parser_xldbse_merox import __get_merox_group_crosslink_position

    with pytest.raises(ValueError, match="is on the modified protein N-terminus"):
        __get_merox_group_crosslink_position(
            position_str, ("Acetylation N-Term", 42.010565)
        )
