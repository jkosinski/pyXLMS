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
# excerpt of public PRIDE deposit PXD044723 (rows copied verbatim)
MEROX8 = "data/merox/PXD044723_XLDID2_01_excerpt.zhrm"
# excerpt of public PRIDE deposit PXD010070 (exported table, rows copied verbatim)
MEROX6 = "data/merox/PXD010070_KMN-C18_F1_excerpt.csv"


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
    with pytest.raises(RuntimeError, match="Could not scan number from"):
        _ = __get_merox_scan_number('Locus:1.1.1.1409.10 File:"run.wiff"')


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
def test12(retention_time, expected):
    from pyXLMS.parser._parser_xldbse_merox import __get_merox_retention_time

    if expected is None:
        assert __get_merox_retention_time(retention_time) is None
    else:
        assert __get_merox_retention_time(retention_time) == pytest.approx(expected)


def test13():
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
def test19(retention_time, expected):
    # retention times with decimal comma
    from pyXLMS.parser._parser_xldbse_merox import __get_merox_retention_time

    assert __get_merox_retention_time(retention_time, ",") == pytest.approx(expected)


def test20(tmp_path):
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


def test16():
    # 3 crosslinks, 2 monolinks (Peptide2 "0") and 2 looplinks (Peptide2 "1")
    from pyXLMS.parser import read_merox

    with pytest.warns(
        RuntimeWarning,
        match=r"Skipped 2 monolink and 2 looplink rows, only crosslinks are read!",
    ):
        pr = read_merox(MEROX8, crosslinker="DSBU")
    csms = pr["crosslink-spectrum-matches"]
    assert len(csms) == 3
    assert [csm["scan_nr"] for csm in csms] == [1581, 1593, 2978]
    assert [csm["beta_peptide"] for csm in csms] == ["LIKSR", "LIKSR", "LIKSR"]
    assert [csm["beta_proteins_crosslink_positions"] for csm in csms] == [
        [364],
        [365],
        [364],
    ]


def test17():
    # no warning if there are no monolinks or looplinks (warnings are errors in pytest.ini)
    from pyXLMS.parser import read_merox

    pr = read_merox(MEROX1, crosslinker="DSS")
    assert len(pr["crosslink-spectrum-matches"]) == 93


def test18():
    # exported table of an older MeroX version: rows end with a separator, no peptide scores
    from pyXLMS.parser import read_merox

    pr = read_merox(MEROX6, crosslinker="DSBU")
    csms = pr["crosslink-spectrum-matches"]
    assert len(csms) == 2
    csm = csms[0]
    assert csm["alpha_peptide"] == "KWYPLRPK"
    assert csm["alpha_proteins"] == ["Nsl1"]
    assert csm["alpha_proteins_crosslink_positions"] == [268]
    assert csm["alpha_score"] is None
    assert csm["beta_peptide"] == "LLGEEIEYLKR"
    assert csm["beta_proteins_crosslink_positions"] == [227]
    assert csm["beta_score"] is None
    assert csm["score"] == pytest.approx(91.0)
    assert csm["scan_nr"] == 10001
    assert csm["retention_time"] == pytest.approx(3379.0815)
    assert csm["charge"] == 4
