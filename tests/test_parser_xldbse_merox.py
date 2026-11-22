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
