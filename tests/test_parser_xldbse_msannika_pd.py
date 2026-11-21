#!/usr/bin/env python3

# pyXLMS - TESTS
# 2024 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

import pytest


PD = "data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1.pdResult"
# excerpt of public PRIDE deposit PXD056680 (Proteome Discoverer 2.4, rows copied verbatim)
PD_PXD056680 = "data/ms_annika/PXD056680_F008604_excerpt.pdResult"
# the same rows of the .msf file of the same search
PD_PXD056680_MSF = "data/ms_annika/PXD056680_F008604_excerpt.msf"


def test1():
    from pyXLMS import parser as p
    from pyXLMS.transform import modifications_to_str as mts

    parser_result = p.read_msannika(PD)
    assert parser_result["data_type"] == "parser_result"
    assert parser_result["completeness"] == "full"
    assert parser_result["search_engine"] == "MS Annika"
    assert parser_result["crosslink-spectrum-matches"] is not None
    assert parser_result["crosslinks"] is not None

    csms = parser_result["crosslink-spectrum-matches"]
    assert len(csms) == 826

    crosslinks = parser_result["crosslinks"]
    assert len(crosslinks) == 300

    csm = csms[822]
    assert csm["data_type"] == "crosslink-spectrum-match"
    assert csm["completeness"] == "full"
    assert csm["alpha_peptide"] == "KIECFDSVEISGVEDR"
    assert (
        mts(csm["alpha_modifications"])
        == "(1:[DSS|138.06808]);(4:[Carbamidomethyl|57.021464])"
    )
    assert csm["alpha_peptide_crosslink_position"] == 1
    assert csm["alpha_proteins"] == ["Cas9"]
    assert csm["alpha_proteins_crosslink_positions"] == [575]
    assert csm["alpha_proteins_peptide_positions"] == [575]
    assert csm["alpha_score"] == pytest.approx(376.15, abs=0.01)
    assert not csm["alpha_decoy"]
    assert csm["beta_peptide"] == "TILDFLKSDGFANR"
    assert mts(csm["beta_modifications"]) == "(7:[DSS|138.06808])"
    assert csm["beta_peptide_crosslink_position"] == 7
    assert csm["beta_proteins"] == ["Cas9"]
    assert csm["beta_proteins_crosslink_positions"] == [688]
    assert csm["beta_proteins_peptide_positions"] == [682]
    assert csm["beta_score"] == pytest.approx(393.87, abs=0.01)
    assert not csm["beta_decoy"]
    assert csm["crosslink_type"] == "intra"
    assert csm["score"] == pytest.approx(376.15, abs=0.01)
    assert csm["spectrum_file"] == "XLpeplib_Beveridge_QEx-HFX_DSS_R1.raw"
    assert csm["scan_nr"] == 23454
    assert csm["charge"] == 4
    assert csm["retention_time"] == pytest.approx(6278.5199999999995, abs=0.01)
    assert csm["ion_mobility"] == pytest.approx(0.0, abs=0.01)

    xl = crosslinks[0]
    assert xl["data_type"] == "crosslink"
    assert xl["completeness"] == "full"
    assert xl["alpha_peptide"] == "GQKNSR"
    assert xl["alpha_peptide_crosslink_position"] == 3
    assert xl["alpha_proteins"] == ["Cas9"]
    assert xl["alpha_proteins_crosslink_positions"] == [779]
    assert not xl["alpha_decoy"]
    assert xl["beta_peptide"] == "GQKNSR"
    assert xl["beta_peptide_crosslink_position"] == 3
    assert xl["beta_proteins"] == ["Cas9"]
    assert xl["beta_proteins_crosslink_positions"] == [779]
    assert not xl["beta_decoy"]
    assert xl["crosslink_type"] == "intra"
    assert xl["score"] == pytest.approx(119.83, abs=0.01)


def test2():
    # pdResult without compensation voltage column, decoy peptides without protein
    from pyXLMS.parser import read_msannika

    pr = read_msannika(PD_PXD056680)
    csms = pr["crosslink-spectrum-matches"]
    assert len(csms) == 5
    assert csms[0]["alpha_proteins"] == ["AT3G50590_TWD40-1_araport11"]
    assert csms[0]["alpha_proteins_crosslink_positions"] == [1494]
    assert csms[0]["ion_mobility"] is None
    csm = csms[3]
    assert csm["alpha_peptide"] == "KLAER"
    assert csm["alpha_peptide_crosslink_position"] == 1
    assert csm["alpha_proteins"] is None
    assert csm["alpha_proteins_crosslink_positions"] is None
    assert csm["alpha_decoy"]
    assert csm["beta_proteins"] == ["recTPLATE_after_TEV_cleavage"]
    assert csm["beta_proteins_crosslink_positions"] == [1174]
    assert not csm["beta_decoy"]

    xls = pr["crosslinks"]
    assert len(xls) == 5
    xl = xls[4]
    assert xl["alpha_peptide"] == "ELVPMKR"
    assert xl["alpha_peptide_crosslink_position"] == 6
    assert xl["alpha_proteins"] is None
    assert xl["beta_proteins_crosslink_positions"] == [894]
    assert xl["alpha_decoy"]
    assert not xl["beta_decoy"]


def test3():
    # .msf file with the same tables as the .pdResult file
    from pyXLMS.parser import read_msannika

    pr = read_msannika(PD_PXD056680_MSF)
    assert pr["search_engine"] == "MS Annika"
    pr_pdresult = read_msannika(PD_PXD056680)
    for key in ["crosslink-spectrum-matches", "crosslinks"]:
        assert len(pr[key]) == 5
        for match, match_pdresult in zip(pr[key], pr_pdresult[key]):
            for attribute in [
                "alpha_peptide",
                "alpha_peptide_crosslink_position",
                "alpha_proteins",
                "alpha_proteins_crosslink_positions",
                "alpha_decoy",
                "beta_peptide",
                "beta_peptide_crosslink_position",
                "beta_proteins",
                "beta_proteins_crosslink_positions",
                "beta_decoy",
                "crosslink_type",
                "score",
            ]:
                assert match[attribute] == match_pdresult[attribute]
    csm = pr["crosslink-spectrum-matches"][3]
    assert csm["alpha_peptide"] == "KLAER"
    assert csm["alpha_proteins"] is None
    assert csm["alpha_decoy"]
    assert csm["scan_nr"] == 17383
