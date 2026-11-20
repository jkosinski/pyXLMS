#!/usr/bin/env python3

# pyXLMS - TESTS
# 2024 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

import pytest


XL_TSV = "data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1_Crosslinks.txt"
XL_XLSX = "data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1_Crosslinks.xlsx"
CSMS_TSV = "data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1_CSMs.txt"
CSMS_XLSX = "data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1_CSMs.xlsx"
# excerpts of public PRIDE deposit PXD043531 (unfiltered, rows copied verbatim)
CSMS_TSV_PXD043531 = "data/ms_annika/PXD043531_DSBSO_HEK_MSAnnika_CSMs_excerpt.txt"
XL_TSV_PXD043531 = "data/ms_annika/PXD043531_DSBSO_HEK_MSAnnika_Crosslinks_excerpt.txt"


def test1():
    from pyXLMS import parser as p

    parser_result = p.read_msannika(XL_TSV)
    assert parser_result["data_type"] == "parser_result"
    assert parser_result["completeness"] == "partial"
    assert parser_result["search_engine"] == "MS Annika"
    assert parser_result["crosslink-spectrum-matches"] is None
    assert parser_result["crosslinks"] is not None

    crosslinks = parser_result["crosslinks"]
    assert len(crosslinks) == 300

    first_crosslink = crosslinks[0]
    last_crosslink = crosslinks[299]

    assert first_crosslink["data_type"] == "crosslink"
    assert first_crosslink["completeness"] == "full"
    assert first_crosslink["alpha_peptide"] == "GQKNSR"
    assert first_crosslink["alpha_peptide_crosslink_position"] == 3
    assert first_crosslink["alpha_proteins"] == ["Cas9"]
    assert first_crosslink["alpha_proteins_crosslink_positions"] == [779]
    assert not first_crosslink["alpha_decoy"]
    assert first_crosslink["beta_peptide"] == "GQKNSR"
    assert first_crosslink["beta_peptide_crosslink_position"] == 3
    assert first_crosslink["beta_proteins"] == ["Cas9"]
    assert first_crosslink["beta_proteins_crosslink_positions"] == [779]
    assert not first_crosslink["beta_decoy"]
    assert first_crosslink["crosslink_type"] == "intra"
    assert first_crosslink["score"] == pytest.approx(119.83)

    assert last_crosslink["data_type"] == "crosslink"
    assert last_crosslink["completeness"] == "full"
    assert last_crosslink["alpha_peptide"] == "MEDESKLHKFKDFK"
    assert last_crosslink["alpha_peptide_crosslink_position"] == 11
    assert last_crosslink["alpha_proteins"] == ["sp"]
    assert last_crosslink["alpha_proteins_crosslink_positions"] == [109]
    assert last_crosslink["alpha_decoy"]
    assert last_crosslink["beta_peptide"] == "SSFEKNPIDFLEAK"
    assert last_crosslink["beta_peptide_crosslink_position"] == 5
    assert last_crosslink["beta_proteins"] == ["Cas9"]
    assert last_crosslink["beta_proteins_crosslink_positions"] == [1180]
    assert last_crosslink["beta_decoy"]
    assert last_crosslink["crosslink_type"] == "inter"
    assert last_crosslink["score"] == pytest.approx(15.89)


def test2():
    from pyXLMS import parser as p

    parser_result = p.read_msannika(XL_XLSX)
    assert parser_result["data_type"] == "parser_result"
    assert parser_result["completeness"] == "partial"
    assert parser_result["search_engine"] == "MS Annika"
    assert parser_result["crosslink-spectrum-matches"] is None
    assert parser_result["crosslinks"] is not None

    crosslinks = parser_result["crosslinks"]
    assert len(crosslinks) == 300

    first_crosslink = crosslinks[0]
    last_crosslink = crosslinks[299]

    assert first_crosslink["data_type"] == "crosslink"
    assert first_crosslink["completeness"] == "full"
    assert first_crosslink["alpha_peptide"] == "GQKNSR"
    assert first_crosslink["alpha_peptide_crosslink_position"] == 3
    assert first_crosslink["alpha_proteins"] == ["Cas9"]
    assert first_crosslink["alpha_proteins_crosslink_positions"] == [779]
    assert not first_crosslink["alpha_decoy"]
    assert first_crosslink["beta_peptide"] == "GQKNSR"
    assert first_crosslink["beta_peptide_crosslink_position"] == 3
    assert first_crosslink["beta_proteins"] == ["Cas9"]
    assert first_crosslink["beta_proteins_crosslink_positions"] == [779]
    assert not first_crosslink["beta_decoy"]
    assert first_crosslink["crosslink_type"] == "intra"
    assert first_crosslink["score"] == pytest.approx(119.83)

    assert last_crosslink["data_type"] == "crosslink"
    assert last_crosslink["completeness"] == "full"
    assert last_crosslink["alpha_peptide"] == "MEDESKLHKFKDFK"
    assert last_crosslink["alpha_peptide_crosslink_position"] == 11
    assert last_crosslink["alpha_proteins"] == ["sp"]
    assert last_crosslink["alpha_proteins_crosslink_positions"] == [109]
    assert last_crosslink["alpha_decoy"]
    assert last_crosslink["beta_peptide"] == "SSFEKNPIDFLEAK"
    assert last_crosslink["beta_peptide_crosslink_position"] == 5
    assert last_crosslink["beta_proteins"] == ["Cas9"]
    assert last_crosslink["beta_proteins_crosslink_positions"] == [1180]
    assert last_crosslink["beta_decoy"]
    assert last_crosslink["crosslink_type"] == "inter"
    assert last_crosslink["score"] == pytest.approx(15.89)


def test3():
    from pyXLMS import parser as p

    parser_result = p.read_msannika([XL_TSV, XL_XLSX])
    assert parser_result["data_type"] == "parser_result"
    assert parser_result["completeness"] == "partial"
    assert parser_result["search_engine"] == "MS Annika"
    assert parser_result["crosslink-spectrum-matches"] is None
    assert parser_result["crosslinks"] is not None

    crosslinks = parser_result["crosslinks"]
    assert len(crosslinks) == 600

    first_crosslink = crosslinks[0]
    last_crosslink = crosslinks[599]

    assert first_crosslink["data_type"] == "crosslink"
    assert first_crosslink["completeness"] == "full"
    assert first_crosslink["alpha_peptide"] == "GQKNSR"
    assert first_crosslink["alpha_peptide_crosslink_position"] == 3
    assert first_crosslink["alpha_proteins"] == ["Cas9"]
    assert first_crosslink["alpha_proteins_crosslink_positions"] == [779]
    assert not first_crosslink["alpha_decoy"]
    assert first_crosslink["beta_peptide"] == "GQKNSR"
    assert first_crosslink["beta_peptide_crosslink_position"] == 3
    assert first_crosslink["beta_proteins"] == ["Cas9"]
    assert first_crosslink["beta_proteins_crosslink_positions"] == [779]
    assert not first_crosslink["beta_decoy"]
    assert first_crosslink["crosslink_type"] == "intra"
    assert first_crosslink["score"] == pytest.approx(119.83)

    assert last_crosslink["data_type"] == "crosslink"
    assert last_crosslink["completeness"] == "full"
    assert last_crosslink["alpha_peptide"] == "MEDESKLHKFKDFK"
    assert last_crosslink["alpha_peptide_crosslink_position"] == 11
    assert last_crosslink["alpha_proteins"] == ["sp"]
    assert last_crosslink["alpha_proteins_crosslink_positions"] == [109]
    assert last_crosslink["alpha_decoy"]
    assert last_crosslink["beta_peptide"] == "SSFEKNPIDFLEAK"
    assert last_crosslink["beta_peptide_crosslink_position"] == 5
    assert last_crosslink["beta_proteins"] == ["Cas9"]
    assert last_crosslink["beta_proteins_crosslink_positions"] == [1180]
    assert last_crosslink["beta_decoy"]
    assert last_crosslink["crosslink_type"] == "inter"
    assert last_crosslink["score"] == pytest.approx(15.89)


def test4():
    from pyXLMS import parser as p

    parser_result = p.read_msannika([XL_TSV, XL_TSV])
    assert parser_result["data_type"] == "parser_result"
    assert parser_result["completeness"] == "partial"
    assert parser_result["search_engine"] == "MS Annika"
    assert parser_result["crosslink-spectrum-matches"] is None
    assert parser_result["crosslinks"] is not None

    crosslinks = parser_result["crosslinks"]
    assert len(crosslinks) == 600

    first_crosslink = crosslinks[0]
    last_crosslink = crosslinks[599]

    assert first_crosslink["data_type"] == "crosslink"
    assert first_crosslink["completeness"] == "full"
    assert first_crosslink["alpha_peptide"] == "GQKNSR"
    assert first_crosslink["alpha_peptide_crosslink_position"] == 3
    assert first_crosslink["alpha_proteins"] == ["Cas9"]
    assert first_crosslink["alpha_proteins_crosslink_positions"] == [779]
    assert not first_crosslink["alpha_decoy"]
    assert first_crosslink["beta_peptide"] == "GQKNSR"
    assert first_crosslink["beta_peptide_crosslink_position"] == 3
    assert first_crosslink["beta_proteins"] == ["Cas9"]
    assert first_crosslink["beta_proteins_crosslink_positions"] == [779]
    assert not first_crosslink["beta_decoy"]
    assert first_crosslink["crosslink_type"] == "intra"
    assert first_crosslink["score"] == pytest.approx(119.83)

    assert last_crosslink["data_type"] == "crosslink"
    assert last_crosslink["completeness"] == "full"
    assert last_crosslink["alpha_peptide"] == "MEDESKLHKFKDFK"
    assert last_crosslink["alpha_peptide_crosslink_position"] == 11
    assert last_crosslink["alpha_proteins"] == ["sp"]
    assert last_crosslink["alpha_proteins_crosslink_positions"] == [109]
    assert last_crosslink["alpha_decoy"]
    assert last_crosslink["beta_peptide"] == "SSFEKNPIDFLEAK"
    assert last_crosslink["beta_peptide_crosslink_position"] == 5
    assert last_crosslink["beta_proteins"] == ["Cas9"]
    assert last_crosslink["beta_proteins_crosslink_positions"] == [1180]
    assert last_crosslink["beta_decoy"]
    assert last_crosslink["crosslink_type"] == "inter"
    assert last_crosslink["score"] == pytest.approx(15.89)


def test5():
    from pyXLMS import parser as p
    from pyXLMS.transform import modifications_to_str as mts

    parser_result = p.read_msannika(CSMS_TSV)
    assert parser_result["data_type"] == "parser_result"
    assert parser_result["completeness"] == "partial"
    assert parser_result["search_engine"] == "MS Annika"
    assert parser_result["crosslink-spectrum-matches"] is not None
    assert parser_result["crosslinks"] is None

    csms = parser_result["crosslink-spectrum-matches"]
    assert len(csms) == 826

    first_csm = csms[0]
    csm = csms[822]
    last_csm = csms[825]

    assert first_csm["data_type"] == "crosslink-spectrum-match"
    assert first_csm["completeness"] == "full"
    assert first_csm["alpha_peptide"] == "GQKNSR"
    assert mts(first_csm["alpha_modifications"]) == "(3:[DSS|138.06808])"
    assert first_csm["alpha_peptide_crosslink_position"] == 3
    assert first_csm["alpha_proteins"] == ["Cas9"]
    assert first_csm["alpha_proteins_crosslink_positions"] == [779]
    assert first_csm["alpha_proteins_peptide_positions"] == [777]
    assert first_csm["alpha_score"] == pytest.approx(119.83)
    assert not first_csm["alpha_decoy"]
    assert first_csm["beta_peptide"] == "GQKNSR"
    assert mts(first_csm["beta_modifications"]) == "(3:[DSS|138.06808])"
    assert first_csm["beta_peptide_crosslink_position"] == 3
    assert first_csm["beta_proteins"] == ["Cas9"]
    assert first_csm["beta_proteins_crosslink_positions"] == [779]
    assert first_csm["beta_proteins_peptide_positions"] == [777]
    assert first_csm["beta_score"] == pytest.approx(119.83)
    assert not first_csm["beta_decoy"]
    assert first_csm["crosslink_type"] == "intra"
    assert first_csm["score"] == pytest.approx(119.83)
    assert first_csm["spectrum_file"] == "XLpeplib_Beveridge_QEx-HFX_DSS_R1.raw"
    assert first_csm["scan_nr"] == 2257
    assert first_csm["charge"] == 3
    assert first_csm["retention_time"] == pytest.approx(733.188)
    assert first_csm["ion_mobility"] == pytest.approx(0.0)

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
    assert csm["alpha_score"] == pytest.approx(376.15)
    assert not csm["alpha_decoy"]
    assert csm["beta_peptide"] == "TILDFLKSDGFANR"
    assert mts(csm["beta_modifications"]) == "(7:[DSS|138.06808])"
    assert csm["beta_peptide_crosslink_position"] == 7
    assert csm["beta_proteins"] == ["Cas9"]
    assert csm["beta_proteins_crosslink_positions"] == [688]
    assert csm["beta_proteins_peptide_positions"] == [682]
    assert csm["beta_score"] == pytest.approx(393.87)
    assert not csm["beta_decoy"]
    assert csm["crosslink_type"] == "intra"
    assert csm["score"] == pytest.approx(376.15)
    assert csm["spectrum_file"] == "XLpeplib_Beveridge_QEx-HFX_DSS_R1.raw"
    assert csm["scan_nr"] == 23454
    assert csm["charge"] == 4
    assert csm["retention_time"] == pytest.approx(6278.5199999999995)
    assert csm["ion_mobility"] == pytest.approx(0.0)

    assert last_csm["data_type"] == "crosslink-spectrum-match"
    assert last_csm["completeness"] == "full"
    assert last_csm["alpha_peptide"] == "MEDESKLHKFKDFK"
    assert mts(last_csm["alpha_modifications"]) == "(11:[DSS|138.06808])"
    assert last_csm["alpha_peptide_crosslink_position"] == 11
    assert last_csm["alpha_proteins"] == ["sp"]
    assert last_csm["alpha_proteins_crosslink_positions"] == [109]
    assert last_csm["alpha_proteins_peptide_positions"] == [99]
    assert last_csm["alpha_score"] == pytest.approx(15.89)
    assert last_csm["alpha_decoy"]
    assert last_csm["beta_peptide"] == "SSFEKNPIDFLEAK"
    assert mts(last_csm["beta_modifications"]) == "(5:[DSS|138.06808])"
    assert last_csm["beta_peptide_crosslink_position"] == 5
    assert last_csm["beta_proteins"] == ["Cas9"]
    assert last_csm["beta_proteins_crosslink_positions"] == [1180]
    assert last_csm["beta_proteins_peptide_positions"] == [1176]
    assert last_csm["beta_score"] == pytest.approx(151.04)
    assert not last_csm["beta_decoy"]
    assert last_csm["crosslink_type"] == "inter"
    assert last_csm["score"] == pytest.approx(15.89)
    assert last_csm["spectrum_file"] == "XLpeplib_Beveridge_QEx-HFX_DSS_R1.raw"
    assert last_csm["scan_nr"] == 27087
    assert last_csm["charge"] == 3
    assert last_csm["retention_time"] == pytest.approx(7421.657999999999)
    assert last_csm["ion_mobility"] == pytest.approx(0.0)


def test6():
    from pyXLMS import parser as p
    from pyXLMS.transform import modifications_to_str as mts

    parser_result = p.read_msannika(CSMS_XLSX)
    assert parser_result["data_type"] == "parser_result"
    assert parser_result["completeness"] == "partial"
    assert parser_result["search_engine"] == "MS Annika"
    assert parser_result["crosslink-spectrum-matches"] is not None
    assert parser_result["crosslinks"] is None

    csms = parser_result["crosslink-spectrum-matches"]
    assert len(csms) == 826

    first_csm = csms[0]
    csm = csms[822]
    last_csm = csms[825]

    assert first_csm["data_type"] == "crosslink-spectrum-match"
    assert first_csm["completeness"] == "full"
    assert first_csm["alpha_peptide"] == "GQKNSR"
    assert mts(first_csm["alpha_modifications"]) == "(3:[DSS|138.06808])"
    assert first_csm["alpha_peptide_crosslink_position"] == 3
    assert first_csm["alpha_proteins"] == ["Cas9"]
    assert first_csm["alpha_proteins_crosslink_positions"] == [779]
    assert first_csm["alpha_proteins_peptide_positions"] == [777]
    assert first_csm["alpha_score"] == pytest.approx(119.83)
    assert not first_csm["alpha_decoy"]
    assert first_csm["beta_peptide"] == "GQKNSR"
    assert mts(first_csm["beta_modifications"]) == "(3:[DSS|138.06808])"
    assert first_csm["beta_peptide_crosslink_position"] == 3
    assert first_csm["beta_proteins"] == ["Cas9"]
    assert first_csm["beta_proteins_crosslink_positions"] == [779]
    assert first_csm["beta_proteins_peptide_positions"] == [777]
    assert first_csm["beta_score"] == pytest.approx(119.83)
    assert not first_csm["beta_decoy"]
    assert first_csm["crosslink_type"] == "intra"
    assert first_csm["score"] == pytest.approx(119.83)
    assert first_csm["spectrum_file"] == "XLpeplib_Beveridge_QEx-HFX_DSS_R1.raw"
    assert first_csm["scan_nr"] == 2257
    assert first_csm["charge"] == 3
    assert first_csm["retention_time"] == pytest.approx(733.188)
    assert first_csm["ion_mobility"] == pytest.approx(0.0)

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
    assert csm["alpha_score"] == pytest.approx(376.15)
    assert not csm["alpha_decoy"]
    assert csm["beta_peptide"] == "TILDFLKSDGFANR"
    assert mts(csm["beta_modifications"]) == "(7:[DSS|138.06808])"
    assert csm["beta_peptide_crosslink_position"] == 7
    assert csm["beta_proteins"] == ["Cas9"]
    assert csm["beta_proteins_crosslink_positions"] == [688]
    assert csm["beta_proteins_peptide_positions"] == [682]
    assert csm["beta_score"] == pytest.approx(393.87)
    assert not csm["beta_decoy"]
    assert csm["crosslink_type"] == "intra"
    assert csm["score"] == pytest.approx(376.15)
    assert csm["spectrum_file"] == "XLpeplib_Beveridge_QEx-HFX_DSS_R1.raw"
    assert csm["scan_nr"] == 23454
    assert csm["charge"] == 4
    assert csm["retention_time"] == pytest.approx(6278.5199999999995)
    assert csm["ion_mobility"] == pytest.approx(0.0)

    assert last_csm["data_type"] == "crosslink-spectrum-match"
    assert last_csm["completeness"] == "full"
    assert last_csm["alpha_peptide"] == "MEDESKLHKFKDFK"
    assert mts(last_csm["alpha_modifications"]) == "(11:[DSS|138.06808])"
    assert last_csm["alpha_peptide_crosslink_position"] == 11
    assert last_csm["alpha_proteins"] == ["sp"]
    assert last_csm["alpha_proteins_crosslink_positions"] == [109]
    assert last_csm["alpha_proteins_peptide_positions"] == [99]
    assert last_csm["alpha_score"] == pytest.approx(15.89)
    assert last_csm["alpha_decoy"]
    assert last_csm["beta_peptide"] == "SSFEKNPIDFLEAK"
    assert mts(last_csm["beta_modifications"]) == "(5:[DSS|138.06808])"
    assert last_csm["beta_peptide_crosslink_position"] == 5
    assert last_csm["beta_proteins"] == ["Cas9"]
    assert last_csm["beta_proteins_crosslink_positions"] == [1180]
    assert last_csm["beta_proteins_peptide_positions"] == [1176]
    assert last_csm["beta_score"] == pytest.approx(151.04)
    assert not last_csm["beta_decoy"]
    assert last_csm["crosslink_type"] == "inter"
    assert last_csm["score"] == pytest.approx(15.89)
    assert last_csm["spectrum_file"] == "XLpeplib_Beveridge_QEx-HFX_DSS_R1.raw"
    assert last_csm["scan_nr"] == 27087
    assert last_csm["charge"] == 3
    assert last_csm["retention_time"] == pytest.approx(7421.657999999999)
    assert last_csm["ion_mobility"] == pytest.approx(0.0)


def test7():
    from pyXLMS import parser as p
    from pyXLMS.transform import modifications_to_str as mts

    parser_result = p.read_msannika([CSMS_TSV, CSMS_XLSX])
    assert parser_result["data_type"] == "parser_result"
    assert parser_result["completeness"] == "partial"
    assert parser_result["search_engine"] == "MS Annika"
    assert parser_result["crosslink-spectrum-matches"] is not None
    assert parser_result["crosslinks"] is None

    csms = parser_result["crosslink-spectrum-matches"]
    assert len(csms) == 826 * 2

    first_csm = csms[0]
    csm = csms[822]
    last_csm = csms[826 * 2 - 1]

    assert first_csm["data_type"] == "crosslink-spectrum-match"
    assert first_csm["completeness"] == "full"
    assert first_csm["alpha_peptide"] == "GQKNSR"
    assert mts(first_csm["alpha_modifications"]) == "(3:[DSS|138.06808])"
    assert first_csm["alpha_peptide_crosslink_position"] == 3
    assert first_csm["alpha_proteins"] == ["Cas9"]
    assert first_csm["alpha_proteins_crosslink_positions"] == [779]
    assert first_csm["alpha_proteins_peptide_positions"] == [777]
    assert first_csm["alpha_score"] == pytest.approx(119.83)
    assert not first_csm["alpha_decoy"]
    assert first_csm["beta_peptide"] == "GQKNSR"
    assert mts(first_csm["beta_modifications"]) == "(3:[DSS|138.06808])"
    assert first_csm["beta_peptide_crosslink_position"] == 3
    assert first_csm["beta_proteins"] == ["Cas9"]
    assert first_csm["beta_proteins_crosslink_positions"] == [779]
    assert first_csm["beta_proteins_peptide_positions"] == [777]
    assert first_csm["beta_score"] == pytest.approx(119.83)
    assert not first_csm["beta_decoy"]
    assert first_csm["crosslink_type"] == "intra"
    assert first_csm["score"] == pytest.approx(119.83)
    assert first_csm["spectrum_file"] == "XLpeplib_Beveridge_QEx-HFX_DSS_R1.raw"
    assert first_csm["scan_nr"] == 2257
    assert first_csm["charge"] == 3
    assert first_csm["retention_time"] == pytest.approx(733.188)
    assert first_csm["ion_mobility"] == pytest.approx(0.0)

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
    assert csm["alpha_score"] == pytest.approx(376.15)
    assert not csm["alpha_decoy"]
    assert csm["beta_peptide"] == "TILDFLKSDGFANR"
    assert mts(csm["beta_modifications"]) == "(7:[DSS|138.06808])"
    assert csm["beta_peptide_crosslink_position"] == 7
    assert csm["beta_proteins"] == ["Cas9"]
    assert csm["beta_proteins_crosslink_positions"] == [688]
    assert csm["beta_proteins_peptide_positions"] == [682]
    assert csm["beta_score"] == pytest.approx(393.87)
    assert not csm["beta_decoy"]
    assert csm["crosslink_type"] == "intra"
    assert csm["score"] == pytest.approx(376.15)
    assert csm["spectrum_file"] == "XLpeplib_Beveridge_QEx-HFX_DSS_R1.raw"
    assert csm["scan_nr"] == 23454
    assert csm["charge"] == 4
    assert csm["retention_time"] == pytest.approx(6278.5199999999995)
    assert csm["ion_mobility"] == pytest.approx(0.0)

    assert last_csm["data_type"] == "crosslink-spectrum-match"
    assert last_csm["completeness"] == "full"
    assert last_csm["alpha_peptide"] == "MEDESKLHKFKDFK"
    assert mts(last_csm["alpha_modifications"]) == "(11:[DSS|138.06808])"
    assert last_csm["alpha_peptide_crosslink_position"] == 11
    assert last_csm["alpha_proteins"] == ["sp"]
    assert last_csm["alpha_proteins_crosslink_positions"] == [109]
    assert last_csm["alpha_proteins_peptide_positions"] == [99]
    assert last_csm["alpha_score"] == pytest.approx(15.89)
    assert last_csm["alpha_decoy"]
    assert last_csm["beta_peptide"] == "SSFEKNPIDFLEAK"
    assert mts(last_csm["beta_modifications"]) == "(5:[DSS|138.06808])"
    assert last_csm["beta_peptide_crosslink_position"] == 5
    assert last_csm["beta_proteins"] == ["Cas9"]
    assert last_csm["beta_proteins_crosslink_positions"] == [1180]
    assert last_csm["beta_proteins_peptide_positions"] == [1176]
    assert last_csm["beta_score"] == pytest.approx(151.04)
    assert not last_csm["beta_decoy"]
    assert last_csm["crosslink_type"] == "inter"
    assert last_csm["score"] == pytest.approx(15.89)
    assert last_csm["spectrum_file"] == "XLpeplib_Beveridge_QEx-HFX_DSS_R1.raw"
    assert last_csm["scan_nr"] == 27087
    assert last_csm["charge"] == 3
    assert last_csm["retention_time"] == pytest.approx(7421.657999999999)
    assert last_csm["ion_mobility"] == pytest.approx(0.0)


def test8():
    from pyXLMS import parser as p
    from pyXLMS.transform import modifications_to_str as mts

    parser_result = p.read_msannika([CSMS_TSV, CSMS_TSV])
    assert parser_result["data_type"] == "parser_result"
    assert parser_result["completeness"] == "partial"
    assert parser_result["search_engine"] == "MS Annika"
    assert parser_result["crosslink-spectrum-matches"] is not None
    assert parser_result["crosslinks"] is None

    csms = parser_result["crosslink-spectrum-matches"]
    assert len(csms) == 826 * 2

    first_csm = csms[0]
    csm = csms[822]
    last_csm = csms[826 * 2 - 1]

    assert first_csm["data_type"] == "crosslink-spectrum-match"
    assert first_csm["completeness"] == "full"
    assert first_csm["alpha_peptide"] == "GQKNSR"
    assert mts(first_csm["alpha_modifications"]) == "(3:[DSS|138.06808])"
    assert first_csm["alpha_peptide_crosslink_position"] == 3
    assert first_csm["alpha_proteins"] == ["Cas9"]
    assert first_csm["alpha_proteins_crosslink_positions"] == [779]
    assert first_csm["alpha_proteins_peptide_positions"] == [777]
    assert first_csm["alpha_score"] == pytest.approx(119.83)
    assert not first_csm["alpha_decoy"]
    assert first_csm["beta_peptide"] == "GQKNSR"
    assert mts(first_csm["beta_modifications"]) == "(3:[DSS|138.06808])"
    assert first_csm["beta_peptide_crosslink_position"] == 3
    assert first_csm["beta_proteins"] == ["Cas9"]
    assert first_csm["beta_proteins_crosslink_positions"] == [779]
    assert first_csm["beta_proteins_peptide_positions"] == [777]
    assert first_csm["beta_score"] == pytest.approx(119.83)
    assert not first_csm["beta_decoy"]
    assert first_csm["crosslink_type"] == "intra"
    assert first_csm["score"] == pytest.approx(119.83)
    assert first_csm["spectrum_file"] == "XLpeplib_Beveridge_QEx-HFX_DSS_R1.raw"
    assert first_csm["scan_nr"] == 2257
    assert first_csm["charge"] == 3
    assert first_csm["retention_time"] == pytest.approx(733.188)
    assert first_csm["ion_mobility"] == pytest.approx(0.0)

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
    assert csm["alpha_score"] == pytest.approx(376.15)
    assert not csm["alpha_decoy"]
    assert csm["beta_peptide"] == "TILDFLKSDGFANR"
    assert mts(csm["beta_modifications"]) == "(7:[DSS|138.06808])"
    assert csm["beta_peptide_crosslink_position"] == 7
    assert csm["beta_proteins"] == ["Cas9"]
    assert csm["beta_proteins_crosslink_positions"] == [688]
    assert csm["beta_proteins_peptide_positions"] == [682]
    assert csm["beta_score"] == pytest.approx(393.87)
    assert not csm["beta_decoy"]
    assert csm["crosslink_type"] == "intra"
    assert csm["score"] == pytest.approx(376.15)
    assert csm["spectrum_file"] == "XLpeplib_Beveridge_QEx-HFX_DSS_R1.raw"
    assert csm["scan_nr"] == 23454
    assert csm["charge"] == 4
    assert csm["retention_time"] == pytest.approx(6278.5199999999995)
    assert csm["ion_mobility"] == pytest.approx(0.0)

    assert last_csm["data_type"] == "crosslink-spectrum-match"
    assert last_csm["completeness"] == "full"
    assert last_csm["alpha_peptide"] == "MEDESKLHKFKDFK"
    assert mts(last_csm["alpha_modifications"]) == "(11:[DSS|138.06808])"
    assert last_csm["alpha_peptide_crosslink_position"] == 11
    assert last_csm["alpha_proteins"] == ["sp"]
    assert last_csm["alpha_proteins_crosslink_positions"] == [109]
    assert last_csm["alpha_proteins_peptide_positions"] == [99]
    assert last_csm["alpha_score"] == pytest.approx(15.89)
    assert last_csm["alpha_decoy"]
    assert last_csm["beta_peptide"] == "SSFEKNPIDFLEAK"
    assert mts(last_csm["beta_modifications"]) == "(5:[DSS|138.06808])"
    assert last_csm["beta_peptide_crosslink_position"] == 5
    assert last_csm["beta_proteins"] == ["Cas9"]
    assert last_csm["beta_proteins_crosslink_positions"] == [1180]
    assert last_csm["beta_proteins_peptide_positions"] == [1176]
    assert last_csm["beta_score"] == pytest.approx(151.04)
    assert not last_csm["beta_decoy"]
    assert last_csm["crosslink_type"] == "inter"
    assert last_csm["score"] == pytest.approx(15.89)
    assert last_csm["spectrum_file"] == "XLpeplib_Beveridge_QEx-HFX_DSS_R1.raw"
    assert last_csm["scan_nr"] == 27087
    assert last_csm["charge"] == 3
    assert last_csm["retention_time"] == pytest.approx(7421.657999999999)
    assert last_csm["ion_mobility"] == pytest.approx(0.0)


def test9():
    from pyXLMS import parser as p
    from pyXLMS.transform import modifications_to_str as mts

    parser_result = p.read_msannika([XL_XLSX, CSMS_XLSX])
    assert parser_result["data_type"] == "parser_result"
    assert parser_result["completeness"] == "full"
    assert parser_result["search_engine"] == "MS Annika"
    assert parser_result["crosslink-spectrum-matches"] is not None
    assert parser_result["crosslinks"] is not None

    crosslinks = parser_result["crosslinks"]
    assert len(crosslinks) == 300

    first_crosslink = crosslinks[0]
    last_crosslink = crosslinks[299]

    csms = parser_result["crosslink-spectrum-matches"]
    assert len(csms) == 826

    first_csm = csms[0]
    csm = csms[822]
    last_csm = csms[825]

    assert first_crosslink["data_type"] == "crosslink"
    assert first_crosslink["completeness"] == "full"
    assert first_crosslink["alpha_peptide"] == "GQKNSR"
    assert first_crosslink["alpha_peptide_crosslink_position"] == 3
    assert first_crosslink["alpha_proteins"] == ["Cas9"]
    assert first_crosslink["alpha_proteins_crosslink_positions"] == [779]
    assert not first_crosslink["alpha_decoy"]
    assert first_crosslink["beta_peptide"] == "GQKNSR"
    assert first_crosslink["beta_peptide_crosslink_position"] == 3
    assert first_crosslink["beta_proteins"] == ["Cas9"]
    assert first_crosslink["beta_proteins_crosslink_positions"] == [779]
    assert not first_crosslink["beta_decoy"]
    assert first_crosslink["crosslink_type"] == "intra"
    assert first_crosslink["score"] == pytest.approx(119.83)

    assert last_crosslink["data_type"] == "crosslink"
    assert last_crosslink["completeness"] == "full"
    assert last_crosslink["alpha_peptide"] == "MEDESKLHKFKDFK"
    assert last_crosslink["alpha_peptide_crosslink_position"] == 11
    assert last_crosslink["alpha_proteins"] == ["sp"]
    assert last_crosslink["alpha_proteins_crosslink_positions"] == [109]
    assert last_crosslink["alpha_decoy"]
    assert last_crosslink["beta_peptide"] == "SSFEKNPIDFLEAK"
    assert last_crosslink["beta_peptide_crosslink_position"] == 5
    assert last_crosslink["beta_proteins"] == ["Cas9"]
    assert last_crosslink["beta_proteins_crosslink_positions"] == [1180]
    assert last_crosslink["beta_decoy"]
    assert last_crosslink["crosslink_type"] == "inter"
    assert last_crosslink["score"] == pytest.approx(15.89)

    assert first_csm["data_type"] == "crosslink-spectrum-match"
    assert first_csm["completeness"] == "full"
    assert first_csm["alpha_peptide"] == "GQKNSR"
    assert mts(first_csm["alpha_modifications"]) == "(3:[DSS|138.06808])"
    assert first_csm["alpha_peptide_crosslink_position"] == 3
    assert first_csm["alpha_proteins"] == ["Cas9"]
    assert first_csm["alpha_proteins_crosslink_positions"] == [779]
    assert first_csm["alpha_proteins_peptide_positions"] == [777]
    assert first_csm["alpha_score"] == pytest.approx(119.83)
    assert not first_csm["alpha_decoy"]
    assert first_csm["beta_peptide"] == "GQKNSR"
    assert mts(first_csm["beta_modifications"]) == "(3:[DSS|138.06808])"
    assert first_csm["beta_peptide_crosslink_position"] == 3
    assert first_csm["beta_proteins"] == ["Cas9"]
    assert first_csm["beta_proteins_crosslink_positions"] == [779]
    assert first_csm["beta_proteins_peptide_positions"] == [777]
    assert first_csm["beta_score"] == pytest.approx(119.83)
    assert not first_csm["beta_decoy"]
    assert first_csm["crosslink_type"] == "intra"
    assert first_csm["score"] == pytest.approx(119.83)
    assert first_csm["spectrum_file"] == "XLpeplib_Beveridge_QEx-HFX_DSS_R1.raw"
    assert first_csm["scan_nr"] == 2257
    assert first_csm["charge"] == 3
    assert first_csm["retention_time"] == pytest.approx(733.188)
    assert first_csm["ion_mobility"] == pytest.approx(0.0)

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
    assert csm["alpha_score"] == pytest.approx(376.15)
    assert not csm["alpha_decoy"]
    assert csm["beta_peptide"] == "TILDFLKSDGFANR"
    assert mts(csm["beta_modifications"]) == "(7:[DSS|138.06808])"
    assert csm["beta_peptide_crosslink_position"] == 7
    assert csm["beta_proteins"] == ["Cas9"]
    assert csm["beta_proteins_crosslink_positions"] == [688]
    assert csm["beta_proteins_peptide_positions"] == [682]
    assert csm["beta_score"] == pytest.approx(393.87)
    assert not csm["beta_decoy"]
    assert csm["crosslink_type"] == "intra"
    assert csm["score"] == pytest.approx(376.15)
    assert csm["spectrum_file"] == "XLpeplib_Beveridge_QEx-HFX_DSS_R1.raw"
    assert csm["scan_nr"] == 23454
    assert csm["charge"] == 4
    assert csm["retention_time"] == pytest.approx(6278.5199999999995)
    assert csm["ion_mobility"] == pytest.approx(0.0)

    assert last_csm["data_type"] == "crosslink-spectrum-match"
    assert last_csm["completeness"] == "full"
    assert last_csm["alpha_peptide"] == "MEDESKLHKFKDFK"
    assert mts(last_csm["alpha_modifications"]) == "(11:[DSS|138.06808])"
    assert last_csm["alpha_peptide_crosslink_position"] == 11
    assert last_csm["alpha_proteins"] == ["sp"]
    assert last_csm["alpha_proteins_crosslink_positions"] == [109]
    assert last_csm["alpha_proteins_peptide_positions"] == [99]
    assert last_csm["alpha_score"] == pytest.approx(15.89)
    assert last_csm["alpha_decoy"]
    assert last_csm["beta_peptide"] == "SSFEKNPIDFLEAK"
    assert mts(last_csm["beta_modifications"]) == "(5:[DSS|138.06808])"
    assert last_csm["beta_peptide_crosslink_position"] == 5
    assert last_csm["beta_proteins"] == ["Cas9"]
    assert last_csm["beta_proteins_crosslink_positions"] == [1180]
    assert last_csm["beta_proteins_peptide_positions"] == [1176]
    assert last_csm["beta_score"] == pytest.approx(151.04)
    assert not last_csm["beta_decoy"]
    assert last_csm["crosslink_type"] == "inter"
    assert last_csm["score"] == pytest.approx(15.89)
    assert last_csm["spectrum_file"] == "XLpeplib_Beveridge_QEx-HFX_DSS_R1.raw"
    assert last_csm["scan_nr"] == 27087
    assert last_csm["charge"] == 3
    assert last_csm["retention_time"] == pytest.approx(7421.657999999999)
    assert last_csm["ion_mobility"] == pytest.approx(0.0)


def test10():
    from pyXLMS import parser as p
    from pyXLMS.transform import modifications_to_str as mts

    parser_result = p.read_msannika([XL_XLSX, CSMS_TSV])
    assert parser_result["data_type"] == "parser_result"
    assert parser_result["completeness"] == "full"
    assert parser_result["search_engine"] == "MS Annika"
    assert parser_result["crosslink-spectrum-matches"] is not None
    assert parser_result["crosslinks"] is not None

    crosslinks = parser_result["crosslinks"]
    assert len(crosslinks) == 300

    first_crosslink = crosslinks[0]
    last_crosslink = crosslinks[299]

    csms = parser_result["crosslink-spectrum-matches"]
    assert len(csms) == 826

    first_csm = csms[0]
    csm = csms[822]
    last_csm = csms[825]

    assert first_crosslink["data_type"] == "crosslink"
    assert first_crosslink["completeness"] == "full"
    assert first_crosslink["alpha_peptide"] == "GQKNSR"
    assert first_crosslink["alpha_peptide_crosslink_position"] == 3
    assert first_crosslink["alpha_proteins"] == ["Cas9"]
    assert first_crosslink["alpha_proteins_crosslink_positions"] == [779]
    assert not first_crosslink["alpha_decoy"]
    assert first_crosslink["beta_peptide"] == "GQKNSR"
    assert first_crosslink["beta_peptide_crosslink_position"] == 3
    assert first_crosslink["beta_proteins"] == ["Cas9"]
    assert first_crosslink["beta_proteins_crosslink_positions"] == [779]
    assert not first_crosslink["beta_decoy"]
    assert first_crosslink["crosslink_type"] == "intra"
    assert first_crosslink["score"] == pytest.approx(119.83)

    assert last_crosslink["data_type"] == "crosslink"
    assert last_crosslink["completeness"] == "full"
    assert last_crosslink["alpha_peptide"] == "MEDESKLHKFKDFK"
    assert last_crosslink["alpha_peptide_crosslink_position"] == 11
    assert last_crosslink["alpha_proteins"] == ["sp"]
    assert last_crosslink["alpha_proteins_crosslink_positions"] == [109]
    assert last_crosslink["alpha_decoy"]
    assert last_crosslink["beta_peptide"] == "SSFEKNPIDFLEAK"
    assert last_crosslink["beta_peptide_crosslink_position"] == 5
    assert last_crosslink["beta_proteins"] == ["Cas9"]
    assert last_crosslink["beta_proteins_crosslink_positions"] == [1180]
    assert last_crosslink["beta_decoy"]
    assert last_crosslink["crosslink_type"] == "inter"
    assert last_crosslink["score"] == pytest.approx(15.89)

    assert first_csm["data_type"] == "crosslink-spectrum-match"
    assert first_csm["completeness"] == "full"
    assert first_csm["alpha_peptide"] == "GQKNSR"
    assert mts(first_csm["alpha_modifications"]) == "(3:[DSS|138.06808])"
    assert first_csm["alpha_peptide_crosslink_position"] == 3
    assert first_csm["alpha_proteins"] == ["Cas9"]
    assert first_csm["alpha_proteins_crosslink_positions"] == [779]
    assert first_csm["alpha_proteins_peptide_positions"] == [777]
    assert first_csm["alpha_score"] == pytest.approx(119.83)
    assert not first_csm["alpha_decoy"]
    assert first_csm["beta_peptide"] == "GQKNSR"
    assert mts(first_csm["beta_modifications"]) == "(3:[DSS|138.06808])"
    assert first_csm["beta_peptide_crosslink_position"] == 3
    assert first_csm["beta_proteins"] == ["Cas9"]
    assert first_csm["beta_proteins_crosslink_positions"] == [779]
    assert first_csm["beta_proteins_peptide_positions"] == [777]
    assert first_csm["beta_score"] == pytest.approx(119.83)
    assert not first_csm["beta_decoy"]
    assert first_csm["crosslink_type"] == "intra"
    assert first_csm["score"] == pytest.approx(119.83)
    assert first_csm["spectrum_file"] == "XLpeplib_Beveridge_QEx-HFX_DSS_R1.raw"
    assert first_csm["scan_nr"] == 2257
    assert first_csm["charge"] == 3
    assert first_csm["retention_time"] == pytest.approx(733.188)
    assert first_csm["ion_mobility"] == pytest.approx(0.0)

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
    assert csm["alpha_score"] == pytest.approx(376.15)
    assert not csm["alpha_decoy"]
    assert csm["beta_peptide"] == "TILDFLKSDGFANR"
    assert mts(csm["beta_modifications"]) == "(7:[DSS|138.06808])"
    assert csm["beta_peptide_crosslink_position"] == 7
    assert csm["beta_proteins"] == ["Cas9"]
    assert csm["beta_proteins_crosslink_positions"] == [688]
    assert csm["beta_proteins_peptide_positions"] == [682]
    assert csm["beta_score"] == pytest.approx(393.87)
    assert not csm["beta_decoy"]
    assert csm["crosslink_type"] == "intra"
    assert csm["score"] == pytest.approx(376.15)
    assert csm["spectrum_file"] == "XLpeplib_Beveridge_QEx-HFX_DSS_R1.raw"
    assert csm["scan_nr"] == 23454
    assert csm["charge"] == 4
    assert csm["retention_time"] == pytest.approx(6278.5199999999995)
    assert csm["ion_mobility"] == pytest.approx(0.0)

    assert last_csm["data_type"] == "crosslink-spectrum-match"
    assert last_csm["completeness"] == "full"
    assert last_csm["alpha_peptide"] == "MEDESKLHKFKDFK"
    assert mts(last_csm["alpha_modifications"]) == "(11:[DSS|138.06808])"
    assert last_csm["alpha_peptide_crosslink_position"] == 11
    assert last_csm["alpha_proteins"] == ["sp"]
    assert last_csm["alpha_proteins_crosslink_positions"] == [109]
    assert last_csm["alpha_proteins_peptide_positions"] == [99]
    assert last_csm["alpha_score"] == pytest.approx(15.89)
    assert last_csm["alpha_decoy"]
    assert last_csm["beta_peptide"] == "SSFEKNPIDFLEAK"
    assert mts(last_csm["beta_modifications"]) == "(5:[DSS|138.06808])"
    assert last_csm["beta_peptide_crosslink_position"] == 5
    assert last_csm["beta_proteins"] == ["Cas9"]
    assert last_csm["beta_proteins_crosslink_positions"] == [1180]
    assert last_csm["beta_proteins_peptide_positions"] == [1176]
    assert last_csm["beta_score"] == pytest.approx(151.04)
    assert not last_csm["beta_decoy"]
    assert last_csm["crosslink_type"] == "inter"
    assert last_csm["score"] == pytest.approx(15.89)
    assert last_csm["spectrum_file"] == "XLpeplib_Beveridge_QEx-HFX_DSS_R1.raw"
    assert last_csm["scan_nr"] == 27087
    assert last_csm["charge"] == 3
    assert last_csm["retention_time"] == pytest.approx(7421.657999999999)
    assert last_csm["ion_mobility"] == pytest.approx(0.0)


def test11():
    from pyXLMS import parser as p
    from pyXLMS.transform import modifications_to_str as mts

    parser_result = p.read_msannika([XL_TSV, CSMS_XLSX])
    assert parser_result["data_type"] == "parser_result"
    assert parser_result["completeness"] == "full"
    assert parser_result["search_engine"] == "MS Annika"
    assert parser_result["crosslink-spectrum-matches"] is not None
    assert parser_result["crosslinks"] is not None

    crosslinks = parser_result["crosslinks"]
    assert len(crosslinks) == 300

    first_crosslink = crosslinks[0]
    last_crosslink = crosslinks[299]

    csms = parser_result["crosslink-spectrum-matches"]
    assert len(csms) == 826

    first_csm = csms[0]
    csm = csms[822]
    last_csm = csms[825]

    assert first_crosslink["data_type"] == "crosslink"
    assert first_crosslink["completeness"] == "full"
    assert first_crosslink["alpha_peptide"] == "GQKNSR"
    assert first_crosslink["alpha_peptide_crosslink_position"] == 3
    assert first_crosslink["alpha_proteins"] == ["Cas9"]
    assert first_crosslink["alpha_proteins_crosslink_positions"] == [779]
    assert not first_crosslink["alpha_decoy"]
    assert first_crosslink["beta_peptide"] == "GQKNSR"
    assert first_crosslink["beta_peptide_crosslink_position"] == 3
    assert first_crosslink["beta_proteins"] == ["Cas9"]
    assert first_crosslink["beta_proteins_crosslink_positions"] == [779]
    assert not first_crosslink["beta_decoy"]
    assert first_crosslink["crosslink_type"] == "intra"
    assert first_crosslink["score"] == pytest.approx(119.83)

    assert last_crosslink["data_type"] == "crosslink"
    assert last_crosslink["completeness"] == "full"
    assert last_crosslink["alpha_peptide"] == "MEDESKLHKFKDFK"
    assert last_crosslink["alpha_peptide_crosslink_position"] == 11
    assert last_crosslink["alpha_proteins"] == ["sp"]
    assert last_crosslink["alpha_proteins_crosslink_positions"] == [109]
    assert last_crosslink["alpha_decoy"]
    assert last_crosslink["beta_peptide"] == "SSFEKNPIDFLEAK"
    assert last_crosslink["beta_peptide_crosslink_position"] == 5
    assert last_crosslink["beta_proteins"] == ["Cas9"]
    assert last_crosslink["beta_proteins_crosslink_positions"] == [1180]
    assert last_crosslink["beta_decoy"]
    assert last_crosslink["crosslink_type"] == "inter"
    assert last_crosslink["score"] == pytest.approx(15.89)

    assert first_csm["data_type"] == "crosslink-spectrum-match"
    assert first_csm["completeness"] == "full"
    assert first_csm["alpha_peptide"] == "GQKNSR"
    assert mts(first_csm["alpha_modifications"]) == "(3:[DSS|138.06808])"
    assert first_csm["alpha_peptide_crosslink_position"] == 3
    assert first_csm["alpha_proteins"] == ["Cas9"]
    assert first_csm["alpha_proteins_crosslink_positions"] == [779]
    assert first_csm["alpha_proteins_peptide_positions"] == [777]
    assert first_csm["alpha_score"] == pytest.approx(119.83)
    assert not first_csm["alpha_decoy"]
    assert first_csm["beta_peptide"] == "GQKNSR"
    assert mts(first_csm["beta_modifications"]) == "(3:[DSS|138.06808])"
    assert first_csm["beta_peptide_crosslink_position"] == 3
    assert first_csm["beta_proteins"] == ["Cas9"]
    assert first_csm["beta_proteins_crosslink_positions"] == [779]
    assert first_csm["beta_proteins_peptide_positions"] == [777]
    assert first_csm["beta_score"] == pytest.approx(119.83)
    assert not first_csm["beta_decoy"]
    assert first_csm["crosslink_type"] == "intra"
    assert first_csm["score"] == pytest.approx(119.83)
    assert first_csm["spectrum_file"] == "XLpeplib_Beveridge_QEx-HFX_DSS_R1.raw"
    assert first_csm["scan_nr"] == 2257
    assert first_csm["charge"] == 3
    assert first_csm["retention_time"] == pytest.approx(733.188)
    assert first_csm["ion_mobility"] == pytest.approx(0.0)

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
    assert csm["alpha_score"] == pytest.approx(376.15)
    assert not csm["alpha_decoy"]
    assert csm["beta_peptide"] == "TILDFLKSDGFANR"
    assert mts(csm["beta_modifications"]) == "(7:[DSS|138.06808])"
    assert csm["beta_peptide_crosslink_position"] == 7
    assert csm["beta_proteins"] == ["Cas9"]
    assert csm["beta_proteins_crosslink_positions"] == [688]
    assert csm["beta_proteins_peptide_positions"] == [682]
    assert csm["beta_score"] == pytest.approx(393.87)
    assert not csm["beta_decoy"]
    assert csm["crosslink_type"] == "intra"
    assert csm["score"] == pytest.approx(376.15)
    assert csm["spectrum_file"] == "XLpeplib_Beveridge_QEx-HFX_DSS_R1.raw"
    assert csm["scan_nr"] == 23454
    assert csm["charge"] == 4
    assert csm["retention_time"] == pytest.approx(6278.5199999999995)
    assert csm["ion_mobility"] == pytest.approx(0.0)

    assert last_csm["data_type"] == "crosslink-spectrum-match"
    assert last_csm["completeness"] == "full"
    assert last_csm["alpha_peptide"] == "MEDESKLHKFKDFK"
    assert mts(last_csm["alpha_modifications"]) == "(11:[DSS|138.06808])"
    assert last_csm["alpha_peptide_crosslink_position"] == 11
    assert last_csm["alpha_proteins"] == ["sp"]
    assert last_csm["alpha_proteins_crosslink_positions"] == [109]
    assert last_csm["alpha_proteins_peptide_positions"] == [99]
    assert last_csm["alpha_score"] == pytest.approx(15.89)
    assert last_csm["alpha_decoy"]
    assert last_csm["beta_peptide"] == "SSFEKNPIDFLEAK"
    assert mts(last_csm["beta_modifications"]) == "(5:[DSS|138.06808])"
    assert last_csm["beta_peptide_crosslink_position"] == 5
    assert last_csm["beta_proteins"] == ["Cas9"]
    assert last_csm["beta_proteins_crosslink_positions"] == [1180]
    assert last_csm["beta_proteins_peptide_positions"] == [1176]
    assert last_csm["beta_score"] == pytest.approx(151.04)
    assert not last_csm["beta_decoy"]
    assert last_csm["crosslink_type"] == "inter"
    assert last_csm["score"] == pytest.approx(15.89)
    assert last_csm["spectrum_file"] == "XLpeplib_Beveridge_QEx-HFX_DSS_R1.raw"
    assert last_csm["scan_nr"] == 27087
    assert last_csm["charge"] == 3
    assert last_csm["retention_time"] == pytest.approx(7421.657999999999)
    assert last_csm["ion_mobility"] == pytest.approx(0.0)


def test12():
    from pyXLMS import parser as p
    from pyXLMS.transform import modifications_to_str as mts

    parser_result = p.read_msannika([XL_TSV, CSMS_TSV])
    assert parser_result["data_type"] == "parser_result"
    assert parser_result["completeness"] == "full"
    assert parser_result["search_engine"] == "MS Annika"
    assert parser_result["crosslink-spectrum-matches"] is not None
    assert parser_result["crosslinks"] is not None

    crosslinks = parser_result["crosslinks"]
    assert len(crosslinks) == 300

    first_crosslink = crosslinks[0]
    last_crosslink = crosslinks[299]

    csms = parser_result["crosslink-spectrum-matches"]
    assert len(csms) == 826

    first_csm = csms[0]
    csm = csms[822]
    last_csm = csms[825]

    assert first_crosslink["data_type"] == "crosslink"
    assert first_crosslink["completeness"] == "full"
    assert first_crosslink["alpha_peptide"] == "GQKNSR"
    assert first_crosslink["alpha_peptide_crosslink_position"] == 3
    assert first_crosslink["alpha_proteins"] == ["Cas9"]
    assert first_crosslink["alpha_proteins_crosslink_positions"] == [779]
    assert not first_crosslink["alpha_decoy"]
    assert first_crosslink["beta_peptide"] == "GQKNSR"
    assert first_crosslink["beta_peptide_crosslink_position"] == 3
    assert first_crosslink["beta_proteins"] == ["Cas9"]
    assert first_crosslink["beta_proteins_crosslink_positions"] == [779]
    assert not first_crosslink["beta_decoy"]
    assert first_crosslink["crosslink_type"] == "intra"
    assert first_crosslink["score"] == pytest.approx(119.83)

    assert last_crosslink["data_type"] == "crosslink"
    assert last_crosslink["completeness"] == "full"
    assert last_crosslink["alpha_peptide"] == "MEDESKLHKFKDFK"
    assert last_crosslink["alpha_peptide_crosslink_position"] == 11
    assert last_crosslink["alpha_proteins"] == ["sp"]
    assert last_crosslink["alpha_proteins_crosslink_positions"] == [109]
    assert last_crosslink["alpha_decoy"]
    assert last_crosslink["beta_peptide"] == "SSFEKNPIDFLEAK"
    assert last_crosslink["beta_peptide_crosslink_position"] == 5
    assert last_crosslink["beta_proteins"] == ["Cas9"]
    assert last_crosslink["beta_proteins_crosslink_positions"] == [1180]
    assert last_crosslink["beta_decoy"]
    assert last_crosslink["crosslink_type"] == "inter"
    assert last_crosslink["score"] == pytest.approx(15.89)

    assert first_csm["data_type"] == "crosslink-spectrum-match"
    assert first_csm["completeness"] == "full"
    assert first_csm["alpha_peptide"] == "GQKNSR"
    assert mts(first_csm["alpha_modifications"]) == "(3:[DSS|138.06808])"
    assert first_csm["alpha_peptide_crosslink_position"] == 3
    assert first_csm["alpha_proteins"] == ["Cas9"]
    assert first_csm["alpha_proteins_crosslink_positions"] == [779]
    assert first_csm["alpha_proteins_peptide_positions"] == [777]
    assert first_csm["alpha_score"] == pytest.approx(119.83)
    assert not first_csm["alpha_decoy"]
    assert first_csm["beta_peptide"] == "GQKNSR"
    assert mts(first_csm["beta_modifications"]) == "(3:[DSS|138.06808])"
    assert first_csm["beta_peptide_crosslink_position"] == 3
    assert first_csm["beta_proteins"] == ["Cas9"]
    assert first_csm["beta_proteins_crosslink_positions"] == [779]
    assert first_csm["beta_proteins_peptide_positions"] == [777]
    assert first_csm["beta_score"] == pytest.approx(119.83)
    assert not first_csm["beta_decoy"]
    assert first_csm["crosslink_type"] == "intra"
    assert first_csm["score"] == pytest.approx(119.83)
    assert first_csm["spectrum_file"] == "XLpeplib_Beveridge_QEx-HFX_DSS_R1.raw"
    assert first_csm["scan_nr"] == 2257
    assert first_csm["charge"] == 3
    assert first_csm["retention_time"] == pytest.approx(733.188)
    assert first_csm["ion_mobility"] == pytest.approx(0.0)

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
    assert csm["alpha_score"] == pytest.approx(376.15)
    assert not csm["alpha_decoy"]
    assert csm["beta_peptide"] == "TILDFLKSDGFANR"
    assert mts(csm["beta_modifications"]) == "(7:[DSS|138.06808])"
    assert csm["beta_peptide_crosslink_position"] == 7
    assert csm["beta_proteins"] == ["Cas9"]
    assert csm["beta_proteins_crosslink_positions"] == [688]
    assert csm["beta_proteins_peptide_positions"] == [682]
    assert csm["beta_score"] == pytest.approx(393.87)
    assert not csm["beta_decoy"]
    assert csm["crosslink_type"] == "intra"
    assert csm["score"] == pytest.approx(376.15)
    assert csm["spectrum_file"] == "XLpeplib_Beveridge_QEx-HFX_DSS_R1.raw"
    assert csm["scan_nr"] == 23454
    assert csm["charge"] == 4
    assert csm["retention_time"] == pytest.approx(6278.5199999999995)
    assert csm["ion_mobility"] == pytest.approx(0.0)

    assert last_csm["data_type"] == "crosslink-spectrum-match"
    assert last_csm["completeness"] == "full"
    assert last_csm["alpha_peptide"] == "MEDESKLHKFKDFK"
    assert mts(last_csm["alpha_modifications"]) == "(11:[DSS|138.06808])"
    assert last_csm["alpha_peptide_crosslink_position"] == 11
    assert last_csm["alpha_proteins"] == ["sp"]
    assert last_csm["alpha_proteins_crosslink_positions"] == [109]
    assert last_csm["alpha_proteins_peptide_positions"] == [99]
    assert last_csm["alpha_score"] == pytest.approx(15.89)
    assert last_csm["alpha_decoy"]
    assert last_csm["beta_peptide"] == "SSFEKNPIDFLEAK"
    assert mts(last_csm["beta_modifications"]) == "(5:[DSS|138.06808])"
    assert last_csm["beta_peptide_crosslink_position"] == 5
    assert last_csm["beta_proteins"] == ["Cas9"]
    assert last_csm["beta_proteins_crosslink_positions"] == [1180]
    assert last_csm["beta_proteins_peptide_positions"] == [1176]
    assert last_csm["beta_score"] == pytest.approx(151.04)
    assert not last_csm["beta_decoy"]
    assert last_csm["crosslink_type"] == "inter"
    assert last_csm["score"] == pytest.approx(15.89)
    assert last_csm["spectrum_file"] == "XLpeplib_Beveridge_QEx-HFX_DSS_R1.raw"
    assert last_csm["scan_nr"] == 27087
    assert last_csm["charge"] == 3
    assert last_csm["retention_time"] == pytest.approx(7421.657999999999)
    assert last_csm["ion_mobility"] == pytest.approx(0.0)


def test13():
    from pyXLMS import parser as p
    from pyXLMS.transform import modifications_to_str as mts

    parser_result = p.read_msannika([XL_XLSX, CSMS_XLSX, CSMS_TSV])
    assert parser_result["data_type"] == "parser_result"
    assert parser_result["completeness"] == "full"
    assert parser_result["search_engine"] == "MS Annika"
    assert parser_result["crosslink-spectrum-matches"] is not None
    assert parser_result["crosslinks"] is not None

    crosslinks = parser_result["crosslinks"]
    assert len(crosslinks) == 300

    first_crosslink = crosslinks[0]
    last_crosslink = crosslinks[299]

    csms = parser_result["crosslink-spectrum-matches"]
    assert len(csms) == 826 * 2

    first_csm = csms[0]
    csm = csms[822]
    last_csm = csms[826 * 2 - 1]

    assert first_crosslink["data_type"] == "crosslink"
    assert first_crosslink["completeness"] == "full"
    assert first_crosslink["alpha_peptide"] == "GQKNSR"
    assert first_crosslink["alpha_peptide_crosslink_position"] == 3
    assert first_crosslink["alpha_proteins"] == ["Cas9"]
    assert first_crosslink["alpha_proteins_crosslink_positions"] == [779]
    assert not first_crosslink["alpha_decoy"]
    assert first_crosslink["beta_peptide"] == "GQKNSR"
    assert first_crosslink["beta_peptide_crosslink_position"] == 3
    assert first_crosslink["beta_proteins"] == ["Cas9"]
    assert first_crosslink["beta_proteins_crosslink_positions"] == [779]
    assert not first_crosslink["beta_decoy"]
    assert first_crosslink["crosslink_type"] == "intra"
    assert first_crosslink["score"] == pytest.approx(119.83)

    assert last_crosslink["data_type"] == "crosslink"
    assert last_crosslink["completeness"] == "full"
    assert last_crosslink["alpha_peptide"] == "MEDESKLHKFKDFK"
    assert last_crosslink["alpha_peptide_crosslink_position"] == 11
    assert last_crosslink["alpha_proteins"] == ["sp"]
    assert last_crosslink["alpha_proteins_crosslink_positions"] == [109]
    assert last_crosslink["alpha_decoy"]
    assert last_crosslink["beta_peptide"] == "SSFEKNPIDFLEAK"
    assert last_crosslink["beta_peptide_crosslink_position"] == 5
    assert last_crosslink["beta_proteins"] == ["Cas9"]
    assert last_crosslink["beta_proteins_crosslink_positions"] == [1180]
    assert last_crosslink["beta_decoy"]
    assert last_crosslink["crosslink_type"] == "inter"
    assert last_crosslink["score"] == pytest.approx(15.89)

    assert first_csm["data_type"] == "crosslink-spectrum-match"
    assert first_csm["completeness"] == "full"
    assert first_csm["alpha_peptide"] == "GQKNSR"
    assert mts(first_csm["alpha_modifications"]) == "(3:[DSS|138.06808])"
    assert first_csm["alpha_peptide_crosslink_position"] == 3
    assert first_csm["alpha_proteins"] == ["Cas9"]
    assert first_csm["alpha_proteins_crosslink_positions"] == [779]
    assert first_csm["alpha_proteins_peptide_positions"] == [777]
    assert first_csm["alpha_score"] == pytest.approx(119.83)
    assert not first_csm["alpha_decoy"]
    assert first_csm["beta_peptide"] == "GQKNSR"
    assert mts(first_csm["beta_modifications"]) == "(3:[DSS|138.06808])"
    assert first_csm["beta_peptide_crosslink_position"] == 3
    assert first_csm["beta_proteins"] == ["Cas9"]
    assert first_csm["beta_proteins_crosslink_positions"] == [779]
    assert first_csm["beta_proteins_peptide_positions"] == [777]
    assert first_csm["beta_score"] == pytest.approx(119.83)
    assert not first_csm["beta_decoy"]
    assert first_csm["crosslink_type"] == "intra"
    assert first_csm["score"] == pytest.approx(119.83)
    assert first_csm["spectrum_file"] == "XLpeplib_Beveridge_QEx-HFX_DSS_R1.raw"
    assert first_csm["scan_nr"] == 2257
    assert first_csm["charge"] == 3
    assert first_csm["retention_time"] == pytest.approx(733.188)
    assert first_csm["ion_mobility"] == pytest.approx(0.0)

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
    assert csm["alpha_score"] == pytest.approx(376.15)
    assert not csm["alpha_decoy"]
    assert csm["beta_peptide"] == "TILDFLKSDGFANR"
    assert mts(csm["beta_modifications"]) == "(7:[DSS|138.06808])"
    assert csm["beta_peptide_crosslink_position"] == 7
    assert csm["beta_proteins"] == ["Cas9"]
    assert csm["beta_proteins_crosslink_positions"] == [688]
    assert csm["beta_proteins_peptide_positions"] == [682]
    assert csm["beta_score"] == pytest.approx(393.87)
    assert not csm["beta_decoy"]
    assert csm["crosslink_type"] == "intra"
    assert csm["score"] == pytest.approx(376.15)
    assert csm["spectrum_file"] == "XLpeplib_Beveridge_QEx-HFX_DSS_R1.raw"
    assert csm["scan_nr"] == 23454
    assert csm["charge"] == 4
    assert csm["retention_time"] == pytest.approx(6278.5199999999995)
    assert csm["ion_mobility"] == pytest.approx(0.0)

    assert last_csm["data_type"] == "crosslink-spectrum-match"
    assert last_csm["completeness"] == "full"
    assert last_csm["alpha_peptide"] == "MEDESKLHKFKDFK"
    assert mts(last_csm["alpha_modifications"]) == "(11:[DSS|138.06808])"
    assert last_csm["alpha_peptide_crosslink_position"] == 11
    assert last_csm["alpha_proteins"] == ["sp"]
    assert last_csm["alpha_proteins_crosslink_positions"] == [109]
    assert last_csm["alpha_proteins_peptide_positions"] == [99]
    assert last_csm["alpha_score"] == pytest.approx(15.89)
    assert last_csm["alpha_decoy"]
    assert last_csm["beta_peptide"] == "SSFEKNPIDFLEAK"
    assert mts(last_csm["beta_modifications"]) == "(5:[DSS|138.06808])"
    assert last_csm["beta_peptide_crosslink_position"] == 5
    assert last_csm["beta_proteins"] == ["Cas9"]
    assert last_csm["beta_proteins_crosslink_positions"] == [1180]
    assert last_csm["beta_proteins_peptide_positions"] == [1176]
    assert last_csm["beta_score"] == pytest.approx(151.04)
    assert not last_csm["beta_decoy"]
    assert last_csm["crosslink_type"] == "inter"
    assert last_csm["score"] == pytest.approx(15.89)
    assert last_csm["spectrum_file"] == "XLpeplib_Beveridge_QEx-HFX_DSS_R1.raw"
    assert last_csm["scan_nr"] == 27087
    assert last_csm["charge"] == 3
    assert last_csm["retention_time"] == pytest.approx(7421.657999999999)
    assert last_csm["ion_mobility"] == pytest.approx(0.0)


def test14():
    from pyXLMS import parser as p

    file_extension = ".tar.gz"
    unsupp_file = XL_XLSX + file_extension
    with pytest.raises(
        ValueError,
        match="Detected file extension .gz is not supported! Input file has to be a valid file with extension '.csv', '.tsv', '.xlsx', or '.pdResult'!",
    ):
        _r = p.read_msannika(unsupp_file)


def test15():
    from pyXLMS import parser as p

    format = "sqlite"
    with pytest.raises(
        ValueError,
        match=f"Provided input format {format} is not supported! Input format has to be of type 'csv', 'tsv', 'xlsx', or 'pdresult'!",
    ):
        _r = p.read_msannika(XL_XLSX, format=format)


def test16():
    from pyXLMS import parser as p

    modifications = {"DSS": 138.06808, "Carbamidomethyl": 57.021464}
    err_str = (
        "Unable to find modification Oxidation in the set of provided modifications. "
        + "Please pass the full set of expected modifications to the parser."
    )
    with pytest.raises(KeyError, match=err_str):
        _r = p.read_msannika(CSMS_XLSX, modifications=modifications)


def test17():
    # unfiltered CSMs: decoy peptides without protein accession and protein position
    from pyXLMS.parser import read_msannika

    pr = read_msannika(CSMS_TSV_PXD043531)
    csms = pr["crosslink-spectrum-matches"]
    assert len(csms) == 5
    assert [csm["scan_nr"] for csm in csms] == [1241, 1243, 1244, 1245, 1246]

    # target-target
    csm = csms[0]
    assert csm["alpha_proteins"] == ["Q14980"]
    assert csm["alpha_proteins_crosslink_positions"] == [708]
    assert csm["beta_proteins_crosslink_positions"] == [715]
    assert not csm["alpha_decoy"] and not csm["beta_decoy"]

    # target-decoy
    csm = csms[1]
    assert csm["alpha_peptide"] == "GGEAQKKENMEWAELIR"
    assert csm["alpha_peptide_crosslink_position"] == 7
    assert csm["alpha_proteins"] is None
    assert csm["alpha_proteins_crosslink_positions"] is None
    assert csm["alpha_proteins_peptide_positions"] is None
    assert csm["alpha_decoy"]
    assert csm["beta_peptide"] == "KSFLMKR"
    assert csm["beta_proteins"] == ["Q9ULH0"]
    assert csm["beta_proteins_crosslink_positions"] == [1448]
    assert csm["beta_proteins_peptide_positions"] == [1448]
    assert not csm["beta_decoy"]

    # decoy-decoy
    csm = csms[2]
    assert csm["alpha_proteins"] is None
    assert csm["beta_proteins"] is None
    assert csm["alpha_decoy"] and csm["beta_decoy"]


def test18():
    # unfiltered crosslinks: decoy peptides without protein have position -1
    from pyXLMS.parser import read_msannika

    pr = read_msannika(XL_TSV_PXD043531)
    xls = pr["crosslinks"]
    assert len(xls) == 4

    xl = xls[1]
    assert xl["alpha_peptide"] == "GGEAQKKENMEWAELIR"
    assert xl["alpha_peptide_crosslink_position"] == 7
    assert xl["alpha_proteins"] is None
    assert xl["alpha_proteins_crosslink_positions"] is None
    assert xl["beta_peptide"] == "KSFLMKR"
    assert xl["beta_peptide_crosslink_position"] == 1
    assert xl["beta_proteins"] == ["Q9ULH0"]
    assert xl["beta_proteins_crosslink_positions"] == [1448]
    # decoy crosslink of a decoy peptide (no protein) and a target peptide, as labelled
    # in the matching crosslink-spectrum-match (scan 1243: KSFLMKR T, GGEAQKKENMEWAELIR D)
    assert xl["alpha_decoy"]
    assert not xl["beta_decoy"]

    xl = xls[2]
    assert xl["alpha_peptide_crosslink_position"] == 13
    assert xl["beta_peptide_crosslink_position"] == 5
    assert xl["alpha_proteins"] is None
    assert xl["beta_proteins"] is None
