#!/usr/bin/env python3

# pyXLMS - TESTS
# 2025 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

import pytest


F1 = "data/ms_annika/XLpeplib_Beveridge_QEx-HFX_DSS_R1.mzid"
F2 = "data/xlinkx/XLpeplib_Beveridge_Lumos_DSSO_MS3.mzid"


def test1():
    from pyXLMS.parser import parse_scan_nr_from_mzid

    assert parse_scan_nr_from_mzid("scan=5321") == 5321


@pytest.mark.slow
def test2():
    from pyXLMS import parser as p

    pr = p.read_mzid(F1, verbose=0)
    assert pr["data_type"] == "parser_result"
    assert pr["completeness"] == "partial"
    assert pr["search_engine"] == "mzIdentML"
    assert pr["crosslink-spectrum-matches"] is not None
    assert pr["crosslinks"] is None

    csms = pr["crosslink-spectrum-matches"]
    assert len(csms) == 786


@pytest.mark.slow
def test3():
    from pyXLMS import parser as p

    pr = p.read_mzid(F2, verbose=0)
    assert pr["data_type"] == "parser_result"
    assert pr["completeness"] == "partial"
    assert pr["search_engine"] == "mzIdentML"
    assert pr["crosslink-spectrum-matches"] is not None
    assert pr["crosslinks"] is None

    csms = pr["crosslink-spectrum-matches"]
    assert len(csms) == 823


def test4():
    from pyXLMS.parser import parse_scan_nr_from_mzid

    with pytest.warns(RuntimeWarning):
        assert parse_scan_nr_from_mzid("index=1") == 1


def test5():
    from pyteomics import mzid
    from pyXLMS import parser as p

    def get_nr_mzid_items(mzid_file: str) -> int:
        i = 0
        with mzid.MzIdentML(mzid_file) as reader:
            for item in reader:
                if "crosslink spectrum identification item" in str(item):
                    i += 1
                elif "cross-link spectrum identification item" in str(item):
                    i += 1
        return i

    test_files = [
        "data/scout2/Cas9_HeLa_Cyt_r1/Cas9_HeLa_Cyt_r1-v1.2.mzid",
        "data/scout2/Cas9_HeLa_Cyt_r1/Cas9_HeLa_Cyt_r1-v1.3.mzid",
        "data/scout2/Cas9_HeLa_Cyt_r2/Cas9_HeLa_Cyt_r2-v1.2.mzid",
        "data/scout2/Cas9_HeLa_Cyt_r2/Cas9_HeLa_Cyt_r2-v1.3.mzid",
        "data/scout2/Cas9_HeLa_Cyt_r3/Cas9_HeLa_Cyt_r3-v1.2.mzid",
        "data/scout2/Cas9_HeLa_Cyt_r3/Cas9_HeLa_Cyt_r3-v1.3.mzid",
    ]

    with pytest.warns(RuntimeWarning):
        for test_file in test_files:
            pr = p.read_mzid(test_file, verbose=0)
            assert pr["data_type"] == "parser_result"
            assert pr["completeness"] == "partial"
            assert pr["search_engine"] == "mzIdentML"
            assert pr["crosslink-spectrum-matches"] is not None
            assert pr["crosslinks"] is None

            csms = pr["crosslink-spectrum-matches"]
            assert len(csms) == get_nr_mzid_items(test_file)


def test6():
    from pyXLMS import parser as p

    test_files = [
        "data/scout2/Cas9_HeLa_Cyt_r1/Cas9_HeLa_Cyt_r1-v1.2.mzid",
        "data/scout2/Cas9_HeLa_Cyt_r1/Cas9_HeLa_Cyt_r1-v1.3.mzid",
    ]

    with pytest.warns(RuntimeWarning):
        for test_file in test_files:
            pr = p.read_mzid(test_file, decoy=None, verbose=0)
            assert pr["data_type"] == "parser_result"
            assert pr["completeness"] == "partial"
            assert pr["search_engine"] == "mzIdentML"
            assert pr["crosslink-spectrum-matches"] is not None
            assert pr["crosslinks"] is None

            csms = pr["crosslink-spectrum-matches"]
            for csm in csms:
                assert csm["alpha_decoy"] is None
                assert csm["beta_decoy"] is None


def test7():
    from pyXLMS import parser as p

    test_files = [
        "data/scout2/Cas9_HeLa_Cyt_r1/Cas9_HeLa_Cyt_r1-v1.2.mzid",
        "data/scout2/Cas9_HeLa_Cyt_r1/Cas9_HeLa_Cyt_r1-v1.3.mzid",
    ]

    with pytest.warns(RuntimeWarning):
        for test_file in test_files:
            pr = p.read_mzid(test_file, decoy=True, verbose=0)
            assert pr["data_type"] == "parser_result"
            assert pr["completeness"] == "partial"
            assert pr["search_engine"] == "mzIdentML"
            assert pr["crosslink-spectrum-matches"] is not None
            assert pr["crosslinks"] is None

            csms = pr["crosslink-spectrum-matches"]
            for csm in csms:
                assert csm["alpha_decoy"]
                assert csm["beta_decoy"]


def test8():
    from pyXLMS import parser as p

    test_files = [
        "data/scout2/Cas9_HeLa_Cyt_r1/Cas9_HeLa_Cyt_r1-v1.2.mzid",
        "data/scout2/Cas9_HeLa_Cyt_r1/Cas9_HeLa_Cyt_r1-v1.3.mzid",
    ]

    with pytest.warns(RuntimeWarning):
        for test_file in test_files:
            pr = p.read_mzid(test_file, decoy=False, verbose=0)
            assert pr["data_type"] == "parser_result"
            assert pr["completeness"] == "partial"
            assert pr["search_engine"] == "mzIdentML"
            assert pr["crosslink-spectrum-matches"] is not None
            assert pr["crosslinks"] is None

            csms = pr["crosslink-spectrum-matches"]
            for csm in csms:
                assert not csm["alpha_decoy"]
                assert not csm["beta_decoy"]


def test9():
    # mzIdentML passThreshold is carried through for every CSM via
    # additional_information (same frozen-record reason as test_parser_xldbse_scout@test15).
    from pyXLMS.parser import read_mzid

    # this Scout mzid has spectrum IDs without "scan=", so read_mzid emits the same
    # RuntimeWarning the other mzid tests expect.
    with pytest.warns(RuntimeWarning):
        result = read_mzid(
            "data/scout2/Cas9_HeLa_Cyt_r2/Cas9_HeLa_Cyt_r2-v1.2.mzid", verbose=0
        )
    csms = result["crosslink-spectrum-matches"]
    assert len(csms) > 0
    for csm in csms:
        assert csm.additional_information["pass_threshold"] is True


@pytest.mark.slow
def test10(tmp_path):
    # pyteomics cannot read a peptide with an empty cvParam (as written for a second
    # crosslinker definition, e.g. in PRIDE deposits PXD066061 and PXD066071), the
    # spectrum is skipped and counted in a warning
    import warnings
    from pyXLMS import parser as p

    with open(F1, "r", encoding="utf-8") as f:
        content = f.read()
    # crosslinker modification of peptide "XL_-132_1_DON" (GQKNSR)
    param = '<cvParam accession="MS:1002509" cvRef="PSI-MS" name="cross-link donor" value="17" />'
    assert content.count(param) == 1
    content = content.replace(param, "<cvParam></cvParam>" + param)
    f_unresolved = tmp_path / "unresolved.mzid"
    f_unresolved.write_text(content, encoding="utf-8")

    with warnings.catch_warnings(record=True) as wl:
        warnings.simplefilter("always")
        pr = p.read_mzid(str(f_unresolved))
    skipped = [str(w.message) for w in wl if str(w.message).startswith("Skipped")]
    assert skipped == [
        "Skipped 1 spectra with crosslink identifications whose peptides could not be "
        "read from the mzIdentML file (unresolved references)!"
    ]
    assert len(pr["crosslink-spectrum-matches"]) == 785

    pr = p.read_mzid(str(f_unresolved), verbose=0)
    assert len(pr["crosslink-spectrum-matches"]) == 785
