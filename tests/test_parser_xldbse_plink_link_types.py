#!/usr/bin/env python3

# pyXLMS - TESTS
# 2026 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

MONO = "data/plink2/link_type_examples/A549_DSS_example.filtered_mono-linked_spectra.csv"
LOOP = "data/plink2/link_type_examples/A549_DSS_example.filtered_loop-linked_spectra.csv"
CROSS = "data/plink2/link_type_examples/A549_DSS_example.filtered_cross-linked_spectra.csv"


def test1():
    # Mono-linked (dead-end) spectra file -> MonoLink objects, no crosslinks/csms.
    from pyXLMS.parser import read_plink

    pr = read_plink(MONO, verbose=0)
    assert pr["data_type"] == "parser_result"
    assert pr.crosslinks is None
    assert pr.crosslink_spectrum_matches is None
    assert pr.mono_links is not None

    mono_links = pr.mono_links
    assert len(mono_links) == 8
    assert all(ml.data_type == "mono-link" for ml in mono_links)
    assert all(isinstance(ml.peptide_crosslink_position, int) for ml in mono_links)

    ml = next(ml for ml in mono_links if ml.peptide == "APKPDGPGGGPGGSHMGGNYGDDR")
    assert ml.peptide_crosslink_position == 3
    assert ml.proteins == ["sp|P35637|FUS_HUMAN"]
    assert ml.proteins_crosslink_positions == [451]
    assert ml.decoy is False


def test2():
    # Loop-linked spectra file -> intra crosslinks (both sites on the same peptide),
    # so alpha and beta refer to the same peptide/protein.
    from pyXLMS.parser import read_plink

    pr = read_plink(LOOP, verbose=0)
    assert pr.mono_links is None
    csms = pr.crosslink_spectrum_matches
    assert csms is not None and len(csms) == 8
    # every loop-link is intra: same peptide on both sides
    assert all(csm.alpha_peptide == csm.beta_peptide for csm in csms)
    assert all(
        csm.additional_information.get("link_type") == "loop-link" for csm in csms
    )

    csm = next(csm for csm in csms if csm.alpha_proteins == ["sp|P26038|MOES_HUMAN"])
    assert csm.alpha_peptide == "DQKKTQEQLALEMAELTAR"
    # the two crosslinked sites of the loop are 3 and 4 (order-independent)
    assert {csm.alpha_peptide_crosslink_position, csm.beta_peptide_crosslink_position} == {3, 4}
    assert {
        csm.alpha_proteins_crosslink_positions[0],
        csm.beta_proteins_crosslink_positions[0],
    } == {411, 412}


def test3():
    # Cross-linked spectra file still reads exactly as before: inter/intra crosslinks
    # between two distinct peptides, and no mono-links.
    from pyXLMS.parser import read_plink

    pr = read_plink(CROSS, verbose=0)
    assert pr.mono_links is None
    csms = pr.crosslink_spectrum_matches
    assert csms is not None and len(csms) == 8
    # at least one genuine two-peptide crosslink
    assert any(csm.alpha_peptide != csm.beta_peptide for csm in csms)


def test4():
    # MonoLink model basics.
    from pyXLMS.data import MonoLink, create_mono_link, create_mono_link_min

    ml = create_mono_link_min("PEPTIDE", 1)
    assert isinstance(ml, MonoLink)
    assert ml.data_type == "mono-link"
    assert ml.completeness == "partial"
    assert ml["peptide"] == "PEPTIDE"

    full = create_mono_link(
        "PEPTIDE", 1, proteins=["PROT"], xl_position_proteins=[10], decoy=False, score=1.0
    )
    assert full.completeness == "full"
    assert "mono-link" == full.model_dump()["data_type"]
