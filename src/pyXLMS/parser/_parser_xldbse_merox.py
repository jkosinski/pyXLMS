#!/usr/bin/env python3

# 2025 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

from __future__ import annotations

import io
import re
import warnings
import zipfile
import pandas as pd
from tqdm import tqdm

from ..data._parser_result import ParserResult
from ..data._util import check_input
from ..data._csm import create_csm
from ..data._parser_result import create_parser_result
from ..constants import AMINO_ACIDS
from ..constants import MODIFICATIONS
from ..constants import MEROX_MODIFICATION_MAPPING
from ._util import __serialize_pandas_series
from ._util import __parse_int, __parse_float

from typing import Optional
from typing import BinaryIO
from typing import Dict
from typing import Any
from typing import Tuple
from typing import List


MEROX_COLNAMES = [
    "Score",
    "m/z",
    "Charge",
    "M+H+",
    "Calculated Mass",
    "Deviation in ppm",
    "Peptide 1",
    "Protein 1",
    "From",
    "To",
    "Peptide2",
    "Protein 2",
    "From.1",
    "To.1",
    "Scan number",
    "is Selected in Table",
    "Candidate identifier",
    "Folder Number",
    "Retention time in sec",
    "miscellaneous",
    "best linkage position peptide 1",
    "best linkage position peptide 2",
    "All linkage positions",
    "Spectrum UUID",
    "local False discovery rate",
    "shortest distance in pdb",
    "Light/Heavy(1/2)",
    "pepScore1",
    "pepScore2",
    "xLinkScore",
    "resultId",
    "MS1intensity",
    "finalScoreComponent",
]


def __read_merox_file(
    file: str | BinaryIO, sep: str = ";", decimal: str = ".", **kwargs
) -> pd.DataFrame:
    r"""Helper function to read MeroX files into pandas DataFrames.

    Reads MeroX files into pandas DataFrames independent of input format. Accepts
    both ``.csv`` and ``.zhrm`` files.

    Parameters
    ----------
    file : str, or file stream
        The name/path of the MeroX result file or a file-like object/stream.
    sep : str, default = ";"
        Seperator used in the ``.csv`` or ``.zhrm`` file.
    decimal : str, default = "."
        Character to recognize as decimal point.
    **kwargs
        Any additional parameters will be passed to ``pandas.read*``.

    Returns
    -------
    pd.DataFrame
        The MeroX result file as a pandas DataFrame.

    Notes
    -----
    This function should not be called directly, it is called from ``read_merox()``.
    """
    # safety
    if not isinstance(file, str):
        file.seek(0)
    # this shifts file pointer
    if zipfile.is_zipfile(file):
        with zipfile.ZipFile(file, "r") as f:
            return pd.read_csv(  # ty: ignore[unsound-return-statement]
                io.BytesIO(f.read("Result.csv")),
                header=None,
                names=MEROX_COLNAMES,
                sep=sep,
                decimal=decimal,
                low_memory=False,
                **kwargs,
            )
    if not isinstance(file, str):
        file.seek(0)
    return pd.read_csv(file, sep=sep, decimal=decimal, low_memory=False, **kwargs)  # ty: ignore[unsound-return-statement]


def __get_merox_sequence(
    sequence: str,
    parse_modifications: bool = True,
    modifications: Dict[str, Dict[str, Any]] = MEROX_MODIFICATION_MAPPING,
) -> str:
    r"""Helper funtion to read the MeroX peptide sequence.

    Reads the MeroX peptide sequence and removes prefix and suffix square brackets,
    then replaces modification symbols with their respective amino acids. If
    ``parse_modifications = True`` non-resolved symbols will raise an error,
    otherwise the symbol will be used in the sequence - even if it is not a valid
    amino acid sequence anymore.

    Parameters
    ----------
    sequence : str
        The MeroX sequence e.g. from column "Peptide 1" or "Peptide2".
    parse_modifications : bool, default = True
        Whether or not post-translational-modifications should be parsed for crosslink-spectrum-matches.
        Requires correct specification of the 'modifications' parameter.
    modifications: dict of str, dict of str, any, default = ``constants.MEROX_MODIFICATION_MAPPING``
        Mapping of modification symbols to their amino acids and modifications. Please refer to
        ``constants.MEROX_MODIFICATION_MAPPING`` for examples.

    Returns
    -------
    str
        The parsed sequence.

    Raises
    ------
    KeyError
        If one of the symbols in the sequence could not be resolved/an unknown modification is encountered.

    Notes
    -----
    This function should not be called directly, it is called from ``read_merox()``.
    """
    parsed_seq = ""
    for amino_acid in sequence.lstrip("[").rstrip("]").strip():
        if amino_acid in modifications:
            parsed_seq += modifications[amino_acid]["Amino Acid"]
        elif amino_acid in AMINO_ACIDS:
            parsed_seq += amino_acid
        else:
            if parse_modifications:
                raise KeyError(
                    f"Key {amino_acid} not found in parameter 'modifications'. Are you missing a modification?"
                )
            else:
                parsed_seq += amino_acid
    return parsed_seq  # ty: ignore[unsound-return-statement]


def __get_merox_modifications(
    sequence: str,
    crosslink_position: int,
    crosslinker: str,
    crosslinker_mass: float,
    modifications: Dict[str, Dict[str, Any]] = MEROX_MODIFICATION_MAPPING,
) -> Dict[int, Tuple[str, float]]:
    r"""Helper function to read modifications from a MeroX sequence.

    Parses post-translational-modifications from the given MeroX sequence.

    Parameters
    ----------
    sequence : str
        The MeroX sequence e.g. from column "Peptide 1" or "Peptide2".
    crosslink_position : int
        Position of the crosslinker in the sequence (1-based).
    crosslinker : str
        Name of the used cross-linking reagent, for example "DSSO".
    crosslinker_mass : float
        Monoisotopic delta mass of the crosslink modification.
    modifications: dict of str, dict of str, any, default = ``constants.MEROX_MODIFICATION_MAPPING``
        Mapping of modification symbols to their amino acids and modifications. Please refer to
        ``constants.MEROX_MODIFICATION_MAPPING`` for examples.

    Returns
    -------
    dict of int, tuple of str, float
        The pyXLMS specfic modification representation of the parsed modifications.

    Raises
    ------
    RuntimeError
        If multiple modifications on the same residue are parsed.
    KeyError
        If an unknown modification is encountered.

    Notes
    -----
    This function should not be called directly, it is called from ``read_merox()``.
    """
    parsed_modifications = {crosslink_position: (crosslinker, crosslinker_mass)}
    for i, amino_acid in enumerate(sequence.lstrip("[").rstrip("]").strip()):
        if amino_acid in modifications:
            if i + 1 not in parsed_modifications:
                parsed_modifications[i + 1] = modifications[amino_acid]["Modification"]
            else:
                raise RuntimeError(f"Modification at position {i + 1} already exists!")
        elif amino_acid in AMINO_ACIDS:
            pass
        else:
            raise KeyError(
                f"Key {amino_acid} not found in parameter 'modifications'. Are you missing a modification?"
            )
    return parsed_modifications  # ty: ignore[unsound-return-statement]


def __get_merox_position(position_str: str) -> int:
    r"""Helper function to extract the peptide crosslink position from MeroX.

    Parameters
    ----------
    position_str : str
        The position string from MeroX e.g. from column "best linkage position peptide 1".

    Returns
    -------
    int
        The parsed peptide crosslink position.

    Raises
    ------
    RuntimeError
        If the position could not be parsed.

    Notes
    -----
    This function should not be called directly, it is called from ``read_merox()``.
    """
    position = None
    try:
        position = __parse_int(position_str[1:])
    except Exception as e:
        raise RuntimeError(f"Could not parse position from {position_str}!") from e
    return position


def __get_merox_protein(proteins: str) -> List[str]:
    r"""Helper function to extract the protein accession from MeroX.

    Parses the leading protein accession from the MeroX protein string. Additional
    proteins in parentheses are not parsed as they do not have corresponding protein
    crosslink positions. A peptide shared by several proteins is written as their
    accessions joined by ``" and "`` followed by the peptide's range, e.g.
    ``">sp|P0DMV9|HS71B_HUMAN and >sp|P0DMV8|HS71A_HUMAN@510-517"``; each protein is
    returned once.

    Parameters
    ----------
    proteins : str
        The protein string from MeroX e.g. from column "Protein 1".

    Returns
    -------
    list of str
        A list containing the parsed protein accession(s).

    Notes
    -----
    This function should not be called directly, it is called from ``read_merox()``.
    """
    accessions = list()
    for protein in __split_merox_shared_proteins(proteins):
        accession = protein.split("(>")[0].strip().lstrip(">")
        if accession not in accessions:
            accessions.append(accession)
    return accessions


def __split_merox_shared_proteins(proteins: str) -> List[str]:
    r"""Helper function to split a MeroX protein string into the proteins of a shared peptide.

    MeroX joins the proteins of a peptide shared by several proteins with ``" and "`` and
    appends the peptide's range, e.g. ``">A and >B@510-517"``. The string is split at each
    ``" and "`` followed by a protein (``">"`` or ``"DEC_>"``) outside parentheses, and the
    range is removed from the last protein.

    Parameters
    ----------
    proteins : str
        The protein string from MeroX e.g. from column "Protein 1".

    Returns
    -------
    list of str
        The protein strings, one for a peptide that is not shared.

    Notes
    -----
    This function should not be called directly, it is called from ``read_merox()``.
    """
    parts = list()
    depth = 0
    start = 0
    for match in re.finditer(r"[()]|\s+and\s+(?=(?:DEC_)?>)", proteins):
        separator = match.group(0)
        if separator == "(":
            depth += 1
        elif separator == ")":
            depth = max(depth - 1, 0)
        elif depth == 0:
            parts.append(proteins[start : match.start()])
            start = match.end()
    parts.append(proteins[start:])
    if len(parts) > 1:
        parts[-1] = re.sub(r"@\d*-?\d*\s*$", "", parts[-1])
    return parts


def __get_merox_scan_number(scan_nr_and_file: str) -> int:
    r"""Helper function to parse the scan number from MeroX.

    Parses the scan number from the MeroX scan number string. Assumes that scan
    number and spectrum file are delimted by the tilde (wave) sign.

    Parameters
    ----------
    scan_nr_and_file : str
        The scan number string from MeroX e.g. from the column "Scan number".

    Returns
    -------
    int
        The parsed scan number.

    Raises
    ------
    RuntimeError
        If the scan number could not be parsed.

    Notes
    -----
    This function should not be called directly, it is called from ``read_merox()``.
    """
    scan_nr = None
    try:
        scan_nr = __parse_int(scan_nr_and_file.split("~")[0])
    except Exception as e:
        raise RuntimeError(f"Could not scan number from {scan_nr_and_file}!") from e
    return scan_nr


def __get_merox_spectrum_file(scan_nr_and_file: str) -> str:
    r"""Helper function to parse the spectrum file from MeroX.

    Parses the spectrum file from the MeroX scan number string. Assumes that scan
    number and spectrum file are delimted by the tilde (wave) sign.

    Parameters
    ----------
    scan_nr_and_file : str
        The scan number string from MeroX e.g. from the column "Scan number".

    Returns
    -------
    str
        The parsed spectrum file.

    Raises
    ------
    RuntimeError
        If the spectrum file could not be parsed.

    Notes
    -----
    This function should not be called directly, it is called from ``read_merox()``.
    """
    spectrum_file = None
    try:
        spectrum_file = scan_nr_and_file.split("~")[1]
    except Exception as e:
        raise RuntimeError(
            f"Could not spectrum file name from {scan_nr_and_file}!"
        ) from e
    return spectrum_file


def parse_site_from_merox(
    peptide: str,
    protein: str,
    peptide_start: Any,
    position: str,
    parse_modifications: bool = True,
    modifications: Dict[str, Dict[str, Any]] = MEROX_MODIFICATION_MAPPING,
) -> Dict[str, Any]:
    r"""Parse a MeroX peptide, protein, start, and linkage position into a linked site.

    Parses the fields of one peptide of a MeroX result row, e.g. ``Peptide 1``,
    ``Protein 1``, ``From``, and ``best linkage position peptide 1``, into its linked
    site, in a representation that is independent of the pyXLMS data model so it can
    be reused to build any downstream object.

    Parameters
    ----------
    peptide : str
        The MeroX sequence e.g. from column "Peptide 1" or "Peptide2".
    protein : str
        The MeroX protein string e.g. from column "Protein 1" or "Protein 2".
    peptide_start : any
        The position of the peptide in the protein e.g. from column "From".
    position : str
        The MeroX linkage position e.g. from column "best linkage position peptide 1".
    parse_modifications : bool, default = True
        Whether or not unknown symbols in the sequence should raise an error.
    modifications: dict of str, dict of str, any, default = ``constants.MEROX_MODIFICATION_MAPPING``
        Mapping of modification symbols to their amino acids and modifications. Please refer to
        ``constants.MEROX_MODIFICATION_MAPPING`` for examples.

    Returns
    -------
    dict of str, any
        The linked site, a dictionary with keys ``peptide`` (str), ``peptide_position`` (int,
        1-based, in-peptide), ``proteins`` (list of str, accessions) and ``protein_positions``
        (list of int, 1-based, absolute; parallel to ``proteins``). For a peptide shared by
        several proteins, or found several times in one (accessions joined by ``" and "``),
        the protein positions are ``None``: MeroX reports one position of the peptide only.

    Raises
    ------
    KeyError
        If one of the symbols in the sequence could not be resolved/an unknown modification is
        encountered (only if ``parse_modifications = True``).
    RuntimeError
        If the linkage position could not be parsed.
    ValueError
        If the peptide start could not be parsed.

    Examples
    --------
    >>> from pyXLMS.parser import parse_site_from_merox
    >>> site = parse_site_from_merox("[VKYVTEGmR]", ">Cas9", 531, "K2")
    >>> site["peptide"]
    'VKYVTEGMR'
    >>> site["protein_positions"]
    [532]
    """
    peptide_position = __get_merox_position(str(position))
    proteins = __get_merox_protein(str(protein))
    return {
        "peptide": __get_merox_sequence(
            str(peptide), parse_modifications, modifications
        ),
        "peptide_position": peptide_position,
        "proteins": proteins,
        "protein_positions": None
        if len(__split_merox_shared_proteins(str(protein))) > 1
        else [__parse_int(peptide_start) + peptide_position - 1],
    }


def read_merox(
    files: str | List[str] | BinaryIO,
    crosslinker: str,
    crosslinker_mass: Optional[float] = None,
    decoy_prefix: str = "REV__",
    parse_modifications: bool = True,
    modifications: Dict[str, Dict[str, Any]] = MEROX_MODIFICATION_MAPPING,
    sep: str = ";",
    decimal: str = ".",
    **kwargs,
) -> ParserResult:
    r"""Read a MeroX result file.

    Reads a MeroX crosslink-spectrum-matches result file in ``.csv`` or ``.zhrm`` format
    and returns a ``parser_result``.

    Parameters
    ----------
    files : str, list of str, or file stream
        The name/path of the MeroX result file(s) or a file-like object/stream.
    crosslinker : str
        Name of the used cross-linking reagent, for example "DSSO".
    crosslinker_mass : float, or None, default = None
        Monoisotopic delta mass of the crosslink modification. If the crosslinker is
        defined in ``constants.MODIFICATIONS`` this can be omitted.
    decoy_prefix : str, default = "REV\_\_"
        The prefix that indicates that a protein is from the decoy database.
    parse_modifications : bool, default = True
        Whether or not post-translational-modifications should be parsed for crosslink-spectrum-matches.
        Requires correct specification of the 'modifications' parameter.
    modifications: dict of str, dict of str, any, default = ``constants.MEROX_MODIFICATION_MAPPING``
        Mapping of modification symbols to their amino acids and modifications. Please refer to
        ``constants.MEROX_MODIFICATION_MAPPING`` for examples.
    sep : str, default = ";"
        Seperator used in the ``.csv`` or ``.zhrm`` file.
    decimal : str, default = "."
        Character to recognize as decimal point.
    **kwargs
        Any additional parameters will be passed to ``pandas.read*``.

    Returns
    -------
    ParserResult
        The ``parser_result`` object containing all parsed information.

    Raises
    ------
    RuntimeError
        If the file(s) could not be read or if the file(s) contain no crosslink-spectrum-matches.
    KeyError
        If the specified crosslinker could not be found/mapped.

    Notes
    -----
    Uses ``pepScore1`` as the score for the alpha peptide, ``pepScore2`` as the score of the
    beta peptide, and ``Score`` as the score of the crosslink-spectrum-match.
    A peptide shared by several proteins (accessions joined by ``" and "``, e.g.
    ``">sp|P0DMV9|HS71B_HUMAN and >sp|P0DMV8|HS71A_HUMAN@510-517"``) is read with all of its
    proteins and without protein positions, since MeroX reports the position in one of the
    proteins only; their number is reported in a warning.

    Warnings
    --------
    MeroX only reports a single protein crosslink position per peptide, for ambiguous peptides
    only the crosslink position of the first matching protein is reported. All matching proteins can be
    retrieved via ``additional_information``, however not their corresponding crosslink positions. For this
    reason it is recommended to use ``transform.reannotate_positions()`` to correctly annotate all crosslink
    positions for all peptides if that is important for downstream analysis. Additionally, please note that
    target and decoy information is derived based off the protein accession and parameter ``decoy_prefix``.
    By default, MeroX only reports target matches that are above the desired FDR.

    Examples
    --------
    >>> from pyXLMS.parser import read_merox
    >>> csms_from_csv = read_merox(
    ...     "data/merox/XLpeplib_Beveridge_QEx-HFX_DSS_R1.csv", crosslinker="DSS"
    ... )

    >>> from pyXLMS.parser import read_merox
    >>> csms_from_zhrm = read_merox(
    ...     "data/merox/XLpeplib_Beveridge_QEx-HFX_DSS_R1.zhrm", crosslinker="DSS"
    ... )
    """
    ## check input
    _ok = check_input(crosslinker, "crosslinker", str)
    _ok = (
        check_input(crosslinker_mass, "crosslinker_mass", float)
        if crosslinker_mass is not None
        else True
    )
    _ok = check_input(decoy_prefix, "decoy_prefix", str)
    _ok = check_input(parse_modifications, "parse_modifications", bool)
    _ok = check_input(modifications, "modifications", dict, dict)
    _ok = check_input(sep, "sep", str)
    _ok = check_input(decimal, "decimal", str)
    if crosslinker_mass is None:
        if crosslinker not in MODIFICATIONS:
            if parse_modifications:
                raise KeyError(
                    "Cannot infer crosslinker mass because crosslinker is unknown. "
                    "Please specify crosslinker mass manually!"
                )
            else:
                crosslinker_mass = 0.0
        else:
            crosslinker_mass = MODIFICATIONS[crosslinker]

    ## data structures
    csms = list()
    shared_peptides = 0

    ## handle input
    if not isinstance(files, list):
        inputs = [files]
    else:
        inputs = files

    for input in inputs:  # noqa: A001
        ## reading data
        data = __read_merox_file(input, sep=sep, decimal=decimal, **kwargs)
        for _i, row in tqdm(
            data.iterrows(),
            total=data.shape[0],
            desc="Reading MeroX CSMs...",
        ):
            # parse linked sites
            site_a = parse_site_from_merox(
                row["Peptide 1"],
                row["Protein 1"],
                row["From"],
                row["best linkage position peptide 1"],
                parse_modifications,
                modifications,
            )
            site_b = parse_site_from_merox(
                row["Peptide2"],
                row["Protein 2"],
                row["From.1"],
                row["best linkage position peptide 2"],
                parse_modifications,
                modifications,
            )
            shared_peptides += int(site_a["protein_positions"] is None) + int(
                site_b["protein_positions"] is None
            )
            # create crosslink
            csm = create_csm(
                peptide_a=site_a["peptide"],
                modifications_a=__get_merox_modifications(
                    str(row["Peptide 1"]),
                    __get_merox_position(str(row["best linkage position peptide 1"])),
                    crosslinker,
                    crosslinker_mass,
                    modifications,
                )
                if parse_modifications
                else None,
                xl_position_peptide_a=site_a["peptide_position"],
                proteins_a=site_a["proteins"],
                xl_position_proteins_a=site_a["protein_positions"],
                pep_position_proteins_a=[__parse_int(row["From"])]
                if site_a["protein_positions"] is not None
                else None,
                score_a=__parse_float(row["pepScore1"]),
                decoy_a=site_a["proteins"][0].startswith(decoy_prefix),
                peptide_b=site_b["peptide"],
                modifications_b=__get_merox_modifications(
                    str(row["Peptide2"]),
                    __get_merox_position(str(row["best linkage position peptide 2"])),
                    crosslinker,
                    crosslinker_mass,
                    modifications,
                )
                if parse_modifications
                else None,
                xl_position_peptide_b=site_b["peptide_position"],
                proteins_b=site_b["proteins"],
                xl_position_proteins_b=site_b["protein_positions"],
                pep_position_proteins_b=[__parse_int(row["From.1"])]
                if site_b["protein_positions"] is not None
                else None,
                score_b=__parse_float(row["pepScore2"]),
                decoy_b=site_b["proteins"][0].startswith(decoy_prefix),
                score=__parse_float(row["Score"]),
                spectrum_file=__get_merox_spectrum_file(str(row["Scan number"])),
                scan_nr=__get_merox_scan_number(str(row["Scan number"])),
                charge=__parse_int(row["Charge"]),
                rt=__parse_float(row["Retention time in sec"]),
                im_cv=None,
                additional_information={
                    "source": __serialize_pandas_series(row),
                    "xLinkScore": row["xLinkScore"],
                    "Protein 1": row["Protein 1"],
                    "Protein 2": row["Protein 2"],
                    "MS1intensity": row["MS1intensity"],
                },
            )
            csms.append(csm)
    ## report peptides shared by several proteins
    if shared_peptides > 0:
        warnings.warn(
            RuntimeWarning(
                f"Read {shared_peptides} peptide(s) shared by several proteins without protein "
                "positions, MeroX reports the position in one of the proteins only!"
            ),
            stacklevel=2,
        )
    ## check results
    if len(csms) == 0:
        raise RuntimeError(
            "No crosslink-spectrum-matches were parsed! If this is unexpected, please file a bug report!"
        )
    ## return parser result
    return create_parser_result(
        search_engine="MeroX",
        csms=csms,
        crosslinks=None,
    )
