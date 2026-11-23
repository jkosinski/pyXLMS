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

# prefix of the decoy proteins MeroX generates, e.g. "DEC_>Cas9"
MEROX_DECOY_PREFIX = "DEC_"


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
    # rows of exported tables can end with a separator, without index_col=False pandas
    # would use the first column as index and shift all columns by one, the empty last
    # field is dropped
    kwargs.setdefault("index_col", False)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", pd.errors.ParserWarning)
        return pd.read_csv(file, sep=sep, decimal=decimal, low_memory=False, **kwargs)  # ty: ignore[unsound-return-statement]


def __get_merox_sequence(
    sequence: str,
    parse_modifications: bool = True,
    modifications: Dict[str, Dict[str, Any]] = MEROX_MODIFICATION_MAPPING,
) -> str:
    r"""Helper funtion to read the MeroX peptide sequence.

    Reads the MeroX peptide sequence and removes prefix and suffix square brackets
    (or curly brackets for protein-terminal peptides), then replaces modification symbols with their respective amino acids. If
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
    for amino_acid in sequence.lstrip("[{").rstrip("]}").strip():
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
    n_terminal_group: Optional[Tuple[str, float]] = None,
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
    n_terminal_group : tuple of str, float, or None, default = None
        Modification of a modified protein N-terminus, set at position 0 (see
        ``__get_merox_protein_n_terminal_group()``).

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
    if n_terminal_group is not None:
        if 0 in parsed_modifications:
            raise RuntimeError("Modification at position 0 already exists!")
        parsed_modifications[0] = n_terminal_group
    for i, amino_acid in enumerate(sequence.lstrip("[{").rstrip("]}").strip()):
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


def __get_merox_crosslink_position(position_str: str) -> int:
    r"""Helper function to extract the peptide crosslink position from MeroX.

    MeroX reports crosslinks to the protein N-terminus at position 0, e.g. ``"{0"``.
    These are assigned to the first residue of the peptide (position 1). All other
    positions are returned as parsed by ``__get_merox_position()``.

    Parameters
    ----------
    position_str : str
        The position string from MeroX e.g. from column "best linkage position peptide 1".

    Returns
    -------
    int
        The parsed peptide crosslink position (1-based).

    Raises
    ------
    RuntimeError
        If the position could not be parsed.

    Notes
    -----
    This function should not be called directly, it is called from ``read_merox()``.
    """
    position = __get_merox_position(position_str)
    if position == 0 and position_str.strip().startswith("{"):
        return 1
    return position


def __get_merox_peptide_position(from_position: Any, sequence: str) -> int:
    r"""Helper function to extract the peptide position in the protein from MeroX.

    MeroX reports the start of protein N-terminal peptides (sequences starting with
    ``"{"``) one position lower than the start of all other peptides, e.g. ``0`` for a
    peptide starting at the first residue of the protein. This function returns the
    1-based position for all peptides.

    Parameters
    ----------
    from_position : any
        The peptide start from MeroX e.g. from column "From".
    sequence : str
        The MeroX sequence e.g. from column "Peptide 1" or "Peptide2".

    Returns
    -------
    int
        The peptide position in the protein (1-based).

    Raises
    ------
    ValueError
        If the position could not be parsed.

    Notes
    -----
    This function should not be called directly, it is called from ``read_merox()``.
    """
    position = __parse_int(from_position)
    if str(sequence).strip().startswith("{"):
        return position + 1
    return position


def __get_merox_placed_peptide_position(from_position: Any, sequence: str) -> int:
    r"""Helper function to extract the position of a peptide placed in a protein from MeroX.

    Returns the position as ``__get_merox_peptide_position()`` and refuses a peptide that
    starts at position 0 without being protein N-terminal (``"{"``): MeroX writes such a
    peptide when its first letter stands for a modified protein N-terminus that is not
    mapped in parameter 'modifications' (see ``__get_merox_protein_n_terminal_group()``).

    Parameters
    ----------
    from_position : any
        The peptide start from MeroX e.g. from column "From".
    sequence : str
        The MeroX sequence e.g. from column "Peptide 1" or "Peptide2".

    Returns
    -------
    int
        The peptide position in the protein (1-based).

    Raises
    ------
    ValueError
        If the position could not be parsed, or the peptide starts at position 0 without
        being protein N-terminal.

    Notes
    -----
    This function should not be called directly, it is called from ``read_merox()``.
    """
    position = __get_merox_peptide_position(from_position, sequence)
    if position < 1:
        raise ValueError(
            f"MeroX peptide {str(sequence).strip()} starts at position {position} of its "
            "protein (column 'From') but is not a protein N-terminal peptide ('{'). If its "
            "first letter stands for a modified protein N-terminus (a static modification "
            "of '{' in the MeroX settings, e.g. '{;X'), map that letter to amino acid '{' "
            "in parameter 'modifications'."
        )
    return position


def __get_merox_protein_n_terminal_group(
    sequence: str,
    modifications: Dict[str, Dict[str, Any]] = MEROX_MODIFICATION_MAPPING,
) -> Tuple[str, Optional[Tuple[str, float]]]:
    r"""Helper function to read a modified protein N-terminus from a MeroX sequence.

    MeroX writes a modified protein N-terminus, a static modification of ``"{"`` (e.g.
    ``"{;X"`` for an acetylated N-terminus), as its own letter after ``"["``, e.g.
    ``"[XMDVFMKGLSK]"``, with the peptide start ("From") of a protein N-terminal peptide
    and the letter counted as peptide position 1. A letter that parameter
    'modifications' maps to amino acid ``"{"`` is read as such a group: the sequence is
    returned as a protein N-terminal peptide without the letter (``"{MDVFMKGLSK]"``),
    together with the modification of the group.

    Parameters
    ----------
    sequence : str
        The MeroX sequence e.g. from column "Peptide 1" or "Peptide2".
    modifications: dict of str, dict of str, any, default = ``constants.MEROX_MODIFICATION_MAPPING``
        Mapping of modification symbols to their amino acids and modifications. Please refer to
        ``constants.MEROX_MODIFICATION_MAPPING`` for examples.

    Returns
    -------
    tuple of str, tuple of str, float or None
        The sequence, and the modification of the protein N-terminus or ``None`` if the
        sequence has no modified protein N-terminus.

    Notes
    -----
    This function should not be called directly, it is called from ``read_merox()``.
    """
    stripped = str(sequence).strip()
    if (
        len(stripped) > 1
        and stripped[0] == "["
        and stripped[1] in modifications
        and modifications[stripped[1]]["Amino Acid"] == "{"
    ):
        return "{" + stripped[2:], modifications[stripped[1]]["Modification"]
    return sequence, None


def __get_merox_group_crosslink_position(
    position_str: str, n_terminal_group: Optional[Tuple[str, float]]
) -> int:
    r"""Helper function to extract the peptide crosslink position from MeroX.

    Returns the position as ``__get_merox_crosslink_position()``. For a peptide with a
    modified protein N-terminus (see ``__get_merox_protein_n_terminal_group()``) MeroX
    counts the letter of the N-terminus as position 1, so the position is one lower.

    Parameters
    ----------
    position_str : str
        The position string from MeroX e.g. from column "best linkage position peptide 1".
    n_terminal_group : tuple of str, float, or None
        The modification of the protein N-terminus of the peptide, or ``None``.

    Returns
    -------
    int
        The parsed peptide crosslink position (1-based).

    Raises
    ------
    RuntimeError
        If the position could not be parsed.
    ValueError
        If the crosslink is on the modified protein N-terminus itself.

    Notes
    -----
    This function should not be called directly, it is called from ``read_merox()``.
    """
    position = __get_merox_crosslink_position(position_str)
    if n_terminal_group is None:
        return position
    if position < 2:
        raise ValueError(
            f"MeroX crosslink position {position_str.strip()} is on the modified protein "
            f"N-terminus {n_terminal_group[0]}, which is not read!"
        )
    return position - 1


def __get_merox_protein(proteins: str) -> List[str]:
    r"""Helper function to extract the protein accession from MeroX.

    Parses the leading protein accession from the MeroX protein string. Additional
    proteins in parentheses are not parsed as they do not have corresponding protein
    crosslink positions. A peptide shared by several proteins is written as their
    accessions joined by ``" and "`` followed by the peptide's range, e.g.
    ``">sp|P0DMV9|HS71B_HUMAN and >sp|P0DMV8|HS71A_HUMAN@510-517"``; each protein is
    returned once. MeroX writes the decoy proteins it generates as ``"DEC_>"`` followed by
    the accession of the target protein, e.g. ``"DEC_>Cas9"``; these are returned as
    ``"DEC_Cas9"``.

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
        accession = protein.split("(>")[0].strip()
        if accession.startswith(MEROX_DECOY_PREFIX + ">"):
            accession = MEROX_DECOY_PREFIX + accession[len(MEROX_DECOY_PREFIX) + 1 :]
        accession = accession.lstrip(">")
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


def __is_merox_decoy(protein: str, decoy_prefix: str) -> bool:
    r"""Helper function to determine if a MeroX protein accession is a decoy.

    A protein is a decoy if its accession starts with ``decoy_prefix`` or with
    ``"DEC_"``, the prefix MeroX gives the decoy proteins it generates.

    Parameters
    ----------
    protein : str
        The protein accession as returned by ``__get_merox_protein()``.
    decoy_prefix : str
        The prefix that indicates that a protein is from the decoy database.

    Returns
    -------
    bool
        If the protein is a decoy.

    Notes
    -----
    This function should not be called directly, it is called from ``read_merox()``.
    """
    return protein.startswith(decoy_prefix) or protein.startswith(MEROX_DECOY_PREFIX)


def __is_merox_unplaced_decoy(
    proteins: str, from_position: Any, sequence: str, decoy_prefix: str
) -> bool:
    r"""Helper function to determine if a MeroX decoy peptide has no protein position.

    MeroX writes peptide start ``0`` (column "From") for a decoy peptide it could not
    place in its decoy protein, e.g. a short peptide like ``[VKR]``. Its protein
    positions are then unknown. A peptide starting with ``"{"`` is a protein N-terminal
    peptide, which MeroX also starts at ``0``.

    Parameters
    ----------
    proteins : str
        The protein string from MeroX e.g. from column "Protein 1".
    from_position : any
        The peptide start from MeroX e.g. from column "From".
    sequence : str
        The MeroX sequence e.g. from column "Peptide 1".
    decoy_prefix : str
        The prefix that indicates that a protein is from the decoy database.

    Returns
    -------
    bool
        If the peptide is a decoy without a protein position.

    Notes
    -----
    This function should not be called directly, it is called from ``read_merox()``.
    """
    return (
        __is_merox_decoy(__get_merox_protein(proteins)[0], decoy_prefix)
        and __parse_int(from_position) == 0
        and not str(sequence).strip().startswith("{")
    )


def __split_merox_scan_number(scan_nr_and_file: str) -> Tuple[Optional[str], str]:
    r"""Helper function to split the MeroX scan number string into scan number and spectrum file.

    Splits the MeroX scan number string into the scan number and the spectrum file. MeroX
    writes the following formats:

    - ``scan~file``, e.g. ``"10061~XLpeplib_Beveridge_QEx-HFX_DSS_R1"``.
    - ``scan`` without a spectrum file, e.g. ``"7755"``.
    - A spectrum title followed by the raw file and native ID, e.g.
      ``'run.10031.10031.4 File:"run.raw", NativeID:"controllerType=0 controllerNumber=1 scan=10031"'``.
      The scan number and the spectrum file are taken from the spectrum title.
      The spectrum title may be followed by annotations, e.g.
      ``'run.1633.1633.4 MSn:"707.85" RTmin:"16.78" File:"run.raw", NativeID:"... scan=1633"'``.
    - ``file|scan`` (StavroX), e.g. ``"3324_AJ_Dimer-Tryp|3953"``.
    - A spectrum title alone, e.g. ``"run.9005.9005.3"``.
    - A spectrum title of a SCIEX file, e.g. ``'Locus:1.1.1.1409.10 File:"run.wiff"'``. It
      contains no scan number: the scan number is ``None`` and the whole title, with its file,
      is returned as the spectrum file, as it identifies the spectrum.
    - A spectrum title of a Bruker timsTOF (PASEF) file, e.g.
      ``"11238 +MS2(650.0308), 32.6eV, 52.2min, 1/K0=0.850 #17538-17540|...|Id=11238|...~run"``.
      It gives the precursor (``Id``), the merged frames (``#``) and a mobility scan, but no
      scan number: as for SCIEX titles, the scan number is ``None`` and the whole title is the
      spectrum file.
    - No spectrum title (``"N/A"``): the scan number is ``None`` and the spectrum file is an
      empty string.

    The scan number may be a range of merged scans, e.g. ``"9512-9606"``, or a list of merged
    scans, e.g. ``"NS=2247,2252,2254-2257~run"``, of which the first scan is returned. A leading
    parenthesized annotation, e.g. ``"(-1) "``, a leading ``"NS="``, and a leading ``"sn"`` of
    the scan number, e.g. ``"NS=sn10318~run"``, are ignored.

    Parameters
    ----------
    scan_nr_and_file : str
        The scan number string from MeroX e.g. from the column "Scan number".

    Returns
    -------
    tuple of str or None, str
        The scan number and the spectrum file. The spectrum file is an empty string if
        it is not given. The scan number is ``None`` for a SCIEX or timsTOF spectrum title,
        and if no spectrum title is given.

    Raises
    ------
    IndexError
        If the spectrum title could not be parsed.

    Notes
    -----
    This function should not be called directly, it is called from ``read_merox()``.
    """
    value = re.sub(r"^\(.*?\)", "", str(scan_nr_and_file).strip()).strip()
    if value in ["", "nan", "N/A"]:
        return None, ""
    if value.startswith("Locus:") or re.match(r"[\d, ]+ \+MS2\(", value):
        return None, value
    value = re.sub(r"^NS=", "", value)
    if "~" in value:
        scan_nr, spectrum_file = value.split("~")[0], value.split("~")[1]
    elif re.fullmatch(r"(.+)\|(\d+)", value):
        spectrum_file, scan_nr = value.rsplit("|", 1)
    elif " File:" in value:
        # the title may be followed by annotations, e.g. 'MSn:"707.85" RTmin:"16.77"'
        title_str = re.sub(r'(\s+\w+:"[^"]*")+$', "", value.split(" File:")[0].strip())
        title = re.fullmatch(r"(.+)\.(\d+)\.\d+\.\d+", title_str)
        if title is None:
            raise IndexError(f"Could not parse spectrum title from {value}!")
        scan_nr, spectrum_file = title.group(2), title.group(1)
        native_scan = re.search(r'NativeID:"[^"]*\bscan=(\d+)', value)
        if native_scan is not None and int(native_scan.group(1)) != int(scan_nr):
            raise IndexError(
                f"Scan number of spectrum title and native ID differ in {value}!"
            )
    elif re.fullmatch(r"(.+)\.(\d+)\.\d+\.\d+", value):
        title = re.fullmatch(r"(.+)\.(\d+)\.\d+\.\d+", value)
        scan_nr, spectrum_file = title.group(2), title.group(1)
    else:
        scan_nr, spectrum_file = value, ""
    scan_nr = re.sub(r"^sn", "", scan_nr.split(",")[0].split("-")[0].strip())
    return scan_nr.strip(), spectrum_file.strip()


def __get_merox_scan_number(scan_nr_and_file: str) -> Optional[int]:
    r"""Helper function to parse the scan number from MeroX.

    Parses the scan number from the MeroX scan number string. Supports the formats
    described in ``__split_merox_scan_number()``, e.g. scan number and spectrum file
    delimted by the tilde (wave) sign. If the scan number is a range of merged scans
    the first scan number is returned. A SCIEX or timsTOF spectrum title contains no scan
    number, ``None`` is returned, as it is if no spectrum title is given.

    Parameters
    ----------
    scan_nr_and_file : str
        The scan number string from MeroX e.g. from the column "Scan number".

    Returns
    -------
    int, or None
        The parsed scan number, or None if the spectrum title contains no scan number.

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
        scan_nr_str = __split_merox_scan_number(scan_nr_and_file)[0]
        if scan_nr_str is not None:
            scan_nr = __parse_int(scan_nr_str)
    except Exception as e:
        raise RuntimeError(f"Could not scan number from {scan_nr_and_file}!") from e
    return scan_nr


def __get_merox_spectrum_file(scan_nr_and_file: str) -> str:
    r"""Helper function to parse the spectrum file from MeroX.

    Parses the spectrum file from the MeroX scan number string. Supports the formats
    described in ``__split_merox_scan_number()``, e.g. scan number and spectrum file
    delimted by the tilde (wave) sign. If the scan number string does not contain a
    spectrum file an empty string is returned.

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
        spectrum_file = __split_merox_scan_number(scan_nr_and_file)[1]
    except Exception as e:
        raise RuntimeError(
            f"Could not spectrum file name from {scan_nr_and_file}!"
        ) from e
    return spectrum_file


def __get_merox_retention_time(
    retention_time: Any, decimal: str = "."
) -> Optional[float]:
    r"""Helper function to parse the retention time from MeroX.

    Parses the retention time in seconds from MeroX. MeroX reports a range of retention
    times, e.g. ``"1854-1886"``, for merged scans, in that case the start of the range
    is returned. Some exported tables leave the retention time empty, e.g. ``" "``, in
    that case ``None`` is returned.

    Parameters
    ----------
    retention_time : any
        The retention time from MeroX e.g. from column "Retention time in sec".
    decimal : str, default = "."
        Character to recognize as decimal point.

    Returns
    -------
    float, or None
        The parsed retention time in seconds, or None if no retention time is given.

    Raises
    ------
    RuntimeError
        If the retention time could not be parsed.

    Notes
    -----
    This function should not be called directly, it is called from ``read_merox()``.
    """
    value = str(retention_time).strip()
    if value == "":
        return None
    # pandas does not apply the decimal character to columns that contain ranges
    if decimal != ".":
        value = value.replace(decimal, ".")
    rt = None
    try:
        rt_range = re.fullmatch(r"(\d+(?:\.\d+)?)-(\d+(?:\.\d+)?)", value)
        rt = __parse_float(rt_range.group(1) if rt_range is not None else value)
    except Exception as e:
        raise RuntimeError(
            f"Could not parse retention time from {retention_time}!"
        ) from e
    return rt


def __is_merox_unplaced(peptide_start: Any, peptide_end: Any) -> bool:
    r"""Helper function to determine if MeroX did not place a peptide in a protein.

    MeroX writes peptide start and end ``0`` (columns "From" and "To") for a peptide it did
    not place in a protein, e.g. one of several proteins in parentheses. A protein
    N-terminal peptide also starts at ``0``, but ends after its last residue.

    Parameters
    ----------
    peptide_start : any
        The peptide start from MeroX e.g. from column "From".
    peptide_end : any
        The peptide end from MeroX e.g. from column "To", or ``None`` if not known.

    Returns
    -------
    bool
        If the peptide has start and end ``0``.

    Notes
    -----
    This function should not be called directly, it is called from ``read_merox()``.
    """
    if peptide_end is None:
        return False
    return __parse_int(peptide_start) == 0 and __parse_int(peptide_end) == 0


def parse_site_from_merox(
    peptide: str,
    protein: str,
    peptide_start: Any,
    position: str,
    parse_modifications: bool = True,
    modifications: Dict[str, Dict[str, Any]] = MEROX_MODIFICATION_MAPPING,
    peptide_end: Any = None,
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
        ``constants.MEROX_MODIFICATION_MAPPING`` for examples. A symbol mapped to amino acid
        ``"{"`` is a modified protein N-terminus that MeroX writes after ``"["`` (see
        ``read_merox()``): the peptide is read as protein N-terminal, without the symbol.
    peptide_end : any, default = None
        The end of the peptide in the protein e.g. from column "To". MeroX writes start and
        end ``0`` for a peptide it did not place in a protein, which then has no protein
        positions.

    Returns
    -------
    dict of str, any
        The linked site, a dictionary with keys ``peptide`` (str), ``peptide_position`` (int,
        1-based, in-peptide), ``proteins`` (list of str, accessions) and ``protein_positions``
        (list of int, 1-based, absolute; parallel to ``proteins``). For a peptide shared by
        several proteins, or found several times in one (accessions joined by ``" and "``),
        the protein positions are ``None``: MeroX reports one position of the peptide only. They
        are also ``None`` for a peptide with start and end ``0`` (see ``peptide_end``).

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
    sequence, n_terminal_group = __get_merox_protein_n_terminal_group(
        str(peptide), modifications
    )
    peptide_position = __get_merox_group_crosslink_position(
        str(position), n_terminal_group
    )
    return {
        "peptide": __get_merox_sequence(sequence, parse_modifications, modifications),
        "peptide_position": peptide_position,
        "proteins": __get_merox_protein(str(protein)),
        "protein_positions": None
        if len(__split_merox_shared_proteins(str(protein))) > 1
        or __is_merox_unplaced(peptide_start, peptide_end)
        else [
            __get_merox_peptide_position(peptide_start, sequence) + peptide_position - 1
        ],
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
        The prefix that indicates that a protein is from the decoy database. Decoy
        proteins generated by MeroX itself (prefix ``"DEC_"``) are always recognised.
    parse_modifications : bool, default = True
        Whether or not post-translational-modifications should be parsed for crosslink-spectrum-matches.
        Requires correct specification of the 'modifications' parameter.
    modifications: dict of str, dict of str, any, default = ``constants.MEROX_MODIFICATION_MAPPING``
        Mapping of modification symbols to their amino acids and modifications. Please refer to
        ``constants.MEROX_MODIFICATION_MAPPING`` for examples. A symbol mapped to amino acid
        ``"{"`` is a modified protein N-terminus, a static modification of ``"{"`` in the
        MeroX settings (e.g. ``"{;X"`` for an acetylated N-terminus) that MeroX writes after
        ``"["`` (e.g. ``"[XMDVFMKGLSK]"``): the peptide is read as protein N-terminal, without
        the symbol, and with its modification at position 0.
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
    beta peptide, and ``Score`` as the score of the crosslink-spectrum-match. Tables exported
    by older versions of MeroX and by StavroX have no ``pepScore1`` and ``pepScore2``, the
    peptide scores are ``None`` for these.
    If MeroX reports a range of merged scans, the first scan number is used. If the scan
    number is reported without a spectrum file, the spectrum file is an empty string.
    Spectrum titles of SCIEX files, e.g. ``'Locus:1.1.1.1409.10 File:"run.wiff"'``, and of
    Bruker timsTOF files, e.g. ``"11238 +MS2(650.0308), ... #17538-17540|...|Id=11238|...~run"``,
    contain no scan number: the scan number is ``None`` and the spectrum file is the whole title,
    which identifies the spectrum. StavroX writes the scan number as ``file|scan``.
    If MeroX reports a range of retention times for merged scans, the start of the range
    is used. Tables without the column "Retention time in sec" are read without retention
    times (``None``).
    Monolinks (``Peptide2`` is ``0``) and looplinks (``Peptide2`` is ``1``) are skipped and their
    number is reported in a warning.
    A peptide shared by several proteins (accessions joined by ``" and "``, e.g.
    ``">sp|P0DMV9|HS71B_HUMAN and >sp|P0DMV8|HS71A_HUMAN@510-517"``) is read with all of its
    proteins and without protein positions, since MeroX reports the position in one of the
    proteins only; their number is reported in a warning.
    A peptide that MeroX did not place in a protein, with start and end ``0`` (columns "From"
    and "To"), is read without protein positions; their number is reported in a warning.

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
    skipped_monolinks = 0
    skipped_looplinks = 0
    unplaced_decoys = 0
    shared_peptides = 0
    unplaced_peptides = 0

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
            # skip monolinks ("0") and looplinks ("1"), pyXLMS only models crosslinks
            if str(row["Peptide2"]).strip() == "0":
                skipped_monolinks += 1
                continue
            if str(row["Peptide2"]).strip() == "1":
                skipped_looplinks += 1
                continue
            # modified protein N-termini, e.g. "[XMDVFMKGLSK]"
            peptide_1, n_terminal_group_1 = __get_merox_protein_n_terminal_group(
                str(row["Peptide 1"]), modifications
            )
            peptide_2, n_terminal_group_2 = __get_merox_protein_n_terminal_group(
                str(row["Peptide2"]), modifications
            )
            # parse linked sites
            site_a = parse_site_from_merox(
                row["Peptide 1"],
                row["Protein 1"],
                row["From"],
                row["best linkage position peptide 1"],
                parse_modifications,
                modifications,
                row["To"] if "To" in row.index else None,
            )
            site_b = parse_site_from_merox(
                row["Peptide2"],
                row["Protein 2"],
                row["From.1"],
                row["best linkage position peptide 2"],
                parse_modifications,
                modifications,
                row["To.1"] if "To.1" in row.index else None,
            )
            # decoy peptides that MeroX could not place in their decoy protein
            unplaced_a = __is_merox_unplaced_decoy(
                str(row["Protein 1"]), row["From"], peptide_1, decoy_prefix
            )
            unplaced_b = __is_merox_unplaced_decoy(
                str(row["Protein 2"]), row["From.1"], peptide_2, decoy_prefix
            )
            unplaced_decoys += int(unplaced_a) + int(unplaced_b)
            shared_peptides += int(
                len(__split_merox_shared_proteins(str(row["Protein 1"]))) > 1
            ) + int(len(__split_merox_shared_proteins(str(row["Protein 2"]))) > 1)
            unplaced_peptides += int(
                not unplaced_a and __is_merox_unplaced(row["From"], row.get("To"))
            ) + int(
                not unplaced_b and __is_merox_unplaced(row["From.1"], row.get("To.1"))
            )
            # a peptide at start 0 that is not protein N-terminal, also one without protein
            # positions (shared): the letter of a modified protein N-terminus is not mapped
            for peptide, start, end, unplaced in (
                (peptide_1, row["From"], row.get("To"), unplaced_a),
                (peptide_2, row["From.1"], row.get("To.1"), unplaced_b),
            ):
                if not unplaced and not __is_merox_unplaced(start, end):
                    __get_merox_placed_peptide_position(start, peptide)
            # create crosslink
            csm = create_csm(
                peptide_a=site_a["peptide"],
                modifications_a=__get_merox_modifications(
                    peptide_1,
                    __get_merox_position(str(row["best linkage position peptide 1"]))
                    if n_terminal_group_1 is None
                    else site_a["peptide_position"],
                    crosslinker,
                    crosslinker_mass,
                    modifications,
                    n_terminal_group_1,
                )
                if parse_modifications
                else None,
                xl_position_peptide_a=site_a["peptide_position"],
                proteins_a=site_a["proteins"],
                xl_position_proteins_a=None
                if unplaced_a
                else site_a["protein_positions"],
                pep_position_proteins_a=None
                if unplaced_a or site_a["protein_positions"] is None
                else [__get_merox_placed_peptide_position(row["From"], peptide_1)],
                score_a=__parse_float(row["pepScore1"])
                if "pepScore1" in row.index
                else None,
                decoy_a=__is_merox_decoy(site_a["proteins"][0], decoy_prefix),
                peptide_b=site_b["peptide"],
                modifications_b=__get_merox_modifications(
                    peptide_2,
                    __get_merox_position(str(row["best linkage position peptide 2"]))
                    if n_terminal_group_2 is None
                    else site_b["peptide_position"],
                    crosslinker,
                    crosslinker_mass,
                    modifications,
                    n_terminal_group_2,
                )
                if parse_modifications
                else None,
                xl_position_peptide_b=site_b["peptide_position"],
                proteins_b=site_b["proteins"],
                xl_position_proteins_b=None
                if unplaced_b
                else site_b["protein_positions"],
                pep_position_proteins_b=None
                if unplaced_b or site_b["protein_positions"] is None
                else [__get_merox_placed_peptide_position(row["From.1"], peptide_2)],
                score_b=__parse_float(row["pepScore2"])
                if "pepScore2" in row.index
                else None,
                decoy_b=__is_merox_decoy(site_b["proteins"][0], decoy_prefix),
                score=__parse_float(row["Score"]),
                spectrum_file=__get_merox_spectrum_file(str(row["Scan number"])),
                scan_nr=__get_merox_scan_number(str(row["Scan number"])),
                charge=__parse_int(row["Charge"]),
                rt=__get_merox_retention_time(row["Retention time in sec"], decimal)
                if "Retention time in sec" in row.index
                else None,
                im_cv=None,
                additional_information={
                    "source": __serialize_pandas_series(row),
                    "xLinkScore": row.get("xLinkScore"),
                    "Protein 1": row["Protein 1"],
                    "Protein 2": row["Protein 2"],
                    "MS1intensity": row.get("MS1intensity"),
                },
            )
            csms.append(csm)
    ## report skipped rows
    if skipped_monolinks + skipped_looplinks > 0:
        warnings.warn(
            RuntimeWarning(
                f"Skipped {skipped_monolinks} monolink and {skipped_looplinks} looplink "
                "rows, only crosslinks are read!"
            ),
            stacklevel=2,
        )
    ## report peptides not placed in a protein and shared by several proteins
    if unplaced_peptides > 0:
        warnings.warn(
            RuntimeWarning(
                f"Read {unplaced_peptides} peptide(s) that MeroX did not place in a protein "
                "(start and end 0) without protein positions!"
            ),
            stacklevel=2,
        )
    if shared_peptides > 0:
        warnings.warn(
            RuntimeWarning(
                f"Read {shared_peptides} peptide(s) shared by several proteins without protein "
                "positions, MeroX reports the position in one of the proteins only!"
            ),
            stacklevel=2,
        )
    ## report decoys without a position
    if unplaced_decoys > 0:
        warnings.warn(
            RuntimeWarning(
                f"Read {unplaced_decoys} decoy peptide(s) that MeroX did not place in "
                "their decoy protein (start 0) without protein positions!"
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
