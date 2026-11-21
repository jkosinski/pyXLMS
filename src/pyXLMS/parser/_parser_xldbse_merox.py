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


def __get_merox_protein(proteins: str) -> List[str]:
    r"""Helper function to extract the protein accession from MeroX.

    Parses the leading protein accession from the MeroX protein string. Additional
    proteins are not parsed as they do not have corresponding protein crosslink
    positions. MeroX writes the decoy proteins it generates as ``"DEC_>"`` followed by
    the accession of the target protein, e.g. ``"DEC_>Cas9"``; these are returned as
    ``"DEC_Cas9"``.

    Parameters
    ----------
    proteins : str
        The protein string from MeroX e.g. from column "Protein 1".

    Returns
    -------
    list of str
        A list containing the single parsed protein accession.

    Notes
    -----
    This function should not be called directly, it is called from ``read_merox()``.
    """
    protein = proteins.split("(>")[0].strip()
    if protein.startswith(MEROX_DECOY_PREFIX + ">"):
        protein = MEROX_DECOY_PREFIX + protein[len(MEROX_DECOY_PREFIX) + 1 :]
    return [protein.lstrip(">")]


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


def __split_merox_scan_number(scan_nr_and_file: str) -> Tuple[str, str]:
    r"""Helper function to split the MeroX scan number string into scan number and spectrum file.

    Splits the MeroX scan number string into the scan number and the spectrum file. MeroX
    writes the following formats:

    - ``scan~file``, e.g. ``"10061~XLpeplib_Beveridge_QEx-HFX_DSS_R1"``.
    - ``scan`` without a spectrum file, e.g. ``"7755"``.
    - A spectrum title followed by the raw file and native ID, e.g.
      ``'run.10031.10031.4 File:"run.raw", NativeID:"controllerType=0 controllerNumber=1 scan=10031"'``.
      The scan number and the spectrum file are taken from the spectrum title.

    The scan number may be a range of merged scans, e.g. ``"9512-9606"``, or a list of merged
    scans, e.g. ``"NS=2247,2252,2254-2257~run"``, of which the first scan is returned. A leading
    parenthesized annotation, e.g. ``"(-1) "``, a leading ``"NS="``, and a leading ``"sn"`` of
    the scan number, e.g. ``"NS=sn10318~run"``, are ignored. Spectrum titles of SCIEX files,
    e.g. ``'Locus:1.1.1.1409.10 File:"run.wiff"'``, contain no scan number and are not
    supported.

    Parameters
    ----------
    scan_nr_and_file : str
        The scan number string from MeroX e.g. from the column "Scan number".

    Returns
    -------
    tuple of str, str
        The scan number and the spectrum file. The spectrum file is an empty string if
        it is not given.

    Raises
    ------
    IndexError
        If the spectrum title could not be parsed.

    Notes
    -----
    This function should not be called directly, it is called from ``read_merox()``.
    """
    value = re.sub(r"^\(.*?\)", "", str(scan_nr_and_file).strip()).strip()
    value = re.sub(r"^NS=", "", value)
    if "~" in value:
        scan_nr, spectrum_file = value.split("~")[0], value.split("~")[1]
    elif " File:" in value:
        if value.startswith("Locus:"):
            raise IndexError(f"Could not parse scan number from SCIEX locus {value}!")
        title = re.fullmatch(r"(.+)\.(\d+)\.\d+\.\d+", value.split(" File:")[0].strip())
        if title is None:
            raise IndexError(f"Could not parse spectrum title from {value}!")
        scan_nr, spectrum_file = title.group(2), title.group(1)
    else:
        scan_nr, spectrum_file = value, ""
    scan_nr = re.sub(r"^sn", "", scan_nr.split(",")[0].split("-")[0].strip())
    return scan_nr.strip(), spectrum_file.strip()


def __get_merox_scan_number(scan_nr_and_file: str) -> int:
    r"""Helper function to parse the scan number from MeroX.

    Parses the scan number from the MeroX scan number string. Supports the formats
    described in ``__split_merox_scan_number()``, e.g. scan number and spectrum file
    delimted by the tilde (wave) sign. If the scan number is a range of merged scans
    the first scan number is returned.

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
        scan_nr = __parse_int(__split_merox_scan_number(scan_nr_and_file)[0])
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
        (list of int, 1-based, absolute; parallel to ``proteins``). MeroX reports a single
        protein per peptide.

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
    peptide_position = __get_merox_crosslink_position(str(position))
    return {
        "peptide": __get_merox_sequence(
            str(peptide), parse_modifications, modifications
        ),
        "peptide_position": peptide_position,
        "proteins": __get_merox_protein(str(protein)),
        "protein_positions": [
            __get_merox_peptide_position(peptide_start, str(peptide))
            + peptide_position
            - 1
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
    beta peptide, and ``Score`` as the score of the crosslink-spectrum-match. Tables exported
    by older versions of MeroX and by StavroX have no ``pepScore1`` and ``pepScore2``, the
    peptide scores are ``None`` for these.
    If MeroX reports a range of merged scans, the first scan number is used. If the scan
    number is reported without a spectrum file, the spectrum file is an empty string.
    If MeroX reports a range of retention times for merged scans, the start of the range
    is used.
    Monolinks (``Peptide2`` is ``0``) and looplinks (``Peptide2`` is ``1``) are skipped and their
    number is reported in a warning.

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
            # decoy peptides that MeroX could not place in their decoy protein
            unplaced_a = __is_merox_unplaced_decoy(
                str(row["Protein 1"]), row["From"], str(row["Peptide 1"]), decoy_prefix
            )
            unplaced_b = __is_merox_unplaced_decoy(
                str(row["Protein 2"]), row["From.1"], str(row["Peptide2"]), decoy_prefix
            )
            unplaced_decoys += int(unplaced_a) + int(unplaced_b)
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
                xl_position_proteins_a=None
                if unplaced_a
                else site_a["protein_positions"],
                pep_position_proteins_a=None
                if unplaced_a
                else [__get_merox_peptide_position(row["From"], str(row["Peptide 1"]))],
                score_a=__parse_float(row["pepScore1"])
                if "pepScore1" in row.index
                else None,
                decoy_a=__is_merox_decoy(site_a["proteins"][0], decoy_prefix),
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
                xl_position_proteins_b=None
                if unplaced_b
                else site_b["protein_positions"],
                pep_position_proteins_b=None
                if unplaced_b
                else [
                    __get_merox_peptide_position(row["From.1"], str(row["Peptide2"]))
                ],
                score_b=__parse_float(row["pepScore2"])
                if "pepScore2" in row.index
                else None,
                decoy_b=__is_merox_decoy(site_b["proteins"][0], decoy_prefix),
                score=__parse_float(row["Score"]),
                spectrum_file=__get_merox_spectrum_file(str(row["Scan number"])),
                scan_nr=__get_merox_scan_number(str(row["Scan number"])),
                charge=__parse_int(row["Charge"]),
                rt=__get_merox_retention_time(row["Retention time in sec"], decimal),
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
