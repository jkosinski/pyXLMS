#!/usr/bin/env python3

# 2025 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

from __future__ import annotations

import io
import re
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
    positions.

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
    return [proteins.split("(>")[0].strip().lstrip(">")]


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
    If MeroX reports a range of merged scans, the first scan number is used. If the scan
    number is reported without a spectrum file, the spectrum file is an empty string.
    Spectrum titles of SCIEX files, e.g. ``'Locus:1.1.1.1409.10 File:"run.wiff"'``, and of
    Bruker timsTOF files, e.g. ``"11238 +MS2(650.0308), ... #17538-17540|...|Id=11238|...~run"``,
    contain no scan number: the scan number is ``None`` and the spectrum file is the whole title,
    which identifies the spectrum. StavroX writes the scan number as ``file|scan``.

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
            # create crosslink
            csm = create_csm(
                peptide_a=__get_merox_sequence(
                    str(row["Peptide 1"]), parse_modifications, modifications
                ),
                modifications_a=__get_merox_modifications(
                    str(row["Peptide 1"]),
                    __get_merox_position(str(row["best linkage position peptide 1"])),
                    crosslinker,
                    crosslinker_mass,
                    modifications,
                )
                if parse_modifications
                else None,
                xl_position_peptide_a=__get_merox_crosslink_position(
                    str(row["best linkage position peptide 1"])
                ),
                proteins_a=__get_merox_protein(str(row["Protein 1"])),
                xl_position_proteins_a=[
                    __get_merox_peptide_position(row["From"], str(row["Peptide 1"]))
                    + __get_merox_crosslink_position(
                        str(row["best linkage position peptide 1"])
                    )
                    - 1
                ],
                pep_position_proteins_a=[
                    __get_merox_peptide_position(row["From"], str(row["Peptide 1"]))
                ],
                score_a=__parse_float(row["pepScore1"]),
                decoy_a=__get_merox_protein(str(row["Protein 1"]))[0].startswith(
                    decoy_prefix
                ),
                peptide_b=__get_merox_sequence(
                    str(row["Peptide2"]), parse_modifications, modifications
                ),
                modifications_b=__get_merox_modifications(
                    str(row["Peptide2"]),
                    __get_merox_position(str(row["best linkage position peptide 2"])),
                    crosslinker,
                    crosslinker_mass,
                    modifications,
                )
                if parse_modifications
                else None,
                xl_position_peptide_b=__get_merox_crosslink_position(
                    str(row["best linkage position peptide 2"])
                ),
                proteins_b=__get_merox_protein(str(row["Protein 2"])),
                xl_position_proteins_b=[
                    __get_merox_peptide_position(row["From.1"], str(row["Peptide2"]))
                    + __get_merox_crosslink_position(
                        str(row["best linkage position peptide 2"])
                    )
                    - 1
                ],
                pep_position_proteins_b=[
                    __get_merox_peptide_position(row["From.1"], str(row["Peptide2"]))
                ],
                score_b=__parse_float(row["pepScore2"]),
                decoy_b=__get_merox_protein(str(row["Protein 2"]))[0].startswith(
                    decoy_prefix
                ),
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
