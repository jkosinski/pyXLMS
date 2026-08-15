#!/usr/bin/env python3

# pyXLMS - DATA - MONO-LINK
# 2026 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

from __future__ import annotations

import sys
import copy
import numpy as np
from pydantic import BaseModel
from pydantic import Field
from pydantic import ConfigDict
from pydantic import computed_field

from ._util import check_input
from ._util import check_indexing

from typing import Optional
from typing import List
from typing import Dict
from typing import Tuple
from typing import Any

if sys.version_info >= (3, 8):
    from typing import Literal
else:  # pragma: no cover
    from typing_extensions import Literal

if sys.version_info >= (3, 9):
    from typing import Annotated
else:  # pragma: no cover
    from typing_extensions import Annotated


class MonoLink(BaseModel):
    r"""Class representing a mono-link (dead-end).

    A mono-link is a single modified residue where only one end of the crosslinker
    reacted with the peptide (the other end is hydrolyzed / a dead-end). Unlike a
    :class:`~pyXLMS.data.Crosslink` it has exactly one linked site and therefore no
    ``beta`` peptide.

    Parameters
    ----------
    peptide : str
        The unmodified amino acid sequence of the mono-linked peptide.
    peptide_crosslink_position : int
        The position of the crosslinker in the sequence of the peptide (1-based).
    proteins : list of str, or None, default = None
        The accessions of proteins that the peptide is associated with.
    proteins_crosslink_positions : list of int, or None, default = None
        Positions of the mono-link in the proteins (1-based).
    decoy : bool, or None, default = None
        Whether the peptide is from the decoy database or not.
    score : float, or None, default = None
        Score of the mono-link.
    additional_information : dict of str, any, or None, default = None
        A dictionary with additional information associated with the mono-link.

    Examples
    --------
    >>> from pyXLMS.data import MonoLink
    >>> ml = MonoLink(
    ...     peptide="PEKP",
    ...     peptide_crosslink_position=3,
    ... )
    """

    peptide: Annotated[
        str,
        Field(
            frozen=True,
            description="The unmodified amino acid sequence of the mono-linked peptide.",
        ),
    ]
    r"""
    The unmodified amino acid sequence of the mono-linked peptide. Amino acids should
    be in upper case. Modifications should not be included in the sequence.
    """
    peptide_crosslink_position: Annotated[
        int,
        Field(
            frozen=True,
            description="The position of the crosslinker in the sequence of the peptide (1-based).",
        ),
    ]
    r"""
    The position of the crosslinker in the sequence of the peptide (1-based).
    """
    proteins: Annotated[
        Optional[List[str]],
        Field(
            frozen=True,
            description="The accessions of proteins that the peptide is associated with.",
        ),
    ] = None
    r"""
    The accessions of proteins that the peptide is associated with.
    """
    proteins_crosslink_positions: Annotated[
        Optional[List[int]],
        Field(
            frozen=True,
            description="Positions of the mono-link in the proteins (1-based).",
        ),
    ] = None
    r"""
    Positions of the mono-link in the proteins (1-based). If given the list should be of
    the same length as ``proteins`` and the position at list index ``i`` should correspond
    to the protein at list index ``i`` in ``proteins``.
    """
    decoy: Annotated[
        Optional[bool],
        Field(
            frozen=True,
            description="Whether the peptide is from the decoy database or not.",
        ),
    ] = None
    r"""
    Whether the peptide is from the decoy database (``True``) or not (``False``).
    """
    score: Annotated[
        Optional[float], Field(frozen=True, description="Score of the mono-link.")
    ] = None
    r"""
    Score of the mono-link.
    """
    additional_information: Annotated[
        Optional[Dict[str, Any]],
        Field(
            frozen=False,
            description="A dictionary with additional information associated with the mono-link.",
        ),
    ] = None
    r"""
    A dictionary with additional information associated with the mono-link.
    """
    model_config = ConfigDict(
        validate_assignment=True, strict=True, str_strip_whitespace=True
    )
    r"""
    Pydantic configuration for the underlying validation model.
    """

    @computed_field(description="Data type of the object.")
    @property
    def data_type(self) -> Literal["mono-link"]:
        r"""
        Data type of the object.
        """
        return "mono-link"

    @computed_field(description="Completeness of the mono-link.")
    @property
    def completeness(self) -> Literal["full", "partial"]:
        r"""
        Completeness of the mono-link, e.g. ``"full"`` if all attributes
        are not ``None`` and else ``"partial"``.
        """
        full = all(
            [
                self.proteins is not None,
                self.proteins_crosslink_positions is not None,
                self.decoy is not None,
                self.score is not None,
            ]
        )
        return "full" if full else "partial"

    def model_post_init(self, context: Any = None) -> None:
        r"""
        Performs extra validation and post init functions.

        Warnings
        --------
        This method should not be called manually!
        """
        # extra validation
        if self.proteins is not None and self.proteins_crosslink_positions is not None:
            if len(self.proteins) != len(self.proteins_crosslink_positions):
                raise ValueError(
                    "Crosslink position has to be given for every protein! Length of proteins and proteins_crosslink_positions has to match!"
                )
        _ok = check_indexing(self.peptide_crosslink_position)
        _ok = (
            check_indexing(self.proteins_crosslink_positions)
            if self.proteins_crosslink_positions is not None
            else True
        )
        # normalize a NaN score to None
        if self.score is not None:
            if np.isnan(self.score):
                self.__dict__["score"] = None
        return

    def __getitem__(self, key: str) -> Any:
        r"""
        Support for dict-like access.
        """
        try:
            return getattr(self, key)
        except AttributeError as e:
            raise KeyError(f"'{key}' is not a valid field!") from e

    def __contains__(self, key: str) -> bool:
        r"""
        Support for ``in`` operator.
        """
        return hasattr(self, key)

    def items(self) -> List[Tuple[str, Any]]:
        r"""
        Support for dict-like read access for backward compatibility.
        """
        return self.model_dump(mode="python").items()

    def keys(self) -> List[str]:
        r"""
        Support for dict-like read access for backward compatibility.
        """
        return self.model_dump(mode="python").keys()  # ty: ignore[unsound-return-statement]

    def values(self) -> List[Any]:
        r"""
        Support for dict-like read access for backward compatibility.
        """
        return self.model_dump(mode="python").values()

    def copy_with_update(self, update: Dict[str, Any] = {}) -> MonoLink:
        r"""Creates a deep copy of the mono-link with optional attribute updates.

        Parameters
        ----------
        update : dict of str, any, default = empty dict
            Dictionary mapping attribute names (str) to their updated values.

        Returns
        -------
        MonoLink
            New mono-link with optionally updated attributes.
        """
        _ok = check_input(update, "update", dict)
        return MonoLink(
            peptide=self.peptide if "peptide" not in update else update["peptide"],
            peptide_crosslink_position=self.peptide_crosslink_position
            if "peptide_crosslink_position" not in update
            else update["peptide_crosslink_position"],
            proteins=copy.deepcopy(self.proteins)
            if "proteins" not in update
            else update["proteins"],
            proteins_crosslink_positions=copy.deepcopy(
                self.proteins_crosslink_positions
            )
            if "proteins_crosslink_positions" not in update
            else update["proteins_crosslink_positions"],
            decoy=self.decoy if "decoy" not in update else update["decoy"],
            score=self.score if "score" not in update else update["score"],
            additional_information=copy.deepcopy(self.additional_information)
            if "additional_information" not in update
            else update["additional_information"],
        )


def create_mono_link(
    peptide: str,
    xl_position_peptide: int,
    proteins: Optional[List[str]],
    xl_position_proteins: Optional[List[int]],
    decoy: Optional[bool],
    score: Optional[float],
    additional_information: Optional[Dict[str, Any]] = None,
) -> MonoLink:
    r"""Creates a mono-link data structure.

    Contains minimal data necessary for representing a single mono-link (dead-end).

    Parameters
    ----------
    peptide : str
        The unmodified amino acid sequence of the mono-linked peptide.
    xl_position_peptide : int
        The position of the crosslinker in the sequence of the peptide (1-based).
    proteins : list of str, or None
        The accessions of proteins that the peptide is associated with.
    xl_position_proteins : list of int, or None
        Positions of the mono-link in the proteins (1-based).
    decoy : bool, or None
        Whether the peptide is from the decoy database or not.
    score : float, or None
        Score of the mono-link.
    additional_information : dict with str keys, or None, default = None
        A dictionary with additional information associated with the mono-link.

    Returns
    -------
    MonoLink
        The data structure representing the mono-link.

    Raises
    ------
    TypeError
        If the parameter is not of the given class.
    ValueError
        If the length of crosslink positions is not equal to the length of proteins.

    Examples
    --------
    >>> from pyXLMS.data import create_mono_link
    >>> minimal_mono_link = create_mono_link(
    ...     peptide="PEPTIDE",
    ...     xl_position_peptide=1,
    ...     proteins=None,
    ...     xl_position_proteins=None,
    ...     decoy=None,
    ...     score=None,
    ... )
    """
    return MonoLink(
        peptide=peptide,
        peptide_crosslink_position=xl_position_peptide,
        proteins=proteins,
        proteins_crosslink_positions=xl_position_proteins,
        decoy=decoy,
        score=score,
        additional_information=additional_information,
    )


def create_mono_link_min(
    peptide: str,
    xl_position_peptide: int,
    **kwargs,
) -> MonoLink:
    r"""Creates a mono-link data structure from minimal input.

    Alias for ``data.create_mono_link()`` that sets all optional parameters to
    ``None`` for convenience.

    Parameters
    ----------
    peptide : str
        The unmodified amino acid sequence of the mono-linked peptide.
    xl_position_peptide : int
        The position of the crosslinker in the sequence of the peptide (1-based).
    **kwargs
        Any additional parameters will be passed to ``data.create_mono_link()``.

    Returns
    -------
    MonoLink
        The data structure representing the mono-link.

    Examples
    --------
    >>> from pyXLMS.data import create_mono_link_min
    >>> minimal_mono_link = create_mono_link_min("PEPTIDE", 1)
    """
    return create_mono_link(
        peptide=peptide,
        xl_position_peptide=xl_position_peptide,
        proteins=kwargs["proteins"] if "proteins" in kwargs else None,
        xl_position_proteins=kwargs["xl_position_proteins"]
        if "xl_position_proteins" in kwargs
        else None,
        decoy=kwargs["decoy"] if "decoy" in kwargs else None,
        score=kwargs["score"] if "score" in kwargs else None,
        additional_information=kwargs["additional_information"]
        if "additional_information" in kwargs
        else None,
    )
