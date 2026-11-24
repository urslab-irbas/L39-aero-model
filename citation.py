"""Citation of the L39-aero-model software.

If you use this code, its results or its polynomial aerodynamic model of the
L-39 aircraft in your work, please cite the article in which the model is
described:

    Madzharov, A. (2026). A six-degree-of-freedom flight-dynamics model of the L-39 aircraft with
    polynomial approximation of the aerodynamic characteristics. Aerospace Research in Bulgaria, XX
    (submitted).

Run  python citation.py  to print the reference in English and Bulgarian.
"""

CITATION_EN = (
    "Madzharov, A. (2026). A six-degree-of-freedom flight-dynamics model of the L-39 "
    "aircraft with polynomial approximation of the aerodynamic characteristics. "
    "Aerospace Research in Bulgaria, XX (submitted)."
)

CITATION_BG = (
    "Маджаров, А. (2026). Модел на динамиката на полета с шест степени на свобода на самолет L-39 "
    "с полиномна апроксимация на аеродинамичните характеристики. "
    "Aerospace Research in Bulgaria, XX (под печат)."
)

BIBTEX = """@article{Madzharov2026L39,
  author  = {Madzharov, Anastas},
  title   = {A six-degree-of-freedom flight-dynamics model of the L-39 aircraft with
             polynomial approximation of the aerodynamic characteristics},
  journal = {Aerospace Research in Bulgaria},
  volume  = {XX},
  year    = {2026},
  note    = {Submitted}
}"""


def print_citation():
    """Print the reference that should be cited when this software is used."""
    print()
    print("If you use this software, please cite:")
    print("  " + CITATION_EN)
    print("При използване на програмата, моля, цитирайте:")
    print("  " + CITATION_BG)
    print()


def handle_help(doc):
    """Print the module help and the citation for -h / --help, then exit."""
    import sys
    if any(a in ("-h", "--help") for a in sys.argv[1:]):
        print(doc)
        print_citation()
        sys.exit(0)


if __name__ == "__main__":
    print_citation()
    print(BIBTEX)
