Weighted reciprocals of the even Lucas numbers

For every fixed positive integer a, the sum of n/L_(2an) over n >= 1
lies outside Q(sqrt(5)). At a = 1 it also equals the alternating sum
of the squared reciprocals of the odd-indexed Lucas numbers.

More generally, let b > 1 be a real algebraic integer whose other
conjugates have modulus at most one. The manuscript proves that the
weight-two value beta_(1/b)(2) lies outside Q(b) by estimating an
integer polynomial at every embedding of that field.

Read the paper at output/pdf/weighted-even-lucas-reciprocals.pdf.
The editable source and bibliography are in paper/weighted-even-lucas-reciprocals.tex.

Run both commands from the directory containing this README.

    python3 tools/build_paper.py
    python3 checks/check_three_poles.py

The build requires Python 3, pdfLaTeX with the Latin Modern fonts,
AMS-LaTeX, mathtools, geometry, microtype, needspace, xcolor and hyperref, and
Poppler's pdfinfo. It writes the PDF and paper/build-record.json after
checking the final LaTeX log.

The checker requires Python 3 and SymPy. At determinant size three, it
compares polynomial division with the pole recurrence and checks all
four normalized coefficients, their field conjugates and integrality,
and the sign of the leading coefficient. It uses exact arithmetic.

