# A weighted reciprocal Lucas sum lies outside ℚ(√5)

For each fixed positive integer $a$, the sum $\sum_{n\ge1} n/L_{2an}$ lies outside $\mathbb Q(\sqrt5)$ and is therefore irrational. The paper proves this from a theorem for the weight-two value $\beta_{1/b}(2)$ at real algebraic integer bases greater than one whose other conjugates have modulus at most one.

[Read the paper (PDF)](output/pdf/weighted-even-lucas-reciprocals.pdf) · [LaTeX source](paper/weighted-even-lucas-reciprocals.tex)

From the repository root, run `python3 -B tools/build_paper.py` to build the paper or `python3 -B checks/check_three_poles.py` to reproduce the exact size-three calculation at both Lucas embeddings. [Requirements and instructions](README.txt) · [Expected output](checks/reference-results.json) · [Review scope](review/audit-scope.json)

The manuscript and explanatory text are available under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/). The software is distributed under [Apache-2.0](LICENSE).
