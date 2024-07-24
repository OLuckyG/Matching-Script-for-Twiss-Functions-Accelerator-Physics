# Matching-Script-for-Twiss-Functions-Accelerator-Physics
This script will match twiss functions and magnetic elements for designing lattices. It uses linear matrix transport of Courant-Snyder parameters(3x3) linear matrix for matching drift lengths, quadrupole focusing strengths, solenoid focusing strengths as well as find the optimal twiss functions. It uses the following formalism:
$$
\begin{equation}
\mathcal{M} = \begin{pmatrix}
M_{11} & M_{12} & 0 & 0 \\
M_{21} & M_{22} &0 & 0 \\
0 & 0 & M_{33} & M_{34} \\
0 & 0 & M_{43} & M_{44}
\end{pmatrix}
\end{equation}
$$
