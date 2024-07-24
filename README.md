# Matching-Script-for-Twiss-Functions-Accelerator-Physics
This script will match twiss functions and magnetic elements for designing lattices. It uses linear matrix transport of Courant-Snyder parameters(3x3) linear matrix for matching drift lengths, quadrupole focusing strengths, solenoid focusing strengths as well as find the optimal twiss functions. It uses the following formalism: 

$$
\begin{equation}
\mathcal{M} = \begin{pmatrix}
M_{11} & M_{12} & M_{13} & M_{14} \\
M_{21} & M_{22} &M_{23} & M_{24} \\
M_{31} & M_{32} & M_{33} & M_{34} \\
M_{41} & M_{42} & M_{43} & M_{44}
\end{pmatrix}
\end{equation}
$$


The transformation of optic functions are given as:


$$
\begin{equation}
\begin{pmatrix}
\beta_{f} \\
\alpha_{f} \\
\gamma_{f}
\end{pmatrix}=\begin{pmatrix}
M_{11}^{2} & -2M_{11}M_{12} & M_{12}^{2} \\
-2M_{11}M_{21} & (M_{11}M_{22} + M_{12}M_{21}) & -2M_{12}M_{22} \\
M_{21}^{2} & -2M_{21}M_{22} & M_{22}^{2}
\end{pmatrix}\cdot \begin{pmatrix}
\beta_{i} \\
\alpha_{i} \\
\gamma_{i}
\end{pmatrix}
\end{equation}
$$

Where, $i,f$ stands for initial and final state. Using this approach we can find the optimized magnetic strengths of the elements. Similarly for Dispersion function we have:

$$
\begin{equation}
\begin{split}
D_{xf} &= M_{11}D_{xi} + M_{12}D_{xi}^{'} + \rho(1-\cos(l/\rho)) + M_{13}D_{yi} + M_{14}D_{yi}^{'}\\
D_{xf}' &= M_{21}D_{xi} + M_{22}D_{xi}^{'} + \sin(l/\rho) + M_{23}D_{yi} + M_{24}D_{yi}^{'}\\
D_{yf} &= M_{31}D_{xi} + M_{32}D_{xi}^{'} + M_{33}D_{yi} + M_{34}D_{yi}^{'}  \\
D_{yf}^{'} &= M_{41}D_{xi} + M_{42}D_{xi}^{'} + M_{43}D_{yi}^{'}
\end{split}
\end{equation}
$$


Where, $\rho$ is bending radius and $l$ is dipole length. Since we know how dispersion gets transported, this script can also be used for matching dispersion. Where the vertical dispersion is also added for solenoid case.

