# Ising spins with Hadamard couplings

**Opus and Teo**, November 2026

A real Hadamard matrix H of order N defines a system of 2N Ising spins, carried by 2 layers of crossed wires, with couplings H<sub>jk</sub>/&radic;N between the layers. This is the switching game of Gale and Berlekamp at finite temperature, and having symmetric orthogonal couplings, it is a deterministic instance of the random orthogonal model of Marinari, Parisi and Ritort. The density of states is the real glow of H, the law of &Omega; = &lang;Hb,a&rang; over the sign vectors a, b. The main results are as follows:

- **Modulo 16.** The law of &Omega; modulo 16 is determined by N and the 2-rank r of H: p<sub>0</sub> &minus; p<sub>8</sub> = (&minus;1)<sup>N/8</sup>(2<sup>&minus;N/2</sup> + 2<sup>&minus;2&minus;r</sup>) for N = 0(8), and (&minus;1)<sup>(N&minus;4)/8</sup>2<sup>&minus;N/2&minus;1</sup> for N = 4(8).
- **Moments.** The moments of &Omega; are universal polynomials in N up to order 14, the first correction, at order 16, being given by the 4-profile of H.
- **Free energy.** The high-temperature expansion of the free energy is universal to leading order, given by a Haar orthogonal computation, and the continued paramagnetic entropy vanishes at &beta; = &radic;240.
- **Order 16.** The 5 densities of states form a 1-parameter family, the parameter being the number 2<sup>14&minus;r</sup> &minus; 128 of crystalline ground states, and the frustrated matrix had<sub>16.3</sub> has 43008 ground states, forming 3 orbits under its symmetry group, with free spins in each of them.
- **Larger orders.** At N = 20 the 3 densities of states coincide. At N = 24 they form a 2-parameter family, 2 arrays having the same density of states precisely when their matrices have the same permanent, up to sign. At N = 28 they form a 4-parameter family, and determine the matrix up to equivalence and transposition. At N = 32, for 360 matrices, they are distinct and lie in an affine space of dimension 8. At N = 36, for 9 regular matrices, the crystallization temperature predicts the maxima of the specific heat within 0.0025.
- **Dynamics.** Metropolis simulations show a crystallization for W<sub>64</sub>, W<sub>256</sub>, while for W<sub>1024</sub>, W<sub>4096</sub> the crystal is missed, the cooling ending at about the same energy as for the Paley matrices of order 252 and 1020.

## Contents

- `021.pdf`, `021.tex`: the paper.
- The programs below, with their outputs (`.log`, `.json`, `.dat`), and `runall.sh`, which reruns them.
- `matrices/`: the matrices of order N &le; 20, from McKay's catalogue, together with the Paley matrix of order 12 and the Walsh matrices.
- `mckay/`: McKay's files of the Hadamard matrices of order 24 and 28, converted into `matrices24/`, `matrices28/` by `mckay.py`.
- `matrices32/`, `matrices36/`: the matrices of order 32 and 36, made by `make32.py`, `make36.py`.
- `matrices24T/`, `matrices28T/`: the transposed matrices, made by `transposes.py`.
- `matrices24_raw/`, `matrices24T_raw/`, `matrices28T_raw/`, and the `_output/` folders: the outputs of `dos.c`, namely the distributions of the multisets of numbers |(Hb)<sub>j</sub>|.
- `hypercube/`: the weighing matrices of the fully frustrated hypercube, for D = 2, ..., 6.
- `dyn/`: the matrices and crystalline states used by `dynamics.py`.

## Programs

| File | Content | Section |
|------|---------|---------|
| `hypercube.py` | The fully frustrated hypercube as the Ising array of a weighing matrix, and its ground states | 2 |
| `hypercubedos.py` | Exact densities of states of the hypercube, for D &le; 6 | 2 |
| `glowlaw.py` | Exact densities of states for N &le; 20 | 3 |
| `twoadic.py` | The residues modulo 16, and the partition function at imaginary temperature | 3 |
| `residues.py` | The law modulo 16 at all the orders, against the 2-rank | 3 |
| `moments.py` | Exact moments, the 4-profile, and the order 16 formula | 4 |
| `universal.py` | The universal polynomials giving the moments and cumulants | 4 |
| `groundstates.py` | Ground states, numbers of lit lamps, free spins | 6 |
| `autgroup.py` | The symmetry group of had<sub>16.3</sub>, of order 86016 | 6 |
| `crystal.py` | Crystallization temperatures, bent functions, comparison with the specific heat | 6, 7, 8 |
| `spectrum16.py` | The 1-parameter family of densities of states of order 16 | 7 |
| `thermo.py` | Exact thermodynamics at N = 16, 20 | 7 |
| `figures.py` | Coordinates of the plots of Section 7 | 7 |
| `multisets.py` | The coincidence at N = 20, at the level of the multisets | 7 |
| `pictures.py` | Data of the pictures and tables of Sections 2 to 8 | 2&ndash;8 |
| `dynamics.py` | Metropolis dynamics of W<sub>64</sub>, W<sub>256</sub> and of the Paley matrix of order 252 | 8 |
| `dynlarge.py` | Summary of the large simulations of `stipe.py` | 8 |
| `mckay.py` | Conversion of McKay's files of orders 24, 28 | 8 |
| `dosall.py` | Exact densities of states of all the matrices of a folder | 8 |
| `structure.py` | Structure of the families of densities of states, first differing moments | 8 |
| `transposes.py`, `transpose.py` | Transposition, and the invariants used to test it | 8 |
| `permanent.py` | Permanents of the matrices of order 24 and 28, against their densities of states | 8 |
| `make32.py`, `make36.py` | The matrices of order 32 and 36 | 8 |
| `dos.c` | Enumeration of the vectors b in Gray code order, counting the multisets of numbers \|(Hb)<sub>j</sub>\| | 2, 8 |
| `mc.c` | Metropolis dynamics, with single spin flips | 8 |
| `perm.c` | Permanents, by the formula of Glynn | 8 |

## Laptop programs

| File | Content | Output |
|------|---------|--------|
| `gaethje.py` | Densities of states of the 487 matrices of order 28 | `gaethje_output/` |
| `golovkin.py` | Densities of states of half of the 360 matrices of order 32 | `golovkin_output/` |
| `valentina.py` | Densities of states of the other half of the matrices of order 32 | `valentina_output/` |
| `frazier.py` | Densities of states of the 9 matrices of order 36 | `frazier_output/` |
| `stipe/stipe.py` | Metropolis dynamics of W<sub>1024</sub>, W<sub>4096</sub> and of the Paley matrix of order 1020 | `stipe/stipe_output/` |

## Running the code

The programs need Python 3 with numpy, and a C compiler. The command `sh runall.sh` reruns all the programs, in about 15 minutes, except for the laptop programs and the permanents of order 28, whose outputs are included. Each laptop program runs on its own, with `python3 gaethje.py` and so on, in about 1 hour, printing its progress and writing its results in its output folder.
