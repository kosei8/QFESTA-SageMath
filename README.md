# QFESTA

A proof-of-concept implementation of the isogeny-based PKE/KEM QFESTA
(Quaternion Fast Encryption/Encapsulation from Supersingular Torsion Attacks),
using [SageMath](https://www.sagemath.org).

Our code is partially based on 
[FESTA-SageMath](https://github.com/FESTA-PKE/FESTA-SageMath/tree/main).
In particular, we use the FESTA-SageMath code for
the computation of isogenies between principally polarized abelian surfaces
and the key compression.
The following files are from FESTA-SageMath:
- compression.py
- divisor_arithmetic.py
- richelot_isogenies.py
- supersingular.py
- utilities.py
