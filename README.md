# QFESTA

A proof-of-concept implementation of the isogeny-based KEM
QFESTA (Quaternion Fast Encapsulation from Supersingular Torsion Attacks)
proposed in [QFESTA: Efficient Algorithms and Parameters for FESTA using Quaternion Algebas](https://eprint.iacr.org/2023/***.pdf)
using [SageMath](https://www.sagemath.org).

Our code is partially based on 
[FESTA-SageMath](https://github.com/FESTA-PKE/FESTA-SageMath/tree/main) [^1]
In particular, we use the FESTA-SageMath code for
the computation of isogenies between principally polarized abelian surfaces
and the key compression.
The following files are from FESTA-SageMath:
- compression.py
- divisor_arithmetic.py
- richelot_isogenies.py
- supersingular.py
- utilities.py

[^1]: committed on Jun 2, 2023; *commit id*: 7bc6c47eb3b87fd483be07fbbb4666174132d1a9.

## Usage
**Requirements**:
for the Fujisaki-Okamoto transformation,
we use SHAKE imported from pycryptodome to extract random bytes. This can be installed using:
```
pip install -r pycryptodome
```

You can execute QFESTA by the following command:
```
sage main.sage {secruty bits}
```
For example,
```
$ sage main.sage 128
Set system parameter: lam=128, a=272, b=162, k=131, f=169. 1.97sec.
ROM version
Keys are generated. Pubkey 174 bytes. 1.68sec.
Encaps. Ciphertext 348 bytes. 1.92sec.
Success Decaps. 2.71sec

QROM version
Keys are generated. Pubkey 174 bytes. 1.68sec.
Encaps. Ciphertext 382 bytes. 1.87sec.
Success Decaps. 2.67sec
```

Bench mark test for security bits 128, 192, and 256 is executed by
```
sage bench.sage
```

