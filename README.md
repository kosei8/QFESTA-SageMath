# QFESTA

A proof-of-concept implementation of the isogeny-based KEM
QFESTA (Quaternion Fast Encapsulation from Supersingular Torsion Attacks)
proposed in [QFESTA: Efficient Algorithms and Parameters for FESTA using Quaternion Algebas](https://eprint.iacr.org/2023/1468)
using [SageMath](https://www.sagemath.org).

Our code is partially based on 
[FESTA-SageMath](https://github.com/FESTA-PKE/FESTA-SageMath/tree/main) [^1]
by
[Basso-Maino-Pope](https://link.springer.com/chapter/10.1007/978-981-99-8739-9_4)
and
[Theta-Sagemath](https://github.com/ThetaIsogenies/two-isogenies)
by
[Dartois-Maino-Pope-Robert](https://eprint.iacr.org/2023/1747).

In particular, we use the FESTA-SageMath code for
the key compression
and
the computation of (2, 2)-isogenies
(the latter is optional).
The following files are from FESTA-SageMath:
- compression.py
- divisor_arithmetic.py
- richelot_isogenies.py
- supersingular.py
- utilities_festa.py

[^1]: committed on Jun 2, 2023; *commit id*: 7bc6c47eb3b87fd483be07fbbb4666174132d1a9.

We use theta model for the computation of (2, 2)-isogenies as default.
For this, we use the files in the following folders from Theta-SageMath:
- theta_isogenies
- theta_structure
- utilities

## Usage
**Requirements**:
for the Fujisaki-Okamoto transformation,
we use SHAKE imported from pycryptodome to extract random bytes. This can be installed using:
```
pip install -r pycryptodome
```

You can execute QFESTA by the following command:
```
sage main.sage {security bits} {"mum" if using Mumford model by FESTA-SageMath}
```
For example,
```
$ sage main.sage 128
Set system parameter: lam=128, a=390, b1=81, f=55. 9.16sec.
Keys are generated. Pubkey 247 bytes. 1.16sec.
Encaps. Ciphertext 494 bytes. 1.80sec.
Success Decaps. 4.72sec
```
```
$ sage main.sage 128 mum
Set system parameter: lam=128, a=390, b1=81, f=55. 9.67sec.
Keys are generated. Pubkey 247 bytes. 3.13sec.
Encaps. Ciphertext 494 bytes. 3.38sec.
Success Decaps. 7.75sec
```

Benchmark test for security bits 128, 192, and 256 is executed by
```
sage bench.sage
```

