import importlib, sys, time
import parameter_generate as param
import elliptic_curve as ec
import QFESTA
import utilities
importlib.reload(ec)
importlib.reload(param)
importlib.reload(QFESTA)

if __name__ == "__main__":
    args = sys.argv
    if len(args) > 1:
        lam = int(args[1])
    else:
        lam = 20

    utilities.speed_up_sagemath()

    t = time.time()
    KEM = QFESTA.QFESTA_KEM(lam)
    print("Set system parameter: lam=%d, a=%d, b=%d, k=%d, f=%d. %.2fsec."
          % (lam, KEM.a, KEM.b, KEM.k, KEM.f, time.time() - t))

    t = time.time()
    sec_key, pub_key = KEM.Gen()
    print("Keys are generated. Pubkey %d bytes. %.2fsec." % (len(pub_key), time.time() - t))

    t = time.time()
    K, ciphertext, d = KEM.Encaps(pub_key)
    print("Encaps. Ciphertext %d bytes. %.2fsec." % (len(ciphertext) + len(d), time.time() - t))

    t = time.time()
    Kd = KEM.Decaps(ciphertext, d, sec_key, pub_key)
    if K == Kd:
        print("Success Decaps. %.2fsec" % (time.time() - t))
    else:
        print("Failed!!!")