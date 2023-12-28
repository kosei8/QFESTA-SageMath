import importlib, sys, time
import parameter_generate as param
import elliptic_curve as ec
import QFESTA
import utilities_festa as utilities
importlib.reload(ec)
importlib.reload(param)
importlib.reload(QFESTA)

if __name__ == "__main__":
    args = sys.argv
    use_theta = False
    if len(args) > 1:
        lam = int(args[1])
        if len(args) > 2 and args[2] == "theta":
            use_theta = True
    else:
        lam = 20

    sys.setrecursionlimit(2000) # for using strategy for lam=256 in theta isogenies
    utilities.speed_up_sagemath()

    t = time.time()
    KEM = QFESTA.QFESTA_KEM(lam, use_theta)
    print("Set system parameter: lam=%d, a=%d, b1=%d, f=%d. %.2fsec."
          % (lam, KEM.a, KEM.b1, KEM.f, time.time() - t))

    t = time.time()
    sec_key, pub_key = KEM.Gen()
    print("Keys are generated. Pubkey %d bytes. %.2fsec." % (len(pub_key), time.time() - t))

    t = time.time()
    K, ciphertext = KEM.Encaps(pub_key)
    print("Encaps. Ciphertext %d bytes. %.2fsec." % (len(ciphertext), time.time() - t))

    t = time.time()
    Kd = KEM.Decaps(ciphertext, sec_key, pub_key)
    if K == Kd:
        print("Success Decaps. %.2fsec" % (time.time() - t))
    else:
        print("Failed!!!")