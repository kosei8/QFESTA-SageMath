import importlib, sys, time
import parameter_generate as param
import elliptic_curve as ec
import QFESTA
import utilities_festa as utilities
importlib.reload(ec)
importlib.reload(param)
importlib.reload(QFESTA)

if __name__ == "__main__":

    utilities.speed_up_sagemath()

    N = 10
    for lam in [128, 192, 256]:
        KEM = QFESTA.QFESTA_KEM(lam)

        t_gen = 0
        t_enc = 0
        t_dec = 0
        for _ in range(N):
            t = time.time()
            sec_key, pub_key = KEM.Gen()
            t_gen += time.time() - t
            t = time.time()
            K, ciphertext = KEM.Encaps(pub_key)
            t_enc += time.time() - t
            t = time.time()
            Kd = KEM.Decaps(ciphertext, sec_key, pub_key)
            t_dec += time.time() - t
            assert K == Kd
        print("lam=%d: %.2f, %.2f, %.2f" % (lam, t_gen/N, t_enc/N, t_dec/N))
