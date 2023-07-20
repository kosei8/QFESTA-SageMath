import importlib, sys, time
import parameter_generate as param
import elliptic_curve as ec
import newFESTA2
importlib.reload(ec)
importlib.reload(param)
importlib.reload(newFESTA2)

if __name__ == "__main__":
    args = sys.argv
    if len(args) > 1:
        lam = int(args[1])
    else:
        lam = 20

    t = time.time()
    KEM = newFESTA2.QFESTA_ROM(lam)
    print("Set system parameter: lam=%d, a=%d, b=%d, k=%d, f=%d. %.2fsec."
          % (lam, KEM.PKE.a, KEM.PKE.b, KEM.PKE.k, KEM.PKE.f, time.time() - t))

    t = time.time()
    sec_key, pub_key = KEM.Gen()
    print("Keys are generated. Pubkey %d bytes. %.2fsec." % (len(pub_key), time.time() - t))

    t = time.time()
    K, ciphertext = KEM.Encaps(pub_key)
    print("Encrypt a message. Ciphertext %d bytes. %.2fsec." % (len(ciphertext), time.time() - t))

    t = time.time()
    Kd = KEM.Decaps(ciphertext, sec_key, pub_key)
    if K == Kd:
        print("Success decryption. %.2fsec" % (time.time() - t))
    else:
        print("Failed!!!")