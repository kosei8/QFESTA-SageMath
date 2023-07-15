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
    sys_param = newFESTA2.setup(lam)
    print("Set system parameter: lam=%d, a=%d, b=%d, k=%d, f=%d. %.2fsec."
          % (lam, sys_param["a"], sys_param["b"], sys_param["k"], sys_param["f"], time.time() - t))

    t = time.time()
    sec_key, pub_key = newFESTA2.key_gen(sys_param)
    print("Keys are generated. Pubkey %d bytes. %.2fsec." % (len(pub_key), time.time() - t))

    t = time.time()
    message = randint(0, 2^(sys_param["a"]-2))
    ciphertext = newFESTA2.encrypt(message, sys_param, pub_key)
    print("Encrypt a message. Ciphertext %d bytes. %.2fsec." % (len(ciphertext), time.time() - t))

    t = time.time()
    m = newFESTA2.decrypt(ciphertext, sys_param, sec_key, pub_key)
    if message == m:
        print("Success decryption. %.2fsec" % (time.time() - t))
    else:
        print("Failed!!!")