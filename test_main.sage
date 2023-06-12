import importlib, sys, time
import parameter_generate as param_gen
import elliptic_curve as ec
import newFESTA
importlib.reload(ec)
importlib.reload(param_gen)
importlib.reload(newFESTA)

if __name__ == "__main__":
    args = sys.argv
    if len(args) > 1:
        lam = int(args[1])
        if len(args) > 2:
            mod_c = float(args[2])
        else:
            mod_c = None
    else:
        lam = 20
    
    t = time.time()
    p, e2, e3, e5, f, D1, D2, Fp4, Fp2, i = param_gen.params(lam, mod_c)

    E0 = EllipticCurve(Fp2, [1, 0])
    E0.set_order((p+1)^2)
    basis2 = ec.basis(E0, 2, e2)
    basis3 = ec.basis(E0, 3, e3)

    sys_param = [Fp4, basis2, basis3, e2, e3, e5, D1, D2, i]
    print("Set system parameter: lam=%d, e2=%d, e3=%d, e5=%d, f=%d. %.2fsec." % (lam, e2, e3, e5, f, time.time() - t))

    t = time.time()
    sec_key, pub_key = newFESTA.key_gen(sys_param)
    print("Keys are generated. %.2fsec." % (time.time() - t))

    t = time.time()
    message = randint(0, 2^(e2-2))
    ciphertext = newFESTA.encryption(message, sys_param, pub_key)
    print("Encrypt a message: %d. %.2fsec." % (message, time.time() - t))

    t = time.time()
    m = newFESTA.decrypt(ciphertext, sys_param, sec_key)
    if message == m:
        print("Success decryption. %.2fsec" % (time.time() - t))
    else:
        print("Failed!!!")
