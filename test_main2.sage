import importlib, sys, time
import parameter_generate as param
import parameter_generate_new as param_new
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
    a, b, f, k, D1, D2 = param_new.SysParam(lam)
    p = 2**a*3*f - 1
    Fp4, Fp2, i = param.calcFields(p)

    E0 = EllipticCurve(Fp2, [1, 0])
    E0.set_order((p+1)^2)
    basis2 = ec.basis(E0, 2, a)

    sys_param = [Fp4, basis2, a, b, k, D1, D2, i]
    print("Set system parameter: lam=%d, a=%d, b=%d, k=%d, f=%d. %.2fsec." % (lam, a, b, k, f, time.time() - t))

    t = time.time()
    sec_key, pub_key = newFESTA2.key_gen(sys_param)
    print("Keys are generated. %.2fsec." % (time.time() - t))

    t = time.time()
    message = randint(0, 2^(a-2))
    ciphertext = newFESTA2.encrypt(message, sys_param, pub_key)
    print("Encrypt a message: %d. %.2fsec." % (message, time.time() - t))

    t = time.time()
    m = newFESTA2.decrypt(ciphertext, sys_param, sec_key, pub_key)
    if message == m:
        print("Success decryption. %.2fsec" % (time.time() - t))
    else:
        print("Failed!!!")
