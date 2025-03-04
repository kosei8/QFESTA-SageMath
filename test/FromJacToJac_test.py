import cProfile
import pstats
import time
import sys
import os

# .env ファイルのパス
env_path = '.env'

# .env ファイルが存在するか確認し、存在する場合は読み込む
if os.path.exists(env_path):
    with open(env_path) as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue  # コメント行または空行をスキップ
            # 環境変数名と値を分割
            key, value = line.strip().split('=', 1)
            os.environ[key] = value  # 環境変数として設定

# 環境変数を使用
path = os.getenv('PATH')

sys.path.append(path)

# Sage Imports
from sage.all import EllipticCurve, GF, ceil, Matrix, ZZ, Zmod, inverse_mod, PolynomialRing

import richelot_isogenies as ri
import utilities.strategy as us
import utilities_festa as uf
import parameter_generate as pm
import elliptic_curve as ec

# ========================= #
#       　   SetUp    　　   #
# ========================= #

# Sage goes vroom!
uf.speed_up_sagemath()

# Default is FESTA_128
SECURITY = "128"

for arg in sys.argv[1:]:
    if arg.lower() in ["--toy", "-t"]:
        SECURITY = "TOY"
    elif arg.lower() in ["--192", "-II"]:
        SECURITY = "192"
    elif arg.lower() in ["--256", "-V"]:
        SECURITY = "256"

# PKE = FESTA
NAME = "QFESTA_" + SECURITY

uf.print_info(f"FromJacToJac Benchmarking {NAME}")
N_Enc = 10000

# set variables
a, b1, b2, f, D1, D2 = pm.SysParam2(int(SECURITY))
p = ZZ(2**a*3*f - 1)
Fp4, Fp2, zeta2 = pm.calcFields(p)

# supersingular elliptic curve E0
Fp2d = GF(p**2, modulus=[1, 0, 1], name="i")
i = Fp2d.gen()
E0 = EllipticCurve(Fp2, [1, 0])
basis2 = ec.basis(E0, 2, a)

# polynomial ring (FromProdToJacで定義されているやつ)
R = PolynomialRing(Fp2d, name="x", implementation="generic")
x = R.gens()[0]
# print(R)
h1, h3, h5 = 0, 0, 0
h0 = 331097855532447697541342307488045815138036440508183133213378496112289420149639990579806419396306486268727010878064653912*i + 140687295187471608261251439694691449581467911545259174431463943959908804617300774770269208381297093393334489626336290844
h2 = 234323807688032733075911540497420562791608660780024524827179874493000568641197702157559957182732337184421912737973688530*i + 86792300056800013181711652320027725549433831141989187034137725760860959359837842642493722628247022180316426256281384657
h4 = 388383452938669061149330672249408011900460325114916241110823683537242333806561789270771636299068370845459022400514190660*i + 297382011647898783423909184312560746341578391037545089162804587715411796691721212962304652822831855307104944090771310145
h6 = 180656699629745285321781431310827035740597172647725440029110025831007248391675250892679055123633374257205747974665813657*i + 100800524209352790314138249857650760902473005292396598966198354254886806508429215511104928984755170690097838778124938083
h_coeffs = [h0,h1,h2,h3,h4,h5,h6]
# print(h_coeffs)
h = R(h_coeffs)


# point
# 128
D1 = (x**2 + (331179519723148188837407421367720507986636386865520153779319236344137997520951231187996303105265849386013731256562949541*i + 266575282061573154374651002342964059446020076195900643270683758128313089839477532223505753339807425050743716819558709321)*x + 1, 0)
D2 = (x**2 + (82672520745883697743190351557808436212681405913850481482011682328209691785650314446631665517736696673943053733809261386*i + 140509532932515728001147123832308703565262584802540271408387553308359819235108815774713733646728724737018924322310145167)*x + 380727083029198096828453705777308579080975505671718848471144022568634758250239690018757914895137039753697396743069410345*i + 233872170612862199922682912663082126801159545409392413855945368958743681339358945042118869853614335671382684793638713570, 0)
# 192
# 256


# ========================= #
#      Start Profiling      #
# ========================= #
# Start the profiler
setup_time = time.time()
pr = cProfile.Profile()
pr.enable()

# ===================== #
#         Main          #
# ===================== #
for _ in range(N_Enc):
    h_new, f = ri.FromJacToJac(h, D1, D2)
    # print(h_new[6])
    # print(h_new[5])
    # print(h_new[4])
    # print(h_new[3])
    # print(h_new[2])
    # print(h_new[1])
    # print(h_new[0])

# h, f = ri.FromJacToJac(h, D1, D2)

# ========================= #
#       End Profiling       #
# ========================= #

pr.disable()
pr.dump_stats("festa_keygen.cProfile")
uf.print_info(f"FromJacToJac took: {(time.time() -  setup_time):.3f} seconds")
p = pstats.Stats("festa_keygen.cProfile")
p.strip_dirs().sort_stats("cumtime").print_stats(50)