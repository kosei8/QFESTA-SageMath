import cProfile
import pstats
import time
import sys
sys.path.append('/Users/konami/utokyo/00_lab/software/QFESTA-SageMath')

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

uf.print_info(f"(2,2)-Isogeny Benchmarking {NAME}")

# set variables
a, b1, b2, f, D1, D2 = pm.SysParam2(int(SECURITY))
p = ZZ(2**a*3*f - 1)
Fp4, Fp2, zeta2 = pm.calcFields(p)

# supersingular elliptic curve E0
Fp2d = GF(p**2, modulus=[1, 0, 1], name="i")
i = Fp2d.gen()
E0 = EllipticCurve(Fp2, [1, 0])
basis2 = ec.basis(E0, 2, a)

# For magic square root in splitting
# 128
b = 390
# 192
# b = 582
# 256
# b = 774


# strategy
strategy = us.optimised_strategy(b - 1)

# 256
# print(strategy)
# strategy = [589, 141, 34, 7, 1, 6, 5, 4, 3, 2, 1, 27, 6, 5, 4, 3, 2, 1, 21, 5, 4, 3, 2, 1, 16, 4, 3, 2, 1, 12, 3, 2, 1, 9, 2, 1, 7, 1, 6, 5, 4, 3, 2, 1, 107, 27, 6, 5, 4, 3, 2, 1, 21, 5, 4, 3, 2, 1, 16, 4, 3, 2, 1, 12, 3, 2, 1, 9, 2, 1, 7, 1, 6, 5, 4, 3, 2, 1, 81, 20, 5, 4, 3, 2, 1, 15, 4, 3, 2, 1, 11, 3, 2, 1, 8, 2, 1, 6, 1, 5, 4, 3, 2, 1, 61, 15, 4, 3, 2, 1, 11, 3, 2, 1, 8, 2, 1, 6, 1, 5, 4, 3, 2, 1, 46, 11, 3, 2, 1, 8, 2, 1, 6, 1, 5, 4, 3, 2, 1, 35, 8, 2, 1, 6, 1, 5, 4, 3, 2, 1, 27, 6, 1, 5, 4, 3, 2, 1, 21, 5, 4, 3, 2, 1, 16, 4, 3, 2, 1, 12, 3, 2, 1, 9, 2, 1, 7, 1, 6, 5, 4, 3, 2, 1, 449, 107, 26, 6, 5, 4, 3, 2, 1, 20, 5, 4, 3, 2, 1, 15, 4, 3, 2, 1, 11, 3, 2, 1, 8, 2, 1, 6, 1, 5, 4, 3, 2, 1, 81, 20, 5, 4, 3, 2, 1, 15, 4, 3, 2, 1, 11, 3, 2, 1, 8, 2, 1, 6, 1, 5, 4, 3, 2, 1, 61, 15, 4, 3, 2, 1, 11, 3, 2, 1, 8, 2, 1, 6, 1, 5, 4, 3, 2, 1, 46, 11, 3, 2, 1, 8, 2, 1, 6, 1, 5, 4, 3, 2, 1, 35, 8, 2, 1, 6, 1, 5, 4, 3, 2, 1, 27, 6, 1, 5, 4, 3, 2, 1, 21, 5, 4, 3, 2, 1, 16, 4, 3, 2, 1, 12, 3, 2, 1, 9, 2, 1, 7, 1, 6, 5, 4, 3, 2, 1, 342, 81, 20, 5, 4, 3, 2, 1, 15, 4, 3, 2, 1, 11, 3, 2, 1, 8, 2, 1, 6, 1, 5, 4, 3, 2, 1, 61, 15, 4, 3, 2, 1, 11, 3, 2, 1, 8, 2, 1, 6, 1, 5, 4, 3, 2, 1, 46, 11, 3, 2, 1, 8, 2, 1, 6, 1, 5, 4, 3, 2, 1, 35, 8, 2, 1, 6, 1, 5, 4, 3, 2, 1, 27, 6, 1, 5, 4, 3, 2, 1, 21, 5, 4, 3, 2, 1, 16, 4, 3, 2, 1, 12, 3, 2, 1, 9, 2, 1, 7, 1, 6, 5, 4, 3, 2, 1, 261, 61, 15, 4, 3, 2, 1, 11, 3, 2, 1, 8, 2, 1, 6, 1, 5, 4, 3, 2, 1, 46, 11, 3, 2, 1, 8, 2, 1, 6, 1, 5, 4, 3, 2, 1, 35, 8, 2, 1, 6, 1, 5, 4, 3, 2, 1, 27, 6, 1, 5, 4, 3, 2, 1, 21, 5, 4, 3, 2, 1, 16, 4, 3, 2, 1, 12, 3, 2, 1, 9, 2, 1, 7, 1, 6, 5, 4, 3, 2, 1, 200, 46, 11, 3, 2, 1, 8, 2, 1, 6, 1, 5, 4, 3, 2, 1, 35, 8, 2, 1, 6, 1, 5, 4, 3, 2, 1, 27, 6, 1, 5, 4, 3, 2, 1, 21, 5, 4, 3, 2, 1, 16, 4, 3, 2, 1, 12, 3, 2, 1, 9, 2, 1, 7, 1, 6, 5, 4, 3, 2, 1, 154, 35, 8, 2, 1, 6, 1, 5, 4, 3, 2, 1, 27, 6, 1, 5, 4, 3, 2, 1, 21, 5, 4, 3, 2, 1, 16, 4, 3, 2, 1, 12, 3, 2, 1, 9, 2, 1, 7, 1, 6, 5, 4, 3, 2, 1, 119, 27, 6, 1, 5, 4, 3, 2, 1, 21, 5, 4, 3, 2, 1, 16, 4, 3, 2, 1, 12, 3, 2, 1, 9, 2, 1, 7, 1, 6, 5, 4, 3, 2, 1, 92, 21, 5, 4, 3, 2, 1, 16, 4, 3, 2, 1, 12, 3, 2, 1, 9, 2, 1, 7, 1, 6, 5, 4, 3, 2, 1, 71, 16, 4, 3, 2, 1, 12, 3, 2, 1, 9, 2, 1, 7, 1, 6, 5, 4, 3, 2, 1, 55, 12, 3, 2, 1, 9, 2, 1, 7, 1, 6, 5, 4, 3, 2, 1, 43, 9, 2, 1, 7, 1, 6, 5, 4, 3, 2, 1, 34, 7, 1, 6, 5, 4, 3, 2, 1, 27, 6, 5, 4, 3, 2, 1, 21, 5, 4, 3, 2, 1, 16, 4, 3, 2, 1, 12, 3, 2, 1, 9, 2, 1, 7, 1, 6, 5, 4, 3, 2, 1]

# point
# 128
glueP1 = (32314166027298155179147420812445594821868890135203684362438075953322320890420053859126669034545056421662113839282046917*i + 204875198158317628090707080901138329537258803315834400847668999577961785852468692151526403357244126238616223666853390478, 376023888332704841032903773824839287032887224444722404874075459743346955918158330515805689684404183952325634672952270779*i + 202559214124421541396086239080159715748750866436143931045820747206496669019616641730759902859010858596871461191164123789)
glueQ1 = (375319323555590268577784155059616891213238044726690980845910619479661168755097870626998963674663616639580798899045337055*i + 24453241683919606556036468683589818707708210318722282836208763788494526709185221886539081983664576835266852514898537817, 187374206753624164154051927203969664850966909537800131673919140880577477764958310467180619568534521632471163016432163218*i + 325001477787759978046708904112398664618254827856621586206813804209730554717526361992698555370081626158901628179594773811)
glueP2 = (147428986658806566494347131078535196195062339071995349685640646854492920297365650962593345296451741824941968934563069203*i + 244074513322811037673961472129964693678446305985139756160122433956136435206810899898446814534882249441840467959920439560, 96026691274428108614990803455299008407097038864923014270146055176263967157485001309704572418084382911291135253369644242*i + 147736276371774519856169444571291219213427509457849323885623583447163897154443884434460704329103789943213922177871393720)
glueQ2 = (265574627140586221098656447589327011688831214844349668614055365661102250908602308901859567526161651736989173974256525971*i + 123183320348321512904392840818275846199275830339105790731856831142830009228868281447305576666750128923936418813725812865, 96884051666487685542290447705872696440365046446507424506463128588755322802603797334433694129890755702879743319130119928*i + 260496637977633215424440621743698717884082233566610212422225118462893874169361658138844148846774611052690904366751839095)
# 192
# glueP1 = (12314147138074671613997623030622518295930193427444427283531256816170421905069127032996033193565044702140748877809627554884991026479412436281581442089133762685170812793873216449611*i + 12029820092072526391633700095497067653431029068206072587819564644097656018846525138548012675731567563570065947426512237435379041936839443735652883355521422780030986124835378427748, 6587796251178914220565103273976418920385629462289790129923028772550849517955217009763613036343026206164871617979935398864132277296064588132010825205243694050403598545971735386067*i + 4511561282068596903978321630429024515231915478218623569183426654784356280984026006721935484228483901640227729251790094057022300163716500134742858283570639465481658071511927091231)
# glueQ1 = (8696664974264538663898621390169187342467984351828559280150199334028221429009110024615605557507890154371213686352092352866219244006535349527088883222886290700017398768308402211267*i + 3523456296114674115518514241321335879543838297715629132429409289759664149662900352924903042998721137611302329679302128238245225005309174704736617124515731805414046912857856368771, 14207201836586605852611735665352620181305430190656159689621399887353553142689335774755077437518001841951050910595566085554456309752581346176838487128710416770207380859682446224445*i + 3081387897490343354452778802100693579617570277925393275821522619930990812951338692177136040775376435574226840056057436176760384326484095700102064143940451398136337810053849652427)
# glueP2 = (1385468960105113126982195510113582062280164097234362553731290747816021745058753842371698427819517847501279532496172215270649410982955640996635576063261575690120689713357277408515*i + 5219880197386071432486447689912716884242216072170750614913160644996455204823256631189173292242681938608681402980316806990112406539115337857097275798378701869532167138352716779443, 12385617968135468371850018462342340202389993752392530762774557659121360500457623040038334980048489641119520198702726318981429842548702111375351246689669098825371049753329317963756*i + 10221851823403445409658628546941216733466194658784835863909674260159813891059484112124026938651945828539014957141028940474611790675402497548556986892608320434639619952067636180929)
# glueQ2 = (9394478192003621387111437737440181097106546011038737986551449795806638175391984587090399772617972122256590962857263657873885177868277804266497164865047899513184840979758766414009*i + 411742843287092035673168271431660533553312010177654798800946974002876126846463380749591666321088099422726416230103966592317121846970155735279486481664702899253582539384257024851, 10097536748975650790146818760588560179673038467755355383760312859058672333427477787933600470642142444370630751422580207828673565552999455544800088163248064486429500480834977561426*i + 1700567247745582792316198457000551755555668336044361756333406947353885479043932811440293638597256229017608008434836510163163208809688935292436587757804992632205963845680252090838)
# 256
# glueP1 = (40115596448416578303493519700826683274395518302345283395927589052141847661228456764651633199118069490721439750477125416807273685619084505086527944743056063955855358643624140769030480612025251293767715219264248261174392956211232487077025*i + 814317225524357075575645219861446708304411390899392021978694551696954044192805085508586652516852159796741087378546655266693086732971887365500349594880558995093185140249640731589471213824544435710335210784706035715929466885639548360262, 8865476727019810015752142979885638158404929450799152728564151438035118763107137520140718216952699077470773409134652960283218379866864464149109762309344190262398628389704532679395894008735203931019017623070465171876809890730338158425945*i + 16981101117425965525812389987268420934807345394828619443130043289390088835206451951584319426336359581358814759066885535995271885798198989667638182028146326349172326924724201885976370960679759827729646961870876421515748593364276604928182)
# glueQ1 = (29227760402716659416960401891052688545958942634044996040094053565952878076317468855909316027896863518455556860152256982262757968180546896765817962207360923949031705482882753079210154757924214316710158338308347561283043028123703027990237*i + 38625977518400337577429572495152361992713569249972393872859070881336411115078877218431081855984814844203816072928800131189196751467984032669930588392268175098839787644009157165493055028490676878000509547780334446289357910861278541712692, 21551887188117867678913587232086707600788292748589101538226045950797054444000648659966395259508904207231902844568350209275576804771597084367411261673367756999782925457362538916021030842525489092546866668577636801224848850375411039509896*i + 34244161023001492426601338543840311295422551018213926263856921085970965784648846040814682565906574986515486151965555225787300453867000909304512060362464710401580115078887307028240195814338798796040083153935299182237437348287893534085878)
# glueP2 = (12664596035661357242515924462670251721967472885562025595274540827807134109002759235412506529716793902419563633980638160291016403502871148180430070307566539699370913070718584364656896923406473702105786414854398518021684775136106228358757*i + 13925456739716686977731622780492031536862772055198036596253213932526994682074978231759356815849119091042981743565449771588733130981722468678430073296001703254907331644640049073801436977126860097728636197359673596458939667429461163284180, 22672383337133004842002079707589290458524625870696293208491331069089514769770772703038856604759647356858225030857588248459315402402403876333701311441297376003586627534879118503943369827141610740577539000030497380991722505872334395234801*i + 1608670984636529003866527452287745645894684037800448869881330776690329282647418530840033952281314732765560004536316501378224891940350279780318507291430904489531891519595372740772524000931353256180812233146446102081777459497499325937120)
# glueQ2 = (30556458906020523175405650040528657812068978227560610712434706380275101966721934242946165442495070076109072122799714894501074896829337799974433958024876686657134847028327183496453765099092381026271668831216493261459096448392144595809940*i + 27652480224181332311977023077894931739586462678045723466790806958110837164461736272203203761586049749307552172715505858342063011994969345454833665826939932271200692747116808010561905659913523856301838107355908546758378525783274151318481, 33590005995657193685482899818520289296493062996220148106545514279564228187669378260420513892394059328270524141829202564341811782125632408879016423099069169810341109827874696804076486533477502393780093124595648965021518997903622192619274*i + 34653786805451160390564543752978739066135821394529319092807070966196834073971427754033000430478752640781218167444606776644125540446237560314273037681751643361900204285546998486591345721645347520953642841942596787759190531792351252606504)

# curve
# 128
E1_x2 = 227260257209177369488637076700303690049635096982597075097756151898043717171662905644421681413318330373113972921600182842*i+291405311519409659121821528686002959378816703097882682681946734184520204746829432708031679330525284294019466789402414572
E2_x2 = 108438894232722577592233910022456017351459189532963321553383775348813645373182493536057443463416923579722660391616224518*i+200475912862667342720963566204592004623656861960914506949076641345035361741560740165713587325566219276940901499546547367
# 192
# E1_x2 = 2087828019606679520884103991304315916221209418899990495129646007240835266893050222551698804185386695500129155224314691216182362641707617831908871054875774309997711236926952880870*i+4777777380432942139236201884669310380344800805588388833176148064574487497218861413549055437747494033520087497387702901035993463990971251205999353750636690425315655010006538257374
# E2_x2 = 13654290957789581172620003523244717044125999368739165612799151055883831020695278476147340810367112889630734332324888271075084092989947102074264468729164568436716640129629739743177*i+11718243093644005006630044873925167292868412015479592226369202956357700276548888651068058473173883516659460374945830055902724954864237728039714549791015932835368900059193079596739
# 256
# E1_x2 = 2507457291512024826402220623143870380862594660835544867702305650590377110751558236354285332891591693230719789795231502630500292972794932097411621706356715082211051085137189250645242053806855560053213270103347080115929822818588892079305*i+11155350524968085878107355247967693562307365088607859307937886805965656909050564153523330409856882270467615806741038183770951324425462166332004362724991164891913284950035303685769998739888116192480730805477424326361061554902027120581671
# E2_x2 = 5986449269543298686326082876393991402318013617642949199113968517600287481654971697758304277397299943907131924927429313745546863311835519979928757907992754458722403795672090356562668509372207451516670553016085884318237146612034855594136*i+13339213860511723519562349216328088240496506059069306511679415573086645943602837595777497977405548169995596329556249459342250357893577294549844198024608601189371746842744223236533645388798814497418508341277583871558928547328445225916535

E1 = EllipticCurve(Fp2d, [0, E1_x2, 0, 1, 0])
E2 = EllipticCurve(Fp2d, [0, E2_x2, 0, 1, 0])
glueP1 = E1(glueP1)
glueQ1 = E1(glueQ1)
glueP2 = E2(glueP2)
glueQ2 = E2(glueQ2)

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
chain, _ = ri.split_richelot_chain(glueP1, glueQ1, glueP2, glueQ2, b, strategy)
phi = chain
# print(phi)


# ========================= #
#       End Profiling       #
# ========================= #

pr.disable()
pr.dump_stats("festa_keygen.cProfile")
uf.print_info(f"(2,2)-Isogeny took: {(time.time() -  setup_time):.3f} seconds")
p = pstats.Stats("festa_keygen.cProfile")
p.strip_dirs().sort_stats("cumtime").print_stats(50)