% parameters shared with the ODE routine
%
global Mg1 Mg2 Mg De1 De2 taup dYbardZbarDifferentialMethod qa Qa qd Qd rhog0 T0 n n10 n20 tbar tbara tbarp tbardep tbard nout nouta noutp noutdep noutd noutcycle noutstep cycle zbar dzbar dtbar miu epsilon epsilonbar tc tt tstep k1 k2 k3 k4 q1 q2 Cs Cg rhob R dp L D A ta td tp tdep y1a y2a y1d y10 y2d Qin P0 Pa Pd Ta Td Pmax Pmin Tmax Tmin Tref Q0 kg MassAxial dtbara dtbarp dtbardep dtbard
%
% Select differential method for dYbardZbar
dYbardZbarDifferentialMethod = 777;
% dYbardZbarDifferentialMethod = 333; 1st order (2 points) upwind differe
% dYbardZbarDifferentialMethod = 444; 2nt order (3 points) upwind differe
% dYbardZbarDifferentialMethod = 777; 2nd order TVD Superbee Flux Limiter
%
% Select Mass Axial Dispersions
%
MassAxial = 881;
% Mass2=Axial = 881; Shut off mass axial dispersion
% MassAxial = 882; Turn on mass axial dispersion, need 2nd BC
%
% Select type of plots
%
choose = 2;
% choose = 1: 2D subplots of y1, y2, n1, n2, T, P
% choose = 2: Individual 2D plots of y1, y2, n1, n2, T, P
% choose = 3: Surface plots of y1, y2, n1, n2, T, P
% choose = 4: The contour plots
%
% Define number of nodes along t and along z
%
tstep = 1; % step time, sec
n = 15; % nodes along z
%
%
% Define cycles
cycle = 5; % Number of cycles
ta = 10; % adsorption time, sec
tp = 5; % pressurization time, sec
td = 10; % desorption time, sec
tdep = 5; % depressurization time, sec
tc = ta+tp+td+tdep; % cycle time, sec
tt = tc*cycle; % total time, sec
%
% Define dimensionless z and dz
%
zbarl = 0.0; % zbar lower limit, zbar = 0/L
zbaru = 1.0; % zbar upper limit, zbar = L/L
dzbar = (zbaru-zbarl)/(n-1); % zbar differential space
zbar = linspace(zbarl, zbaru, n);
%
% Define dimensionless t and dt
% Independent variable, tbar, for ODE integration
%
% tbar = t/tfeed: tbarl = tl/tfeed = 0/tfeed = 0, tbaru = tu/tfeed =
% tfeed/tfeed = 1;
tbarl = 0.0; % tbar lower limit, tbar = 0/tfeed
tbaru = 1.0; % tbar upper limit, tbar = tfeed
nouta = ta/tstep+1; % nodes of adsorption step
noutd = td/tstep+1; % nodes of desorption step
noutp = tp/tstep+1; % nodes of pressurization step
noutdep = tdep/tstep+1; % nodes of depressurization step
noutcycle = tc/tstep+1; % nodes of t along a cycle
nout = tt/tstep+1; % total nodes along t
dtbara = (tbaru-tbarl)/(nouta-1); % tbar differential space when adsorption
dtbard = (tbaru-tbarl)/(noutd-1); % tbar differential space when desorption
dtbarp = (tbaru-tbarl)/(noutp-1); % tbar differential space when pressurization
dtbardep = (tbaru-tbarl)/(noutdep-1); % tbar differential space when depressurization
dtbar = (tbaru-tbarl)/(nout-1); % tbar differential space for all the process
tbara = linspace(tbarl,dtbar*(nouta-1),nouta); % horizontal tbar when adsorption
tbard = linspace(tbarl,dtbar*(noutd-1),noutd); % horizontal tbar when desorption
tbarp = linspace(tbarl,dtbar*(noutp-1),noutp); % horizontal tbar when pressurization
tbardep = linspace(tbarl,dtbar*(noutdep-1),noutdep); % horizontal tbar when depressurization
%
% Packed bed dimensions
L = 30; %30.48; % bed length, cm % I chnaged from 30.38 to 30
D = 2.54; % bed diameter, cm
A = pi*(D/2)^2; % bed cross-sectional area, cm2
rhob = 1.574; % bulk density, g/cm3
%
% Adsorbent properties
%
dp = 0.100; % particle diameter, cm
epsilon = 0.51; % total column (helium) void fraction = ep x (1-e)+e
epsilonp = 0.3; % internal/intraparticle void fraction
epsilonbar = 0.3; % external/interparticle void fraction
taup = 1; % particle tortuosity factor
Cs = 0.28; % adsorbent (solid phase) heat capacity, cal/g/K
%
% Feed conditions at i = 1 (adsorption)
% at i = n (purge)
%
y1a = 0.21; % feed gas composition, 21% O2 + 79% N2
y2a = 0.79; % feed gas composition, 21% O2 + 79% N2
Ta = 873; % feed temperature, K
qa = 2.08; % feed mass flow rate, mmol/sec
Pa = 1.36*1013250; % pressure at i=n during adsorption, g/cm/sec2
y1d = 0; % feed gas composition, 0% O2 + 100% N2
y2d = 1; % feed gas composition, 0% O2 + 100% N2
Td = 873; % feed temperature, K
qd = 0.409; % feed mass flow rate, mmol/sec
Pd = 0.136*1013250; % pressure at i = 1 during desorption, g/cm/sec2
% Bed initial conditions
y10 = 0; % saturated with 0% IgG + 100% N2
y20 = 0;%1;
P0 = 10*101325; % pressure after pressurization step, g/cm/sec2
T0 = 873; % temperature after pressurization step, K
%
% Langmuir isotherms of O2 and N2 on adsorbent
% n1eqm0 = 0.1895*exp(925.0/ta)*2.534*exp(2865.0/ta)*(y10)*P0/(1+2.534*exp(2865/ta)*y10*P0);
n1eqm0 = 0;
n2eqm0 = 0;
n10 = n1eqm0; % O2 adsorbed, mmol O2/g
n20 = n2eqm0; % N2 adsorbed, mmol N2/g
%
% Gas peoperties
%
R = 8.314472*10^4; % gas constant, g.cm2/k/mmol/sec2
Cg = 0.00687; % gas phase heat capacity, cal/mmol/K
kg = 6.44e-5; % gas phase (Air at 1 atm 295K) thermal conductivity
% cal/cm/sec/K
miu = 18.2385e-5; % air dynamic viscosity at 25 C and 1 atm, g/cm/sec
Mg1 = 0.032; % O2 molecular weight, g/mmol
Mg2 = 0.028; % N2 molecular weight, g/mmol
q1 = 3.16; % isosteric heat of adsorption of O2 on adsorbent, cal/mmol
q2 = 5.60; % isosteric heat of adsorption of N2 on adsorbent, cal/mmol
%
% Mass flux and superficial velocity
Qa = qa/A; % mass flux, mmol/cm2/sec
Qd = qd/A;
% usfeed = Qa/rhog; % superficial velocity. cm/sec
%
Q0 = 0; % without gas flux, mmol/cm2/sec
rhog0 = P0/(R*T0); % gas phase density, mmol/cm3
us0 = Qin/rhog0; % superficial velocity, cm/sec
% Mass transfer kinetic
% Gas-solid mass transfer
E1 = 0.45*q1*4184; % J/mol
E2 = 0.45*q2*4184; % J/mol
De1 = 1.61*10^(-6)*exp(-E1/8.314/298)*10^4; % O2 effective diffusivity in adsorbent cm2/s
De2 = 1.61*10^(-6)*exp(-E2/8.314/298)*10^4; % N2 effective diffusivity in adsorbent cm2/s
k1 = 15*De1/(dp/2)^2; % O2 overall mass transfer coefficient, sec-1
k2 = 15*De2/(dp/2)^2; % N2 overall mass transfer coefficient, sec-1
%
% Additional parameters required for proper nondimensionalization of
% dependent variables
%
Pmax = Pa; % maximum pressure, g/cm/sec2
Pmin = 0; % minimum pressure, g/cm/sec2
Tref = T0;
Tmax = Ta; % maximum temperature, K
Tmin = 273; % minimum temperature, K