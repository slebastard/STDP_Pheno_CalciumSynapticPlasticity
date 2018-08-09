% function [t, varNames, y] = FullDynamics_NMDA(CaInputs, initVals, paramSetName, showPlots, savePlots, varParamName, varParamValue)

% Naming conventions
% k --- reaction constant
% r --- reaction rate ------------------------------- r_{11} = k_{11}*[C_p]*[PP1_{free}]
% K --- dissociation constant ----------------------- K_5 = k_{-5} / k_{5}
% KM -- Michaelis constant for catalytic reaction --- KM_Dpu = (k_{-11} + k_{12})/k_{11}

% All time units in s
% All concentrations in microMolars = 1mumol/L

% Vector variable
% ID 	Name 		Description
% 1		Ca 			Free calcium
% 2		C 			Free Ca2+/CaM
% 3		gam_u 		Rate of CaM-bound monomers among T-286 unphosphorilated CaMKII
% 4		gam_p		Rate of CaM-bound monomers among T-286 phosphorilated CaMKII
% 5		zet_u		Rate of T-306 phosphorylated monomers among T-286 unphosphorilated CaMKII
% 6 	zet_p		Rate of T-306 phosphorylated monomers among T-286 phosphorilated CaMKII
% 7		k10 		Dephosphorylation rate
% 8		B0 			Concentration of olygomer 0 (fully dephosphorylated)
% 9		B1 			Concentration of olygomer 1
% 10	B2 			Concentration of olygomer 2
% 11	B3 			Concentration of olygomer 3
% 12	B4 			Concentration of olygomer 4
% 13	B5 			Concentration of olygomer 5
% 14	B6 			Concentration of olygomer 6
% 15	B7 			Concentration of olygomer 7
% 16	B8 			Concentration of olygomer 8
% 17	B9 			Concentration of olygomer 9
% 18	B10 		Concentration of olygomer 10
% 19	B11			Concentration of olygomer 11
% 20	B12			Concentration of olygomer 12 
% 21	B13			Concentration of olygomer 13
% 22	chi 		Intermediate expression for comp of B1'..B13'
% 23	nu 			Intermediate expression for comp of B1'..B13'
% 24	rr 			Sum(B1..B13)
% 25	Sp 			Total concentration of T-286 phosphorylated CaMKII monomers
% 26	Su 			Total concentration of T-286 unphosphorylated CaMKII monomers
% 27	AMPA_bnd	Concentration of CamKII monomer available for AMPA phosphorylation
% 28	vPKA_I1		Phosphorylation rate of I1 by PKA
% 29	vPKA_phos	Phosphorylation rate of AMPA by PKA
% 30	vCaN_I1		Phosphatase rate of I1 by CaN
% 31	vCaN_endo	Endocytosis rate of AMPA receptors by CaN
% 32	vPP1_pase	Phosphatase rate of AMPA receptors by PP1
% 33	vCK2_exo	Exocytosis rate of AMPA receptors by CaMKII-NMDA
% 34 	PP1 		Concentration of active PP1
% 35 	I1 			Concentration of phosphorylated PP inhibitor 1
% 36	mu 			Concentration of CaMKII-NMDA complexes
% 37 	U 			Phosphorylation potential from CaMKII-NMDA complexes
% 26	n0 			Number of AMPA in phosphorylation state 0 (fully dephosphorylated)
% 27	n1 			Number of AMPA in phosphorylation state 1

CaInputs = 0.2;
initVals = [];
paramSetName = 'Graupner';
showPlots = true;
savePlots = true;

% if nargin~=5 && nargin~=7
%     error('Incorrect number of input arguments');
% end

% RECOMMENDED INPUTS
% paramSetName = 'Graupner';
% in_CaInit = 0.1;

t0=0; tf_p1=1; tf_p2=5; tfinal=50; 
step_p1=0.001; step_p2 = 0.003; step_p3 = 0.01;

% VARIABLES AND PARAMETERS

%Variables
global t_ Ca C B0 B1 B2 B3 B4 B5 B6 B7 B8 B9 B10 B11 B12 B13 chi nu rr Sp Su
global AMPA_bnd Cb vPKA_I1 vPKA_phos vCaN_I1 vCaN_endo vPP1_pase vCK2_exo
global PP1 I1P mu U NR2BbvC NR2BbvP Bpsd0 Bpsd1 Bpsd2 Bpsd3 Bpsd4 Bpsd5
global Bpsd6 Bpsd7 Bpsd8 Bpsd9 Bpsd10 Bpsd11 Bpsd12 Bpsd13 Uvar tauCa CaBas
global Stot CaM K5 K9 L1 L2 L3 L4 k6 k7 k8 k19 k17 k18 KM k12 k11 km11 I10
global PP10 Kdcan ncan kcan0_I1 kcan_I1 kcan0_endo kcan_endo kPP10_pase
global kPP1_pase Kdpka npka kpka0_I1 kpka_I1 kpka0_phos kpka_phos kCK2_exo
global kNMDA_bind N g0 g1 g2 M gam_u gam_p zet_u zet_p k10 kbc kbpc kbp
global kbpp hgt rad_spn rad_psd NA rate_PSD

syms Ca(t_) C(t_)
syms B0(t_) B1(t_) B2(t_) B3(t_) B4(t_) B5(t_) B6(t_) B7(t_) B8(t_) B9(t_) B10(t_) B11(t_) B12(t_) B13(t_) chi(t_) nu(t_) rr(t_)
syms Sp(t_) Su(t_) AMPA_bnd(t_) Cb(t_) 
syms vPKA_I1(t_) vPKA_phos(t_) vCaN_I1(t_) vCaN_endo(t_) vPP1_pase(t_) vCK2_exo(t_)
syms PP1(t_) I1P(t_) mu(t_) U(t_)

syms NR2BbvC(t_) NR2BbvP(t_) Bpsd0(t_) Bpsd1(t_) Bpsd2(t_) Bpsd3(t_) Bpsd4(t_)
syms Bpsd5(t_) Bpsd6(t_) Bpsd7(t_) Bpsd8(t_) Bpsd9(t_) Bpsd10(t_) Bpsd11(t_) Bpsd12(t_) Bpsd13(t_)
syms Uvar(t_)
%Params
syms tauCa CaBas Stot CaM K5 K9 L1 L2 L3 L4 k6 k7 k8 k19 k17 k18 KM k12 k11 km11 I10 PP10 Kdcan ncan
syms kcan0_I1 kcan_I1 kcan0_endo kcan_endo kPP10_pase kPP1_pase Kdpka npka kpka0_I1 kpka_I1 kpka0_phos kpka_phos
syms kCK2_exo kNMDA_bind N g0 g1 g2 M
syms gam_u gam_p zet_u zet_p k10

syms kbc kbpc kbp kbpp
syms hgt rad_spn rad_psd NA rate_PSD

% QUANTITIES AND EQUATIONS

% Membrane AMPA dynamics
% ode_n0 = diff(n0) == vCK2_exo*(N-n0-n1) - vCaN_endo*n0 - vPKA_phos*n0*U + vPP1_pase*n1*PP1;
% ode_n1 = iff(n1) == vPKA_phos*n0*U - vPP1_pase*n1*PP1;

% W = g0*n0 + g1*n1;

daesCaMKII = [
    rr(t_) == B1(t_) + B2(t_) + B3(t_) + B4(t_) + B5(t_) + B6(t_) + B7(t_) + B8(t_) + B9(t_) + B10(t_) + B11(t_) + B12(t_) + B13(t_);
    B0(t_) == Stot - rr(t_);

    Sp(t_) == B1(t_) + 2*(B2(t_) + B3(t_) + B4(t_)) + 3*(B5(t_) + B6(t_) + B7(t_) + B8(t_)) + 4*(B9(t_) + B10(t_) + B11(t_)) + 5*B12(t_) + 6*B13(t_);
    Su(t_) == 6*Stot - Sp(t_);

    C(t_) == CaM/(1 + L4/CaBas + L3*L4/(CaBas^2) + L2*L3*L4/(CaBas^3) + L1*L2*L3*L4/(CaBas^4));
    Cb(t_) == gam_u*Su(t_) + gam_p*Sp(t_);

    chi(t_) == k7*gam_p + k8*(1 - gam_p - zet_p) + k19*zet_p;
    nu(t_) == k10*PP1(t_)*(1+zet_p);
    
    AMPA_bnd(t_) == 2*(B2(t_)+B3(t_)+B4(t_)) + 6*(B5(t_)+B6(t_)+B7(t_)+B8(t_)) + 12*(B9(t_)+B10(t_)+B11(t_)) + 20*B12(t_) + 30*B13(t_);

    vPKA_I1(t_) == (heaviside((C(t_)-Cb(t_))^2 - 1e-6))*(kpka0_I1 + kpka_I1/(1 + (Kdpka/(C(t_)-Cb(t_)))^npka));
    vPKA_phos(t_) == (heaviside((C(t_)-Cb(t_))^2 - 1e-6))*(kpka0_phos + kpka_phos/(1 + (Kdpka/(C(t_)-Cb(t_)))^npka));

    vCaN_I1(t_) == (heaviside((C(t_)-Cb(t_))^2 - 1e-6))*(kcan0_I1 + kcan_I1/(1 + (Kdcan/(C(t_)-Cb(t_)))^ncan));
    vCaN_endo(t_) == (heaviside((C(t_)-Cb(t_))^2 - 1e-6))*(kcan0_endo + kcan_endo/(1 + (Kdcan/(C(t_)-Cb(t_)))^ncan));

    vPP1_pase(t_) == (heaviside((C(t_)-Cb(t_))^2 - 1e-6))*(kPP10_pase + kPP1_pase/(1 + (Kdcan/(C(t_)-Cb(t_)))^ncan));

    vCK2_exo(t_) == kCK2_exo*Sp(t_);
];

odesCaMKII = [
    diff(Ca(t_), t_) == -1/(tauCa)*(Ca(t_) - CaBas);
	diff(B1(t_), t_) == 6*k6*gam_u^2*B0(t_) - 4*k6*gam_u^2*B1(t_) - chi(t_)*gam_u*B1(t_) + nu(t_)*(2*(B2(t_)+B3(t_)+B4(t_))-B1(t_)) - kNMDA_bind*(M-mu(t_))*B1(t_);
	diff(B2(t_), t_) == k6*gam_u^2*B1(t_) + chi(t_)*gam_u*B1(t_) + nu(t_)*(3*(B5(t_)+B6(t_)+B7(t_)+B8(t_))-2*B2(t_)) - 3*k6*gam_u^2*B2(t_) - chi(t_)*gam_u*B2(t_) - 2*kNMDA_bind*(M-mu(t_))*B2(t_);
	diff(B3(t_), t_) == 2*k6*gam_u^2*B1(t_) + nu(t_)*(3*(B5(t_)+B6(t_)+B7(t_)+B8(t_))-2*B3(t_)) - 3*k6*gam_u^2*B3(t_) - chi(t_)*gam_u*B3(t_) - 2*kNMDA_bind*(M-mu(t_))*B3(t_);
	diff(B4(t_), t_) == k6*gam_u^2*B1(t_) + nu(t_)*(3*(B5(t_)+B6(t_)+B7(t_)+B8(t_))-2*B4(t_)) - 2*k6*gam_u^2*B4(t_) - 2*chi(t_)*gam_u*B4(t_) - 2*kNMDA_bind*(M-mu(t_))*B4(t_);
	diff(B5(t_), t_) == k6*gam_u^2*(B2(t_)+B3(t_)) + chi(t_)*gam_u*B2(t_) + nu(t_)*(4*(B9(t_)+B10(t_)+B11(t_))-3*B5(t_)) - 2*k6*gam_u^2*B5(t_) - chi(t_)*gam_u*B5(t_) - 3*kNMDA_bind*(M-mu(t_))*B5(t_);
	diff(B6(t_), t_) == k6*gam_u^2*(B2(t_)+B3(t_)) + 2*chi(t_)*gam_u*B4(t_) + nu(t_)*(4*(B9(t_)+B10(t_)+B11(t_))-3*B6(t_)) - k6*gam_u^2*B6(t_) - 2*chi(t_)*gam_u*B6(t_) - 3*kNMDA_bind*(M-mu(t_))*B6(t_);
	diff(B7(t_), t_) == k6*gam_u^2*(B2(t_)+2*B4(t_)) + chi(t_)*gam_u*B3(t_) + nu(t_)*(4*(B9(t_)+B10(t_)+B11(t_))-3*B7(t_)) - k6*gam_u^2*B7(t_) - 2*chi(t_)*gam_u*B7(t_) - 3*kNMDA_bind*(M-mu(t_))*B7(t_);
	diff(B8(t_), t_) == k6*gam_u^2*B3(t_) + nu(t_)*(4*(B9(t_)+B10(t_)+B11(t_))-3*B8(t_)) - 3*chi(t_)*gam_u*B8(t_) - 3*kNMDA_bind*(M-mu(t_))*B8(t_);
	diff(B9(t_), t_) == k6*gam_u^2*B5(t_) + chi(t_)*gam_u*(B6(t_)+B7(t_)) + nu(t_)*(5*B12(t_)-4*B9(t_)) - k6*gam_u^2*B9(t_) - chi(t_)*gam_u*B9(t_) - 4*kNMDA_bind*(M-mu(t_))*B9(t_);
	diff(B10(t_), t_) == k6*gam_u^2*(B5(t_)+B6(t_)) + chi(t_)*gam_u*(B7(t_)+B8(t_)) + nu(t_)*(5*B12(t_)-4*B10(t_)) - 2*chi(t_)*gam_u*B10(t_) - 4*kNMDA_bind*(M-mu(t_))*B10(t_);
	diff(B11(t_), t_) == k6*gam_u^2*B7(t_) + chi(t_)*gam_u*B6(t_) + nu(t_)*(5*B12(t_)-4*B11(t_)) - 2*chi(t_)*gam_u*B11(t_) - 4*kNMDA_bind*(M-mu(t_))*B11(t_);
	diff(B12(t_), t_) == k6*gam_u^2*B9(t_) + chi(t_)*gam_u*(B9(t_)+2*B10(t_)+2*B11(t_)) + nu(t_)*(6*B13(t_)-5*B12(t_)) - chi(t_)*gam_u*B12(t_) - 5*kNMDA_bind*(M-mu(t_))*B11(t_);
	diff(B13(t_), t_) == chi(t_)*gam_u*B12(t_) - nu(t_)*6*B13(t_) - 6*kNMDA_bind*(M-mu(t_))*B13(t_);
	diff(PP1(t_), t_) == -k11*I1P(t_)*PP1(t_) + km11*(PP10 - PP1(t_));
	diff(I1P(t_), t_) == -k11*I1P(t_)*PP1(t_) + km11*(PP10 - PP1(t_)) + vPKA_I1(t_)*(I10-I1P(t_)) - vCaN_I1(t_)*I1P(t_);
	diff(mu(t_), t_) == kNMDA_bind*(M-mu(t_))*Sp(t_);
       % TOTAL 18 ODEs
];

eqsCaMKII = [odesCaMKII;daesCaMKII];

params = [
	tauCa; CaBas;
	Stot; CaM;
	K5; K9;
	L1; L2; L3; L4;
	k6; k7; k8; k19; k17; k18;
	KM; k12;
	k11; km11; I10; PP10;
	Kdcan; ncan; kcan0_I1; kcan_I1; kcan0_endo; kcan_endo;
	kPP10_pase; kPP1_pase;
	Kdpka; npka; kpka0_I1; kpka_I1; kpka0_phos; kpka_phos;
	kCK2_exo; kNMDA_bind;
	N; g0; g1; g2; M;
    kbc; kbpc; kbp; kbpp;
    hgt; rad_spn; rad_psd; NA; rate_PSD;
    gam_u; gam_p; zet_u; zet_p; k10
];

paramNames = strings([length(params),1]);
for i=1:length(paramNames)
    paramNames(i) = char(params(i));
end

paramVals = getParams(paramSetName);

% if nargin==7
%     paramInd = (paramNames==varParamName);
%     paramVals(paramInd) = varParamValue;
% end

CaInit = CaInputs(1);

C0 = getC(CaInit,paramVals(4),paramVals(7),paramVals(8),paramVals(9),paramVals(10));
[p0,i0] = getP0(C0, paramVals(23:26), paramVals(31:34), paramVals(19:22));

opt_fsolve = optimoptions('fsolve','Display','off');
[gu,gp,zu,zp,lk10] = getRates(C0, 199.8, 0, paramVals(5), paramVals(6), paramVals(15), paramVals(16), paramVals(18), paramVals(17), p0, opt_fsolve);

paramVals = [
    paramVals;
    gu; gp; zu; zp; lk10
];

varsCaMKII = [
	C(t_);
	B0(t_); B1(t_); B2(t_); B3(t_); B4(t_); B5(t_); B6(t_); B7(t_); B8(t_);
    B9(t_); B10(t_); B11(t_); B12(t_); B13(t_); chi(t_); nu(t_); rr(t_);
	Sp(t_); Su(t_); AMPA_bnd(t_); Cb(t_);
	vPKA_I1(t_); vPKA_phos(t_); vCaN_I1(t_); vCaN_endo(t_); vPP1_pase(t_); vCK2_exo(t_);
	PP1(t_); I1P(t_);
    mu(t_);
    Ca(t_)
];

% Initiating variables from default values, then overwriting with values
% provided via arg(1) of input
y0est = [
	C0;
	33.3; 0; 0; 0; 0; 0; 0; 0; 0;
    0; 0; 0; 0; 0; 0; 0; 0;
	0; 199.8; 0; 0;
	0; 0; 0; 0; 0; 0;
	p0; i0;
    0;
    CaInit
];
if ~isempty(initVals)
    y0est = initVals;
    y0est(1)=C0;
    y0est(29)=p0;
    y0est(30)=i0;
    y0est(32) = CaInit;
end

y0fix = [
	0;
	0; 0; 0; 0; 1; 1; 1; 1; 1;
    1; 1; 1; 1; 1; 0; 0; 0;
	0; 0; 0; 0;
	0; 0; 0; 0; 0; 0;
	1; 1;
    0;
    1
];

[mass,F] = massMatrixForm(eqsCaMKII, varsCaMKII);
mass = odeFunction(mass, varsCaMKII);
F = odeFunction(F, varsCaMKII, params);
f = @(t, y) F(t, y, paramVals);

nsteps_p1 = (tf_p1 - t0)/step_p1;
nsteps_p2 = (tf_p2 - tf_p1)/step_p2;
nsteps_p3 = (tfinal - tf_p2)/step_p3;

opt_ode = odeset('Mass', mass,...
'AbsTol',1e-6,...
'RelTol',1e-3);

implicitDAE = @(t,y,yp) mass(t,y)*yp - f(t,y);
[y0, yp0] = decic(implicitDAE, t0, y0est, y0fix, zeros(32,1), [], opt_ode);

yi = y0';
t = []; yCaMKII = [];
ratesHist=[];
stop = 0;

opt_ode = odeset('Mass', mass,...
'AbsTol',1e-12, ...
'RelTol',1e-3);

nstep=1;
while ~stop && nstep~=nsteps_p1+1

    ratesSyst =  @(y) root(y, yi, paramVals);
    a = fsolve(ratesSyst,[paramVals(53),paramVals(54),paramVals(55),paramVals(56),paramVals(57)], opt_fsolve);

    paramVals(53:57) = a(1:5);

    f = @(t, y) F(t, y, paramVals);
    
    tstep = t0 + (nstep-1)*step_p1;
    [ti,yi] = ode15s(f, [tstep, tstep+step_p1], yi(end,:), opt_ode);
    t = [t;ti];
    yCaMKII = [yCaMKII;yi];  
    ratesHist = [ratesHist;repmat(paramVals(53:57).',length(ti),1)];
    
    stop = (tstep~=t0 && abs(yCaMKII(end,19)-yCaMKII(end-1,19))<1.5e-6*yCaMKII(end,19));
    nstep = nstep + 1;
end

nstep=1;
while ~stop && nstep~=nsteps_p2+1

    ratesSyst =  @(y) root(y, yi, paramVals);
    a = fsolve(ratesSyst,[paramVals(53),paramVals(54),paramVals(55),paramVals(56),paramVals(57)], opt_fsolve);

    paramVals(53:57) = a(1:5);

    f = @(t, y) F(t, y, paramVals);
    
    tstep = tf_p1 + (nstep-1)*step_p2;
    [ti,yi] = ode15s(f, [tstep, tstep+step_p2], yi(end,:), opt_ode);
    t = [t;ti];
    yCaMKII = [yCaMKII;yi];
    ratesHist = [ratesHist;repmat(paramVals(53:57).',length(ti),1)];
    
    stop = (tstep~=t0 && abs(yCaMKII(end,19)-yCaMKII(end-1,19))<1.5e-6*yCaMKII(end,19));
    nstep = nstep + 1;
end

nstep=1;
while ~stop && nstep~=nsteps_p3+1

    ratesSyst =  @(y) root(y, yi, paramVals);
    a = fsolve(ratesSyst,[paramVals(53),paramVals(54),paramVals(55),paramVals(56),paramVals(57)], opt_fsolve);

    paramVals(53:57) = a(1:5);

    f = @(t, y) F(t, y, paramVals);
    
    tstep = tf_p2 + (nstep-1)*step_p3;
    [ti,yi] = ode15s(f, [tstep, tstep+step_p3], yi(end,:), opt_ode);
    t = [t;ti];
    yCaMKII = [yCaMKII;yi];
    ratesHist = [ratesHist;repmat(paramVals(53:57).',length(ti),1)];
    
    stop = (tstep~=t0 && abs(yCaMKII(end,19)-yCaMKII(end-1,19))<1.5e-6*yCaMKII(end,19));
    nstep = nstep + 1;
end

%% NMDA binding

Bp0 = floor(max(0,1e-21* paramVals(52) * paramVals(51) * paramVals(48)*pi*paramVals(49)^2 * yCaMKII(:,2)));
Bp1 = floor(max(0,1e-21* paramVals(52) * paramVals(51) * paramVals(48)*pi*paramVals(49)^2 * yCaMKII(:,3)));
Bp2 = floor(max(0,1e-21* paramVals(52) * paramVals(51) * paramVals(48)*pi*paramVals(49)^2 * yCaMKII(:,4)));
Bp3 = floor(max(0,1e-21* paramVals(52) * paramVals(51) * paramVals(48)*pi*paramVals(49)^2 * yCaMKII(:,5)));
Bp4 = floor(max(0,1e-21* paramVals(52) * paramVals(51) * paramVals(48)*pi*paramVals(49)^2 * yCaMKII(:,6)));
Bp5 = floor(max(0,1e-21* paramVals(52) * paramVals(51) * paramVals(48)*pi*paramVals(49)^2 * yCaMKII(:,7)));
Bp6 = floor(max(0,1e-21* paramVals(52) * paramVals(51) * paramVals(48)*pi*paramVals(49)^2 * yCaMKII(:,8)));
Bp7 = floor(max(0,1e-21* paramVals(52) * paramVals(51) * paramVals(48)*pi*paramVals(49)^2 * yCaMKII(:,9)));
Bp8 = floor(max(0,1e-21* paramVals(52) * paramVals(51) * paramVals(48)*pi*paramVals(49)^2 * yCaMKII(:,10)));
Bp9 = floor(max(0,1e-21* paramVals(52) * paramVals(51) * paramVals(48)*pi*paramVals(49)^2 * yCaMKII(:,11)));
Bp10 = floor(max(0,1e-21* paramVals(52) * paramVals(51) * paramVals(48)*pi*paramVals(49)^2 * yCaMKII(:,12)));
Bp11 = floor(max(0,1e-21* paramVals(52) * paramVals(51) * paramVals(48)*pi*paramVals(49)^2 * yCaMKII(:,13)));
Bp12 = floor(max(0,1e-21* paramVals(52) * paramVals(51) * paramVals(48)*pi*paramVals(49)^2 * yCaMKII(:,14)));
Bp13 = floor(max(0,1e-21* paramVals(52) * paramVals(51) * paramVals(48)*pi*paramVals(49)^2 * yCaMKII(:,15)));
NR2C = 5*Bp1 + 8*(Bp2+Bp3+Bp4) + 9*(Bp5+Bp6+Bp7+Bp8) + 8*(Bp9+Bp10+Bp11) + 5*Bp12;
NR2P = 2*(Bp2+Bp3+Bp4) + 6*(Bp5+Bp6+Bp7+Bp8) + 12*(Bp9+Bp10+Bp11) + 20*Bp12 + 30*Bp13;

y = [
    yCaMKII, ...
    Bp0, Bp1, Bp2, Bp3, Bp4, Bp5, Bp6, ...
    Bp7, Bp8, Bp9, Bp10, Bp11, Bp12, Bp13, ...
    NR2C, NR2P ...
];

varsNMDA = [
    Bpsd0(t_); Bpsd1(t_); Bpsd2(t_); Bpsd3(t_); Bpsd4(t_); Bpsd5(t_); Bpsd6(t_);
    Bpsd7(t_); Bpsd8(t_); Bpsd9(t_); Bpsd10(t_); Bpsd11(t_); Bpsd12(t_); Bpsd13(t_);
    NR2BbvC(t_); NR2BbvP(t_)
];

vars = [varsCaMKII; varsNMDA; Uvar(t_); gam_u; gam_p; zet_u; zet_p; k10];

U = [0];
nstep = length(t);
for stp=1:nstep-1
    U = [
        U; U(end) + (t(stp+1) - t(stp))*((paramVals(43)-y(stp,31))*paramVals(44)*ratesHist(stp,1)*y(stp,20) + (paramVals(43)-y(stp,31))*(paramVals(45)*ratesHist(stp,2) + paramVals(46)*(1-ratesHist(stp,2))*y(stp,19)))
    ];
end

y = [y, U, ratesHist];
varNames = strings([size(y,2),1]);
for i=1:length(varNames)
    varNames(i) = char(vars(i));
end

%%
plt_h=4; plt_l=4;
subt = sprintf('Response at initial Ca=%0.1f, basal Ca=%0.1f for total CaM=%0.2f', CaInputs(1), paramVals(2), paramVals(4));
figName = sprintf('Output_%s_initCa%0.1f_steadyCa%0.1f_CaM%0.2f', paramSetName, CaInputs(1), paramVals(2), paramVals(4));
subplot = @(m,n,p) subtightplot (m, n, p, [0.05 0.035], [0.1 0.01], [0.1 0.01]);

if showPlots
    for idx = 1:size(y,2)
        if mod(idx,plt_h*plt_l)==1
            ax = subtitle(subt);
            axes(ax);
            if savePlots && idx > 1
                saveas(fig, strcat(figName, sprintf('_%0d.png',fix(idx/(plt_h*plt_l)))));
            end
            fig = figure(1 + fix(idx/(plt_h*plt_l)));
            set(gcf, 'Position', get(0, 'Screensize'));
        end
        h = subplot(plt_h,plt_l,1+mod(idx-1,plt_h*plt_l));
        p = get(h, 'Position');
        set(h,'pos',p+[0 -0.04 0 -0.04]);

        plot(t(:,1),y(:,idx), 'x')
        title(varNames(idx))

    end

    ax = subtitle(subt);
    axes(ax);
    if savePlots
        saveas(fig, strcat(figName, sprintf('_%0d.png',1 + fix(idx/(plt_h*plt_l)))));
    end
end
function K = getC(CaBas, CaM, L1, L2, L3, L4)
    K = double(CaM/(1 + L4/CaBas + L3*L4/(CaBas^2) + L2*L3*L4/(CaBas^3) + L1*L2*L3*L4/(CaBas^4)));
end

function [p0,i0] = getP0(C, kin_CaN, kin_PKA, par_PP1)
    vCaN = kin_CaN(3) + kin_CaN(4)/(1 + (kin_CaN(1)/C)^kin_CaN(2));
    vPKA = kin_PKA(3) + kin_PKA(4)/(1 + (kin_PKA(1)/C)^kin_PKA(2));

    i0 = min(par_PP1(3)*vPKA/vCaN,par_PP1(3));
    p0 = par_PP1(4)/(1+i0*par_PP1(1)/par_PP1(2));
end


function [gu, gp, zu, zp, lk10] = getRates(C, Su, Sp, K5, K9, k17, k18, k12, KM, PP1, opt)
    
    function r = root(x, C, Su, Sp, K5, K9, k17, k18, k12, KM, PP1)
        r(1) = (1-x(1)-x(3))*(C - x(1)*Su - x(2)*Sp) - K5*x(1);
        r(2) = (1-x(2)-x(4))*(C - x(1)*Su - x(2)*Sp) - K9*x(2);
        r(3) = k18*(1-x(1)-x(3)) - x(5)*x(3)*PP1;
        r(4) = k17*(1-x(2)-x(4)) - x(5)*x(4)*PP1;
        r(5) = x(5) - k12/(KM + (1+x(4))*Sp + x(3)*Su);
    end
    
    fun = @(x) root(x, C, Su, Sp, K5, K9, k17, k18, k12, KM, PP1);

    a = fsolve(fun,[0.5; 0.5; 0.5; 0.5; 1000], opt);
    gu = double(a(1));
    gp = double(a(2));
    zu = double(a(3));
    zp = double(a(4));
    lk10 = double(a(5));
end

function r = root(x, yi, paramVals)
        r(1) = (1-x(1)-x(3))*(yi(end,1) - x(1)*yi(end,20) - x(2)*yi(end,19)) - paramVals(5)*x(1);
        r(2) = (1-x(2)-x(4))*(yi(end,1) - x(1)*yi(end,20) - x(2)*yi(end,19)) - paramVals(6)*x(2);
        r(3) = paramVals(16)*(1-x(1)-x(3)) - x(5)*x(3)*yi(end,29);
        r(4) = paramVals(15)*(1-x(2)-x(4)) - x(5)*x(4)*yi(end,29);
        r(5) = x(5) - paramVals(18)/(paramVals(17) + (1+x(4))*yi(end,19) + x(3)*yi(end,20));
end

function paramVals = getParams(author)
    if strcmp(author,'Graupner')
        paramVals = [
            0.012; 0.1;
            33.3; 10;
            0.1; 0.0001;
            0.1; 0.025; 0.32; 0.40;
            6; 6; 6; 6; 10; 0.0005;
            0.4; 6000;
            500; 0.1; 1; 0.2;
            0.053; 3; 0.1; 18; 0.1; 18;
            0.1; 18;
            0.11; 8; 0.00359; 100; 0.00359; 100;
            0.0005; 0; %kNMDA_bind temporarilly set to 0
            1000; 0.0010; 0.0017; 0.0024; 100;
            0.01; 0.01; 0.01; 0.01; % NMDA binding is ON
            % 0; 0; 0; 0; % NMDA binding is OFF
            0.5; 0.2; 0.1; 6.02e23; 0.05
        ];
    else
        paramVals = [
            0.012; 0.1;
            33.3; 10;
            0.1; 0.0001;
            20; 0.57; 100; 5;
            6; 6; 4.8; 4.8; 10; 0.0005;
            11; 1.72;
            500; 0.1; 1; 0.2;
            0.053; 3; 0.1; 18; 0.1; 18;
            0.1; 18;
            0.11; 8; 0.00359; 100; 0.00359; 100;
            0.0005; 0; %kNMDA_bind temporarilly set to 0
            1000; 0.0010; 0.0017; 0.0024; 100;
            0.01; 0.01; 0.01; 0.01; % NMDA binding is ON
            % 0; 0; 0; 0; % NMDA binding is OFF
            % Params 48 to 52
            0.5; 0.2; 0.1; 6.02e23; 0.05
        ];
    end
end

% end