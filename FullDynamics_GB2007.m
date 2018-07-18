clear all;

% Naming conventions
% k --- reaction constant
% r --- reaction rate ------------------------------- r_{11} = k_{11}*[C_p]*[PP1_{free}]
% K --- dissociation constant ----------------------- K_5 = k_{-5} / k_{5}
% KM -- Michaelis constant for catalytic reaction --- KM_Dpu = (k_{-11} + k_{12})/k_{11}

% All time units in s
% All concentrations in microMolars = 1mumol/L

% Vector variable

% Variables for cytosolic CaMKII pathways
% ID 	Name 		Description
% 1		Ca 			Free calcium concentration
% 2		C 			Free Ca2+/CaM concentration
% 3		B0 			Concentration of olygomer 0 (fully dephosphorylated)
% 4		B1 			Concentration of olygomer 1
% 5     B2 			Concentration of olygomer 2
% 6     B3 			Concentration of olygomer 3
% 7     B4 			Concentration of olygomer 4
% 8     B5 			Concentration of olygomer 5
% 9     B6 			Concentration of olygomer 6
% 10	B7 			Concentration of olygomer 7
% 11	B8 			Concentration of olygomer 8
% 12	B9 			Concentration of olygomer 9
% 13	B10 		Concentration of olygomer 10
% 14	B11			Concentration of olygomer 11
% 15	B12			Concentration of olygomer 12 
% 16	B13			Concentration of olygomer 13
% 17	chi 		Intermediate expression for comp of B1'..B13'
% 18	nu 			Intermediate expression for comp of B1'..B13'
% 19	S0 			Sum(Bi,i=0..13)
% 20	Sp 			Total concentration of T-286 phosphorylated CaMKII monomers
% 21	Su 			Total concentration of T-286 unphosphorylated CaMKII monomers
% 22	AMPA_bnd	Concentration of CamKII monomer available for AMPA phosphorylation
% 23    Cb          Bound Ca2+/calmodulin concentration
% 24	vPKA_I1		Phosphorylation rate of I1 by PKA
% 25	vPKA_phos	Phosphorylation rate of AMPA by PKA
% 26	vCaN_I1		Phosphatase rate of I1 by CaN
% 27	vCaN_endo	Endocytosis rate of AMPA receptors by CaN
% 28	vPP1_pase	Phosphatase rate of AMPA receptors by PP1
% 29	vCK2_exo	Exocytosis rate of AMPA receptors by CaMKII-NMDA
% 30 	PP1 		Concentration of active PP1
% 31 	I1P 		Concentration of phosphorylated PP inhibitor 1
% 32    gam_u       Rate of subunits that are bound to Ca2+/CaM among
% unphosphorylated subunits
% 33    gam_p       Rate of subunits that are bound to Ca2+/CaM among
% phosphorylated subunits
% 34    k10         Rate of dephosphorylation from phosphatase
% 35	mu 			Concentration of CaMKII-NMDA complexes

% Variables for PSD dynamics
% 36	Bp0 		Number of CaMKII-NR2B in state 0 (fully dephosphorylated)
% 37	Bp1 		Number of CaMKII-NR2B in state 1
% 38    Bp2 		Number of CaMKII-NR2B in state 2
% 39    Bp3 		Number of CaMKII-NR2B in state 3
% 40    Bp4 		Number of CaMKII-NR2B in state 4
% 41    Bp5 		Number of CaMKII-NR2B in state 5
% 42    Bp6 		Number of CaMKII-NR2B in state 6
% 43	Bp7 		Number of CaMKII-NR2B in state 7
% 44	Bp8 		Number of CaMKII-NR2B in state 8
% 45	Bp9 		Number of CaMKII-NR2B in state 9
% 46	Bp10 		Number of CaMKII-NR2B in state 10
% 47	Bp11		Number of CaMKII-NR2B in state 11
% 48	Bp12		Number of CaMKII-NR2B in state 12 
% 49	Bp13		Number of CaMKII-NR2B in state 13
% 50    S1          Sum(m_i*B_i,i=0..13)
% 51    S2          Sum(m_i^2*B_i,i=0..13)
% 52 	U 			Phosphorylation potential from CaMKII-NMDA complexes
U; U(end) + (tSimu(stp+1) - tSimu(stp))*((paramVals(28)-y(end,35))*(6*paramVals(29)*y(end,32) - paramVals(30)*y(end,33) - paramVals(31)*(1-y(end,33)))*S1 + (M-mu)*(paramVals(30)*y(end,33) + paramVals(31)*(1-y(end,33)) - paramVals(29)*y(end,32))*S2
% Parameters
% ID    Name        Description
% 1     tauCa       Ca2+ exflux time constant
% 2     CaBas       Basal Ca2+ concentration
% 3     Stot        Total concentration of CaMKII hexamers
% 4     CaM         Total concentration of calmodulin
% 5     K5          Dissociation constant for CaM binding to unphosphorylated CaMKII
% 6     K9          Dissociation constant for CaM binding to phosphorylated CaMKII
% 7     L1          Dissociation constant for 1st Ca2+ binding to CaM
% 8     L2          Dissociation constant for 2nd Ca2+ binding to CaM
% 9     L3          Dissociation constant for 3rd Ca2+ binding to CaM
% 10    L4          Dissociation constant for 4rd Ca2+ binding to CaM
% 11    k6          Reaction constant for phosphorylation of subunit by CaM-bound, unphosphorylated neighbor
% 12    k7          Reaction constant for phosphorylation of subunit by CaM-bound, phosphorylated neighbor
% 13    k8          Reaction constant for phosphorylation of subunit by CaM-free, phosphorylated neighbor
% 14    KM          Michaelis-Menten constant for dephosphorylation by PP1
% 15    k12         Reaction constant for dissociation of PP1 from CaMKII
% 16    k11         Reaction constant for inhibition of PP1 by phosphorylated I1
% 17    km11        Reaction constant for reactivation of PP1, freeing one phosphorylated I1
% 18    I10         Total concentration of inhibitor 1
% 19    PP10        Total concentration of PP1
% 20    Kdcan       Calneurin half activity concentration
% 21    ncan        Calneurin Hill coefficient
% 22    kcan0_I1    Calneurin base activity
% 23    kcan_I1     Calneurin maximum Ca2+/CaM-dependent activity
% 24    Kdpka       PKA half activity concentration
% 25    npka        PKA Hill coefficient
% 26    kpka0_I1    PKA base activity
% 27    kpka_I1     PKA maximum Ca2+/CaM-dependent activity
% 28    M           Effective number of NR2B subunits (or scaffold binding sites) at PSD
% 29    kbc         Reaction constant for binding of CaM-bound, unphosphorylated CaMKII subunit to NR2B
% 30    kbpc        Reaction constant for binding of CaM-bound, phosphorylated CaMKII subunit to NR2B
% 31    kbp         Reaction constant for binding of CaM-free, phosphorylated CaMKII subunit to NR2B
% 32    hgt         Spine height
% 33    rad_spn     Spine radius
% 34    rad_psd     PSD radius
% 35    NA          Avogadro number
% 36    rate_PSD    Rate of CaMKII hexamers concentrated at the PSD

% INPUTS
paramSetName = 'Graupner';
in_CaInit = 0.1;
in_CaM = 0.1;
savePlots = true;

t0=0; tfinal=500;

% VARIABLES AND PARAMETERS

%Variables
syms C(t)
syms B0(t) B1(t) B2(t) B3(t) B4(t) B5(t) B6(t) B7(t) B8(t) B9(t) B10(t) B11(t) B12(t) B13(t) chi(t) nu(t) S0(t)
syms Sp(t) Su(t) AMPA_bnd(t) Cb(t) 
syms vPKA_I1(t) vPKA_phos(t) vCaN_I1(t) vCaN_endo(t) vPP1_pase(t) vCK2_exo(t)
syms PP1(t) I1P(t) mu(t)
syms gam_u(t) gam_p(t) k10(t)
syms Ca(t)
syms S1(t) S2(t) Bpsd0(t) Bpsd1(t) Bpsd2(t) Bpsd3(t) Bpsd4(t)
syms Bpsd5(t) Bpsd6(t) Bpsd7(t) Bpsd8(t) Bpsd9(t) Bpsd10(t) Bpsd11(t) Bpsd12(t) Bpsd13(t)
syms Uvar(t)
%Params
syms tauCa CaBas Stot CaM K5 K9 L1 L2 L3 L4 k6 k7 k8 k19 k17 k18 KM k12 k11 km11 I10 PP10 Kdcan ncan
syms kcan0_I1 kcan_I1 kcan0_endo kcan_endo kPP10_pase kPP1_pase Kdpka npka kpka0_I1 kpka_I1 kpka0_phos kpka_phos
syms kCK2_exo kNMDA_bind N g0 g1 g2 M
syms kbc kbpc kbp kbpp
syms hgt rad_spn rad_psd NA rate_PSD

% QUANTITIES AND EQUATIONS


daesCaMKII = [
    S0(t) == B1(t) + B2(t) + B3(t) + B4(t) + B5(t) + B6(t) + B7(t) + B8(t) + B9(t) + B10(t) + B11(t) + B12(t) + B13(t);
    B0(t) == Stot - S0(t);

    k10(t) == k12/(KM + Sp(t));
    
    gam_u(t) == C(t)/(K5 + C(t));
    gam_p(t) == C(t)/(K9 + C(t));
    Sp(t) == B1(t) + 2*(B2(t) + B3(t) + B4(t)) + 3*(B5(t) + B6(t) + B7(t) + B8(t)) + 4*(B9(t) + B10(t) + B11(t)) + 5*B12(t) + 6*B13(t);
    Su(t) == 6*Stot - Sp(t);
    
    C(t) == CaM/(1 + L4/Ca(t) + L3*L4/(Ca(t)^2) + L2*L3*L4/(Ca(t)^3) + L1*L2*L3*L4/(Ca(t)^4));
    Cb(t) == gam_u(t)*Su(t) + gam_p(t)*Sp(t);

    chi(t) == k7*gam_p(t) + k8*(1 - gam_p(t));
    nu(t) == k10*PP1(t);
    
    AMPA_bnd(t) == 2*(B2(t)+B3(t)+B4(t)) + 6*(B5(t)+B6(t)+B7(t)+B8(t)) + 12*(B9(t)+B10(t)+B11(t)) + 20*B12(t) + 30*B13(t);

    vPKA_I1(t) == (heaviside((C(t)-Cb(t))^2 - 1e-6))*(kpka0_I1 + kpka_I1/(1 + (Kdpka/(C(t)-Cb(t)))^npka));
    vCaN_I1(t) == (heaviside((C(t)-Cb(t))^2 - 1e-6))*(kcan0_I1 + kcan_I1/(1 + (Kdcan/(C(t)-Cb(t)))^ncan));
    % 14 DAEs
];
odesCaMKII = [
	diff(Ca(t), t) == -1/(tauCa)*(Ca(t) - CaBas);
    diff(B1(t), t) == 6*k6*gam_u^2*B0(t) - 4*k6*gam_u^2*B1(t) - chi(t)*gam_u*B1(t) + nu(t)*(2*(B2(t)+B3(t)+B4(t))-B1(t)) - (M-mu(t))*(5*kbc*gam_u(t) + kbpc*gam_p(t) + kbp*(1-gam_p(t)))*B1(t);
	diff(B2(t), t) == k6*gam_u^2*B1(t) + chi(t)*gam_u*B1(t) + nu(t)*(3*(B5(t)+B6(t)+B7(t)+B8(t))-2*B2(t)) - 3*k6*gam_u^2*B2(t) - chi(t)*gam_u*B2(t) - (M-mu(t))*(4*kbc*gam_u(t) + 2*kbpc*gam_p(t) + 2*kbp*(1-gam_p(t)))*B2(t);
	diff(B3(t), t) == 2*k6*gam_u^2*B1(t) + nu(t)*(3*(B5(t)+B6(t)+B7(t)+B8(t))-2*B3(t)) - 3*k6*gam_u^2*B3(t) - chi(t)*gam_u*B3(t) - (M-mu(t))*(4*kbc*gam_u(t) + 2*kbpc*gam_p(t) + 2*kbp*(1-gam_p(t)))*B3(t);
	diff(B4(t), t) == k6*gam_u^2*B1(t) + nu(t)*(3*(B5(t)+B6(t)+B7(t)+B8(t))-2*B4(t)) - 2*k6*gam_u^2*B4(t) - 2*chi(t)*gam_u*B4(t) - (M-mu(t))*(4*kbc*gam_u(t) + 2*kbpc*gam_p(t) + 2*kbp*(1-gam_p(t)))*B4(t);
	diff(B5(t), t) == k6*gam_u^2*(B2(t)+B3(t)) + chi(t)*gam_u*B2(t) + nu(t)*(4*(B9(t)+B10(t)+B11(t))-3*B5(t)) - 2*k6*gam_u^2*B5(t) - chi(t)*gam_u*B5(t) - (M-mu(t))*(3*kbc*gam_u(t) + 3*kbpc*gam_p(t) + 3*kbp*(1-gam_p(t)))*B5(t);
	diff(B6(t), t) == k6*gam_u^2*(B2(t)+B3(t)) + 2*chi(t)*gam_u*B4(t) + nu(t)*(4*(B9(t)+B10(t)+B11(t))-3*B6(t)) - k6*gam_u^2*B6(t) - 2*chi(t)*gam_u*B6(t) - (M-mu(t))*(3*kbc*gam_u(t) + 3*kbpc*gam_p(t) + 3*kbp*(1-gam_p(t)))*B6(t);
	diff(B7(t), t) == k6*gam_u^2*(B2(t)+2*B4(t)) + chi(t)*gam_u*B3(t) + nu(t)*(4*(B9(t)+B10(t)+B11(t))-3*B7(t)) - k6*gam_u^2*B7(t) - 2*chi(t)*gam_u*B7(t) - (M-mu(t))*(3*kbc*gam_u(t) + 3*kbpc*gam_p(t) + 3*kbp*(1-gam_p(t)))*B7(t);
	diff(B8(t), t) == k6*gam_u^2*B3(t) + nu(t)*(4*(B9(t)+B10(t)+B11(t))-3*B8(t)) - 3*chi(t)*gam_u*B8(t) - (M-mu(t))*(3*kbc*gam_u(t) + 3*kbpc*gam_p(t) + 3*kbp*(1-gam_p(t)))*B8(t);
	diff(B9(t), t) == k6*gam_u^2*B5(t) + chi(t)*gam_u*(B6(t)+B7(t)) + nu(t)*(5*B12(t)-4*B9(t)) - k6*gam_u^2*B9(t) - chi(t)*gam_u*B9(t) - (M-mu(t))*(2*kbc*gam_u(t) + 4*kbpc*gam_p(t) + 4*kbp*(1-gam_p(t)))*B9(t);
	diff(B10(t), t) == k6*gam_u^2*(B5(t)+B6(t)) + chi(t)*gam_u*(B7(t)+B8(t)) + nu(t)*(5*B12(t)-4*B10(t)) - 2*chi(t)*gam_u*B10(t) - (M-mu(t))*(2*kbc*gam_u(t) + 4*kbpc*gam_p(t) + 4*kbp*(1-gam_p(t)))*B10(t);
	diff(B11(t), t) == k6*gam_u^2*B7(t) + chi(t)*gam_u*B6(t) + nu(t)*(5*B12(t)-4*B11(t)) - 2*chi(t)*gam_u*B11(t) - (M-mu(t))*(2*kbc*gam_u(t) + 4*kbpc*gam_p(t) + 4*kbp*(1-gam_p(t)))*B11(t);
	diff(B12(t), t) == k6*gam_u^2*B9(t) + chi(t)*gam_u*(B9(t)+2*B10(t)+2*B11(t)) + nu(t)*(6*B13(t)-5*B12(t)) - chi(t)*gam_u*B12(t) - (M-mu(t))*(kbc*gam_u(t) + 5*kbpc*gam_p(t) + 5*kbp*(1-gam_p(t)))*B12(t);
	diff(B13(t), t) == chi(t)*gam_u*B12(t) - nu(t)*6*B13(t) - (M-mu(t))*(kbc*gam_u(t) + 6*kbpc*gam_p(t) + 6*kbp*(1-gam_p(t)))*B13(t);
	diff(PP1(t), t) == -k11*I1P(t)*PP1(t) + km11*(PP10 - PP1(t));
	diff(I1P(t), t) == -k11*I1P(t)*PP1(t) + km11*(PP10 - PP1(t)) + vPKA_I1(t)*(I10-I1P(t)) - vCaN_I1(t)*I1P(t);
    diff(mu(t), t) == 6*(M-mu(t))*(kbc*gam_u(t))*NR0(t) + (M-mu(t))*(kbpc*gam_p(t) + kbp*(1-gam_p(t)) - kbc*gam_u(t))*NR1(t);
       % 17 ODEs
];

eqsCaMKII = [odesCaMKII;daesCaMKII];

params = [
	tauCa; CaBas;
	Stot; CaM;
	K5; K9;
	L1; L2; L3; L4;
	k6; k7; k8;
	KM; k12;
	k11; km11; I10; PP10;
	Kdcan; ncan; kcan0_I1; kcan_I1;
	Kdpka; npka; kpka0_I1; kpka_I1;
	M;
    kbc; kbpc; kbp;
    hgt; rad_spn; rad_psd; NA; rate_PSD
];

paramVals = getParams(paramSetName, in_CaM);

C0 = getC(in_CaInit,paramVals(4),paramVals(7),paramVals(8),paramVals(9),paramVals(10));
[p0,i0] = getP0(C0, paramVals(23:26), paramVals(31:34), paramVals(19:22));

varsCaMKII = [
	Ca(t); C(t);
	B0(t); B1(t); B2(t); B3(t); B4(t); B5(t); B6(t); B7(t); B8(t);
    B9(t); B10(t); B11(t); B12(t); B13(t); chi(t); nu(t); S0(t);
	Sp(t); Su(t); AMPA_bnd(t); Cb(t);
	vPKA_I1(t); vCaN_I1(t);
	PP1(t); I1P(t);
    gam_u(t); gam_p(t); k10(t);
    mu(t)
];

y0est = [
	in_CaInit; C0;
	33.3; 0; 0; 0; 0; 0; 0; 0; 0;
    0; 0; 0; 0; 0; 0; 0; 0;
	0; 199.8; 0; 0;
	0; 0;
	p0; i0;
    0.5; 0.5; 10000;
    0
];

y0fix = [
	1; 0;
	1; 0; 0; 0; 0; 0; 0; 0; 0;
    0; 0; 0; 0; 0; 0; 0; 0;
	0; 0; 0; 0;
	0; 0;
	1; 1;
    0; 0; 0;
    0
];

[mass,F] = massMatrixForm(eqsCaMKII, varsCaMKII);
mass = odeFunction(mass, varsCaMKII);
F = odeFunction(F, varsCaMKII, params);
f = @(t, y) F(t, y, paramVals);
implicitDAE = @(t,y,yp) mass(t,y)*yp - f(t,y);

opt_ode = odeset('Mass', mass,...
'AbsTol',1e-5,...
'RelTol',1e-3);

% Solving system of eqs for CaMKII
[y0, yp0] = decic(implicitDAE, t0, y0est, y0fix, zeros(35,1), [], opt_ode);

opt_ode = odeset('Mass', mass,...
'AbsTol',1e-10, ...
'RelTol',1e-3);
    
[tSimu,yCaMKII] = ode15s(f, [t0, tfinal], y0', opt_ode);

% Compute NMDA variables

Bp0 = floor(max(0,1e-21* paramVals(52) * paramVals(51) * paramVals(48)*pi*paramVals(49)^2 * yCaMKII(:,3)));
Bp1 = floor(max(0,1e-21* paramVals(52) * paramVals(51) * paramVals(48)*pi*paramVals(49)^2 * yCaMKII(:,4)));
Bp2 = floor(max(0,1e-21* paramVals(52) * paramVals(51) * paramVals(48)*pi*paramVals(49)^2 * yCaMKII(:,5)));
Bp3 = floor(max(0,1e-21* paramVals(52) * paramVals(51) * paramVals(48)*pi*paramVals(49)^2 * yCaMKII(:,6)));
Bp4 = floor(max(0,1e-21* paramVals(52) * paramVals(51) * paramVals(48)*pi*paramVals(49)^2 * yCaMKII(:,7)));
Bp5 = floor(max(0,1e-21* paramVals(52) * paramVals(51) * paramVals(48)*pi*paramVals(49)^2 * yCaMKII(:,8)));
Bp6 = floor(max(0,1e-21* paramVals(52) * paramVals(51) * paramVals(48)*pi*paramVals(49)^2 * yCaMKII(:,9)));
Bp7 = floor(max(0,1e-21* paramVals(52) * paramVals(51) * paramVals(48)*pi*paramVals(49)^2 * yCaMKII(:,10)));
Bp8 = floor(max(0,1e-21* paramVals(52) * paramVals(51) * paramVals(48)*pi*paramVals(49)^2 * yCaMKII(:,11)));
Bp9 = floor(max(0,1e-21* paramVals(52) * paramVals(51) * paramVals(48)*pi*paramVals(49)^2 * yCaMKII(:,12)));
Bp10 = floor(max(0,1e-21* paramVals(52) * paramVals(51) * paramVals(48)*pi*paramVals(49)^2 * yCaMKII(:,13)));
Bp11 = floor(max(0,1e-21* paramVals(52) * paramVals(51) * paramVals(48)*pi*paramVals(49)^2 * yCaMKII(:,14)));
Bp12 = floor(max(0,1e-21* paramVals(52) * paramVals(51) * paramVals(48)*pi*paramVals(49)^2 * yCaMKII(:,15)));
Bp13 = floor(max(0,1e-21* paramVals(52) * paramVals(51) * paramVals(48)*pi*paramVals(49)^2 * yCaMKII(:,16)));
S1v = Bp1 + 2*(Bp2 + Bp3 + Bp4) + 3*(Bp5 + Bp6 + Bp7 + Bp8) + 4*(Bp9 + Bp10 + Bp11) + 5*Bp12 + 6*Bp13;
S2v = 2*(Bp2+Bp3+Bp4) + 6*(Bp5+Bp6+Bp7+Bp8) + 12*(Bp9+Bp10+Bp11) + 20*Bp12 + 30*Bp13;

y = [
    yCaMKII, ...
    Bp0, Bp1, Bp2, Bp3, Bp4, Bp5, Bp6, ...
    Bp7, Bp8, Bp9, Bp10, Bp11, Bp12, Bp13, ...
    S1v, S2v ...
];

varsNMDA = [
    Bpsd0(t); Bpsd1(t); Bpsd2(t); Bpsd3(t); Bpsd4(t); Bpsd5(t); Bpsd6(t);
    Bpsd7(t); Bpsd8(t); Bpsd9(t); Bpsd10(t); Bpsd11(t); Bpsd12(t); Bpsd13(t);
    S1(t); S2(t)
];

vars = [varsCaMKII; varsNMDA; Uvar(t)];

U = [0];
nstep = length(tSimu);
for stp=1:nstep-1
    U = [
        U; U(end) + (tSimu(stp+1) - tSimu(stp))*((M-mu)*(6*kb*gam_u - kb*gam_p - kb*(1-gam_p))*NR1 + (M-mu)*(kb*gam_p + kb*(1-gam_p) - kb*gam_u)*NR2
    ];
end

y = [y, U];

%%
plt_h=4; plt_l=4;
subt = sprintf('Response for impulse CaInit=%0.1f for total CaM=%0.2f', in_CaInit, in_CaM);
figName = sprintf('OutputGB_%s_CaInit%0.1f_CaM%0.2f', paramSetName, in_CaInit, in_CaM);
subplot = @(m,n,p) subtightplot (m, n, p, [0.05 0.035], [0.1 0.01], [0.1 0.01]);

for idx = 1:length(vars)
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
    plot(tSimu(:,1),y(:,idx), 'x')
    title(char(vars(idx))) 
end

ax = subtitle(subt);
axes(ax);
if savePlots
    saveas(fig, strcat(figName, sprintf('_%0d.png',1 + fix(idx/(plt_h*plt_l)))));
end

function K = getC(Ca, CaM, L1, L2, L3, L4)
    K = double(CaM/(1 + L4/Ca + L3*L4/(Ca^2) + L2*L3*L4/(Ca^3) + L1*L2*L3*L4/(Ca^4)));
end

function [p0,i0] = getP0(C, kin_CaN, kin_PKA, par_PP1)
    vCaN = kin_CaN(3) + kin_CaN(4)/(1 + (kin_CaN(1)/C)^kin_CaN(2));
    vPKA = kin_PKA(3) + kin_PKA(4)/(1 + (kin_PKA(1)/C)^kin_PKA(2));

    i0 = min(par_PP1(3)*vPKA/vCaN,par_PP1(3));
    p0 = par_PP1(4)/(1+i0*par_PP1(1)/par_PP1(2));
end

function paramVals = getParams(author,CaM)
    if strcmp(author,'Graupner')
        paramVals = [
            0.012; 0.05;
            33.3; CaM;
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
            0.012; 0.05;
            33.3; CaM;
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
    
    function y0 = getInitCond(eqs)
        
    end
end