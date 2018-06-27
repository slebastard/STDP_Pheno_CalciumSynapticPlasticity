% Naming conventions
% k --- reaction constant
% r --- reaction rate ------------------------------- r_{11} = k_{11}*[C_p]*[PP1_{free}]
% K --- dissociation constant ----------------------- K_5 = k_{-5} / k_{5}
% KM -- Michaelis constant for catalytic reaction --- KM_Dpu = (k_{-11} + k_{12})/k_{11}

% Vector variable
% ID 	Init 	Name 	Description
% 1		1		Ca 		Free calcium
% 2				C 		Free Ca2+/CaM
% 3				gam_u 	Rate of CaM-bound monomers among T-286 unphosphorilated CaMKII
% 4				gam_p	Rate of CaM-bound monomers among T-286 phosphorilated CaMKII
% 5				zet_u	Rate of T-306 phosphorylated monomers among T-286 unphosphorilated CaMKII
% 6 			zet_p	Rate of T-306 phosphorylated monomers among T-286 phosphorilated CaMKII
% 7				Su 		Total concentration of T-286 unphosphorylated CaMKII monomers	
% 8				Sp 		Total concentration of T-286 phosphorylated CaMKII monomers
% 9		2		B1 		Concentration of olygomer 1
% 10	3		B2 		Concentration of olygomer 2
% 11	4		B3 		Concentration of olygomer 3
% 12	5		B4 		Concentration of olygomer 4
% 13	6		B5 		Concentration of olygomer 5
% 14	7		B6 		Concentration of olygomer 6
% 15	8		B7 		Concentration of olygomer 7
% 16	9		B8 		Concentration of olygomer 8
% 17	10		B9 		Concentration of olygomer 9
% 18	11		B10 	Concentration of olygomer 10
% 19	12		B11		Concentration of olygomer 11
% 20	13		B12		Concentration of olygomer 12 
% 21	14		B13		Concentration of olygomer 13
% 22 	15		PP1 	Concentration of active PP1
% 23 	16		I1 		Concentration of phosphorylated PP inhibitor 1
% 24	17		mu 		Concentration of CaMKII-NMDA complexes
% 25 	18		U 		Phosphorylation potential from CaMKII-NMDA complexes
% 26	19		n0
% 27	20		n1


% PARAMS

% CALCIUM PARAMETERS
tauCa=0.012, CaBas=0.2;

%% CaMKII PATHWAY PARAMETERS
% total CaMKII conc
Stot=33.3;
% total CaM conc
CaM=10;
% Ca2+/CaM binding dissociation const, from G&B 2007. Assume K9=K5
K5=0.1, K9=0.0001;

% Intermediate Ca to CaM association rates, from G&B 2007
L1=0.1, L2=0.025, L3=0.32, L4=0.40;
% Autophosphorilation by unphosphorilated subunit - reaction constant from G&B 2007
k6=6;

% Autophosphorylation by phosphorylated subunit: CaM-bound, free, T305-phosphorylated - from G&B 2007
% Lets settle on these values
k7=6, k8=6, k19=6;

% T-305 autophosphorylation kinetics
k17=10, k18=0.0005;

% PP1-CaMKII to dephosCaMKII constants + Michaelis: Dpu -> Duu (G&B 2007), Dpp -> Dpu, Dup -> Du
KM=0.4;
k12=6000;

% I1-PP1 kinetics
k11=500, km11=0.1;
I10=1;
PP10=0.2;

Kdcan=0.053, ncan=3;
kcan0_I1=0.1, kcan_I1=18;
% Values for the next two lines need to be assessed
kcan0_endo=0.1, kcan_endo=18;
kPP10_pase=0.1, kPP1_pase=18;

Kdpka=0.11, npka=8;
kpka0_I1=0.00359, kpka_I1=100;
% Values for the next line need to be assessed
kpka0_phos=0.00359, kpka_phos=100;

% Values for the next two lines need to be assessed
kCK2_exo=0.0005;
kNMDA_bind=0.1;

% AMPA TRAFFICKING PARAMETERS
% conductance of one channel is assumed to be 20pS, and we assume a total of 1000 AMPAr
N=1000;
g0=0.0010, g1=0.0017, g2=0.0024;

% NMDA PARAMETERS
M=100;



%%%%%%% CONSTRAINTS & DEFINITIONS %%%%%%%

% CaMKII CONSERVATION
rr = sum(0,12)of(shift(B1,i'))
% p0 is whats left from total
B0 = Stot-rr

Sp = B1 + 2*(B2 + B3 + B4) + 3*(B5 + B6 + B7 + B8) + 4*(B9 + B10 + B11) + 5*B12 + 6*B13
Su = 6*b0i-Sp
AMPA_bnd = 2*(B2+B3+B4) + 6*(B5+B6+B7+B8) + 12*(B9+B10+B11) + 20*B12 + 30*B13

k10 == k12*PP1/(KM + (1+zet_p)*Sp + zet_u*Su);
C == CaM/(1 + L4/Ca + L3*L4/(Ca^2) + L2*L3*L4/(Ca^3) + L1*L2*L3*L4/(Ca^4));


%%%%%% KINETICS %%%%%%

daes = [
	(1-gam_u-zet_u)*(C - gam_u*Su - gam_p*Sp) - K5*gam_u == 0;
	(1-gam_p-zet_p)*(C - gam_u*Su - gam_p*Sp) - K9*gam_p == 0;
	k18*(1-gam_u-zet_u) - k10*zet_u*PP1 == 0;
	k17*(1-gam_p-zet_p) - k10*zet_p*PP1 == 0;
]

% function out = ratesDAE(y)
% 	idx_C = 2;
% 	idx_gamu = 3;
% 	idx_Su = 7;
% 	idx_Sp = 8;
% 	idx_PP1 = 22;
% 
% 	k10 = k12*y(idx_PP1)/(KM + (1+y(idx_gamu+3))*y(idx_Sp) + y(idx_gamu+2)*y(idx_Su));
% 	out = [
% 		(1-y(idx_gamu)-y(idx_gamu+2))*(y(idx_C) - y(idx_gamu)*y(idx_Su) - y(idx_gamu+1)*y(idx_Sp)) - K5*y(idx_gamu);
% 		(1-y(idx_gamu+1)-y(idx_gamu+3))*(y(idx_C) - y(idx_gamu)*y(idx_Su) - y(idx_gamu+1)*y(idx_Sp)) - K9*y(idx_gamu+1);
% 		k18*(1-y(idx_gamu)-y(idx_gamu+2)) - k10*y(idx_gamu+2)*y(idx_PP1);
% 		k17*(1-y(idx_gamu+1)-y(idx_gamu+3)) - k10*y(idx_gamu+3)*y(idx_PP1);
% 	];
% end

% PKA activity
vPKA_I1=kpka0_I1 + kpka_I1/(1 + (Kdpka/C)^npka)
vPKA_phos=kpka0_phos + kpka_phos/(1 + (Kdpka/C)^npka)

% CaN activity
vCaN_I1=kcan0_I1 + kcan_I1/(1 + (Kdcan/C)^ncan)
vCaN_endo=kcan0_endo + kcan_endo/(1 + (Kdcan/C)^ncan)

% PP1 activity
vPP1_pase=kPP10_pase + kPP1_pase/(1 + (Kdcan/C)^ncan)

% CaMKII activity
vCK2_exo=kCK2_exo*Sp


%%%%%% DYNAMICS %%%%%%

% Ca2+ equation
ode_Ca = diff(Ca)==(CaBas-Ca)/tauCa;

% Useful qties for following equations
chi = k7*gam_p + k8*(1-gam_p-zet_p) + k19*zet_p;
nu = k10*(1+zet_p);

% CaMKII concentration equations
ode_B1 = diff(B1) == 6*k6*gam_u^2*B0 - 4*k6*gam_u^2*B1 - chi*gam_u*B1 + nu*(2*(B2+B3+B4)-B1) - kNMDA_bind*(M-mu)*B1;
%
ode_B2 = diff(B2) == k6*gam_u^2*B1 + chi*gam_u*B1 + nu*(3*(B5+B6+B7+B8)-2*B2) - 3*k6*gam_u^2*B2 - chi*gam_u*B2 - 2*kNMDA_bind*(M-mu)*B2;
ode_B3 = diff(B3) == 2*k6*gam_u^2*B1 + nu*(3*(B5+B6+B7+B8)-2*B3) - 3*k6*gam_u^2*B3 - chi*gam_u*B3 - 2*kNMDA_bind*(M-mu)*B3;
ode_B4 = diff(B4) == k6*gam_u^2*B1 + nu*(3*(B5+B6+B7+B8)-2*B4) - 2*k6*gam_u^2*B4 - 2*chi*gam_u*B4 - 2*kNMDA_bind*(M-mu)*B4;
%
ode_B5 = diff(B5) == k6*gam_u^2*(B2+B3) + chi*gam_u*B2 + nu*(4*(B9+B10+B11)-3*B5) - 2*k6*gam_u^2*B5 - chi*gam_u*B5 - 3*kNMDA_bind*(M-mu)*B5;
ode_B6 = diff(B6) == k6*gam_u^2*(B2+B3) + 2*chi*gam_u*B4 + nu*(4*(B9+B10+B11)-3*B6) - k6*gam_u^2*B6 - 2*chi*gam_u*B6 - 3*kNMDA_bind*(M-mu)*B6;
ode_B7 = diff(B7) == k6*gam_u^2*(B2+2*B4) + chi*gam_u*B3 + nu*(4*(B9+B10+B11)-3*B7) - k6*gam_u^2*B7 - 2*chi*gam_u*B7 - 3*kNMDA_bind*(M-mu)*B7;
ode_B8 = diff(B8) == k6*gam_u^2*B3 + nu*(4*(B9+B10+B11)-3*B8) - 3*chi*gam_u*B8 - 3*kNMDA_bind*(M-mu)*B8;
%
ode_B9 = diff(B9) == k6*gam_u^2*B5 + chi*gam_u*(B6+B7) + nu*(5*B12-4*B9) - k6*gam_u^2*B9 - chi*gam_u*B9 - 4*kNMDA_bind*(M-mu)*B9;
ode_B10 = diff(B10) == k6*gam_u^2*(B5+B6) + chi*gam_u*(B7+B8) + nu*(5*B12-4*B10) - 2*chi*gam_u*B10 - 4*kNMDA_bind*(M-mu)*B10;
ode_B11 = diff(B11) == k6*gam_u^2*B7 + chi*gam_u*B6 + nu*(5*B12-4*B11) - 2*chi*gam_u*B11 - 4*kNMDA_bind*(M-mu)*B11;
%
ode_B12 = diff(B12) == k6*gam_u^2*B9 + chi*gam_u*(B9+2*B10+2*B11) + nu*(6*B13-5*B12) - chi*gam_u*B12 - 5*kNMDA_bind*(M-mu)*B11;
%
ode_B13 = diff(B13) == chi*gam_u*B12 - nu*6*B13 - 6*kNMDA_bind*(M-mu)*B13;

% Phosphatase dynamics
ode_PP1 = diff(PP1) == -k11*I1P*PP1 + km11*(PP10 - PP1);
% ORIGINAL: I1P'= -k11*I1P*PP1 + km11*(PP10 - PP1) + vPKA*I10 - vCaN*I1P
ode_I1P = diff(I1P) == -k11*I1P*PP1 + km11*(PP10 - PP1) + vPKA_I1*(I10-I1P) - vCaN_I1*I1P;

% NMDA binding dynamics - Version 1
ode_mu = diff(mu) == kNMDA_bind*(M-mu)*Sp;
ode_U = diff(U) == kNMDA_bind*(M-mu)*AMPA_bnd;

% Membrane AMPA dynamics
%n0'=vCK2_exo*(N-n0-n1) - vCaN_endo*n0 - vPKA_phos*n0*U + vPP1_pase*n1*PP1
%n1'=vPKA_phos*n0*U - vPP1_pase*n1*PP1

%W=g0*n0 + g1*n1

%Variables
syms Ca(t), gam_u(t), gam_p(t), zet_u(t), zet_p(t), B1(t), B2(t), B3(t), B4(t), B5(t), B6(t), B7(t), B8(t), B9(t), B10(t), B11(t), B12(t), B13(t), PP1(t), I1P(t), mu(t), U(t)
%Params
syms tauCa, CaBas, Stot,CaM, K5, K9, L1, L2, L3, L4, k6, k7, k8, k19, k17, k18, KM, k12, k11, km11, I10, PP10, Kdcan, ncan
syms kcan0_I1, kcan_I1, kcan0_endo, kcan_endo, kPP10_pase, kPP1_pase, Kdpka, npka, kpka0_I1, kpka_I1, kpka0_phos kpka_phos
syms kCK2_exo, kNMDA_bind, N, g0, g1, g2, M

% Putting together the ODE/DAE system

vars = [
	Ca(t),
	gam_u(t), gam_p(t), zet_u(t), zet_p(t),
	B1(t), B2(t), B3(t), B4(t), B5(t), B6(t), B7(t), B8(t), B9(t), B10(t), B11(t), B12(t), B13(t),
	PP1(t), I1P(t),
	mu(t), U(t)
];

initCond =
[
	1,
	0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0,
	100, 0
];

initCondFixed = 
[
	1,
	0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0,
	0, 0
];

params = [
	tauCa, CaBas,
	Stot, CaM,
	K5, K9,
	L1, L2, L3, L4,
	k6, k7, k8, k19, k17, k18,
	KM, k12,
	k11, km11, I10, PP10,
	Kdcan, ncan, kcan0_I1, kcan_I1, kcan0_endo, kcan_endo,
	kPP10_pase, kPP1_pase,
	Kdpka, npka, kpka0_I1, kpka_I1, kpka0_phos, kpka_phos,
	kCK2_exo, kNMDA_bind,
	N, g0, g1, g2, M
];

paramVals = [
	0.012, 0.2,
	33.3, 10,
	0.1, 0.0001,
	0.1, 0.025, 0.32, 0.40,
	6, 6, 6, 6,
	10, 0.0005,
	0.4, 6000,
	500, 0.1, 1, 0.2,
	0.053, 3, 0.1, 18, 0.1, 18,
	0.1, 18,
	0.11, 8, 0.00359, 100, 0.00359, 100,
	0.0005, 0.1,
	1000, 0.0010, 0.0017, 0.0024, 100
];

odes = [
	ode_Ca,
	ode_B1,
	ode_B2,
	ode_B3,
	ode_B4,
	ode_B5,
	ode_B6,
	ode_B7,
	ode_B8,
	ode_B9,
	ode_B10,
	ode_B11,
	ode_B12,
	ode_B13,
	ode_PP1,
	ode_I1P,
	ode_mu,
	ode_U
];

eqs = [odes;daes];
F = daeFunction(eqs, vars);

f = @(t, y, yp)  F(t, y, yp, paramVals)
t0=0;
[y0,yp0] = decic(f, t0, initCond, initCondFixed, zeros(22,1), zeros(22,1), opt)

f(t0, y0, yp0)
ode15i(f, [t0, tfinal], y0, yp0, opt)