% Naming conventions
% k --- reaction constant
% r --- reaction rate ------------------------------- r_{11} = k_{11}*[C_p]*[PP1_{free}]
% K --- dissociation constant ----------------------- K_5 = k_{-5} / k_{5}
% KM -- Michaelis constant for catalytic reaction --- KM_Dpu = (k_{-11} + k_{12})/k_{11}

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

% VARIABLES AND PARAMETERS

%Variables
syms Ca(t) C(t)
syms gam_u(t) gam_p(t) zet_u(t) zet_p(t) k10(t)
syms B0(t) B1(t) B2(t) B3(t) B4(t) B5(t) B6(t) B7(t) B8(t) B9(t) B10(t) B11(t) B12(t) B13(t) chi(t) nu(t) rr(t)
syms Sp(t) Su(t) AMPA_bnd(t) 
syms vPKA_I1(t) vPKA_phos(t) vCaN_I1(t) vCaN_endo(t) vPP1_pase(t) vCK2_exo(t)
syms PP1(t) I1P(t) mu(t) U(t)
%Params
syms tauCa CaBas Stot CaM K5 K9 L1 L2 L3 L4 k6 k7 k8 k19 k17 k18 KM k12 k11 km11 I10 PP10 Kdcan ncan
syms kcan0_I1 kcan_I1 kcan0_endo kcan_endo kPP10_pase kPP1_pase Kdpka npka kpka0_I1 kpka_I1 kpka0_phos kpka_phos
syms kCK2_exo kNMDA_bind N g0 g1 g2 M


% QUANTITIES AND EQUATIONS


% DAEs
daes = [
	(1-gam_u(t)-zet_u(t))*(C(t) - gam_u(t)*Su(t) - gam_p(t)*Sp(t)) - K5*gam_u(t) == 0;
	(1-gam_p(t)-zet_p(t))*(C(t) - gam_u(t)*Su(t) - gam_p(t)*Sp(t)) - K9*gam_p(t) == 0;
	k18*(1-gam_u(t)-zet_u(t)) - k10(t)*zet_u(t)*PP1(t) == 0;
	k17*(1-gam_p(t)-zet_p(t)) - k10(t)*zet_p(t)*PP1(t) == 0;
    
    rr(t) == B1(t) + B2(t) + B3(t) + B4(t) + B5(t) + B6(t) + B7(t) + B8(t) + B9(t) + B10(t) + B11(t) + B12(t) + B13(t);
    B0(t) == Stot - rr(t);

    Sp(t) == B1(t) + 2*(B2(t) + B3(t) + B4(t)) + 3*(B5(t) + B6(t) + B7(t) + B8(t)) + 4*(B9(t) + B10(t) + B11(t)) + 5*B12(t) + 6*B13(t);
    Su(t) == 6*Stot - Sp(t);

    k10(t) == k12*PP1(t)/(KM + (1+zet_p(t))*Sp(t) + zet_u(t)*Su(t));
    C(t) == CaM/(1 + L4/Ca(t) + L3*L4/(Ca(t)^2) + L2*L3*L4/(Ca(t)^3) + L1*L2*L3*L4/(Ca(t)^4));

    chi(t) == k7*gam_p(t) + k8*(1 - gam_p(t) - zet_p(t)) + k19*zet_p(t);
    nu(t) == k10(t)*(1+zet_p(t));
    
    AMPA_bnd(t) == 2*(B2(t)+B3(t)+B4(t)) + 6*(B5(t)+B6(t)+B7(t)+B8(t)) + 12*(B9(t)+B10(t)+B11(t)) + 20*B12(t) + 30*B13(t);

    vPKA_I1(t) == kpka0_I1 + kpka_I1/(1 + (Kdpka/C(t))^npka);
    vPKA_phos(t) == kpka0_phos + kpka_phos/(1 + (Kdpka/C(t))^npka);

    vCaN_I1(t) == kcan0_I1 + kcan_I1/(1 + (Kdcan/C(t))^ncan);
    vCaN_endo(t) == kcan0_endo + kcan_endo/(1 + (Kdcan/C(t))^ncan);

    vPP1_pase(t) == kPP10_pase + kPP1_pase/(1 + (Kdcan/C(t))^ncan);

    vCK2_exo(t) == kCK2_exo*Sp(t);
                        % TOTAL 19 DAEs
];

% ODEs
% Ca2+ equation
ode_Ca = diff(Ca(t),t)==(CaBas-Ca(t))/tauCa;

% CaMKII concentration equations
ode_B1 = diff(B1(t), t) == 6*k6*gam_u(t)^2*B0(t) - 4*k6*gam_u(t)^2*B1(t) - chi(t)*gam_u(t)*B1(t) + nu(t)*(2*(B2(t)+B3(t)+B4(t))-B1(t)) - kNMDA_bind*(M-mu(t))*B1(t);
%
ode_B2 = diff(B2(t), t) == k6*gam_u(t)^2*B1(t) + chi(t)*gam_u(t)*B1(t) + nu*(3*(B5(t)+B6(t)+B7(t)+B8(t))-2*B2(t)) - 3*k6*gam_u(t)^2*B2(t) - chi(t)*gam_u(t)*B2(t) - 2*kNMDA_bind*(M-mu(t))*B2(t);
ode_B3 = diff(B3(t), t) == 2*k6*gam_u(t)^2*B1(t) + nu(t)*(3*(B5(t)+B6(t)+B7(t)+B8(t))-2*B3(t)) - 3*k6*gam_u(t)^2*B3(t) - chi(t)*gam_u(t)*B3(t) - 2*kNMDA_bind*(M-mu(t))*B3(t);
ode_B4 = diff(B4(t), t) == k6*gam_u(t)^2*B1(t) + nu(t)*(3*(B5(t)+B6(t)+B7(t)+B8(t))-2*B4(t)) - 2*k6*gam_u(t)^2*B4(t) - 2*chi(t)*gam_u(t)*B4(t) - 2*kNMDA_bind*(M-mu(t))*B4(t);
%
ode_B5 = diff(B5(t), t) == k6*gam_u(t)^2*(B2(t)+B3(t)) + chi(t)*gam_u(t)*B2(t) + nu(t)*(4*(B9(t)+B10(t)+B11(t))-3*B5(t)) - 2*k6*gam_u(t)^2*B5(t) - chi(t)*gam_u(t)*B5(t) - 3*kNMDA_bind*(M-mu(t))*B5(t);
ode_B6 = diff(B6(t), t) == k6*gam_u(t)^2*(B2(t)+B3(t)) + 2*chi(t)*gam_u(t)*B4(t) + nu(t)*(4*(B9(t)+B10(t)+B11(t))-3*B6(t)) - k6*gam_u(t)^2*B6(t) - 2*chi(t)*gam_u(t)*B6(t) - 3*kNMDA_bind*(M-mu(t))*B6(t);
ode_B7 = diff(B7(t), t) == k6*gam_u(t)^2*(B2(t)+2*B4(t)) + chi(t)*gam_u(t)*B3(t) + nu(t)*(4*(B9(t)+B10(t)+B11(t))-3*B7(t)) - k6*gam_u(t)^2*B7(t) - 2*chi(t)*gam_u(t)*B7(t) - 3*kNMDA_bind*(M-mu(t))*B7(t);
ode_B8 = diff(B8(t), t) == k6*gam_u(t)^2*B3(t) + nu*(4*(B9(t)+B10(t)+B11(t))-3*B8(t)) - 3*chi(t)*gam_u(t)*B8(t) - 3*kNMDA_bind*(M-mu(t))*B8(t);
%
ode_B9 = diff(B9(t), t) == k6*gam_u(t)^2*B5(t) + chi(t)*gam_u(t)*(B6(t)+B7(t)) + nu(t)*(5*B12(t)-4*B9(t)) - k6*gam_u(t)^2*B9(t) - chi(t)*gam_u(t)*B9(t) - 4*kNMDA_bind*(M-mu(t))*B9(t);
ode_B10 = diff(B10(t), t) == k6*gam_u(t)^2*(B5(t)+B6(t)) + chi(t)*gam_u(t)*(B7(t)+B8(t)) + nu(t)*(5*B12(t)-4*B10(t)) - 2*chi(t)*gam_u(t)*B10(t) - 4*kNMDA_bind*(M-mu(t))*B10(t);
ode_B11 = diff(B11(t), t) == k6*gam_u(t)^2*B7(t) + chi(t)*gam_u(t)*B6(t) + nu*(5*B12(t)-4*B11(t)) - 2*chi(t)*gam_u(t)*B11(t) - 4*kNMDA_bind*(M-mu(t))*B11(t);
%
ode_B12 = diff(B12(t), t) == k6*gam_u(t)^2*B9(t) + chi(t)*gam_u(t)*(B9(t)+2*B10(t)+2*B11(t)) + nu(t)*(6*B13(t)-5*B12(t)) - chi(t)*gam_u(t)*B12(t) - 5*kNMDA_bind*(M-mu(t))*B11(t);
%
ode_B13 = diff(B13(t), t) == chi(t)*gam_u(t)*B12(t) - nu(t)*6*B13(t) - 6*kNMDA_bind*(M-mu(t))*B13(t);

% Phosphatase dynamics
ode_PP1 = diff(PP1(t), t) == -k11*I1P(t)*PP1(t) + km11*(PP10 - PP1(t));
ode_I1P = diff(I1P(t), t) == -k11*I1P(t)*PP1(t) + km11*(PP10 - PP1(t)) + vPKA_I1(t)*(I10-I1P(t)) - vCaN_I1(t)*I1P(t);

% NMDA binding dynamics - Version 1
ode_mu = diff(mu(t), t) == kNMDA_bind*(M-mu(t))*Sp(t);
ode_U = diff(U(t), t) == kNMDA_bind*(M-mu(t))*AMPA_bnd(t);

% Membrane AMPA dynamics
% ode_n0 = diff(n0) == vCK2_exo*(N-n0-n1) - vCaN_endo*n0 - vPKA_phos*n0*U + vPP1_pase*n1*PP1;
% ode_n1 = iff(n1) == vPKA_phos*n0*U - vPP1_pase*n1*PP1;

% W = g0*n0 + g1*n1;

M = [
	1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
	0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
	0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
	0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
	0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
	0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
	0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
	0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
	0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
	0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
	0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
	0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
	0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
	0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
	0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0;
	0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0;
	0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0;
	0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1;
	0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
	0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
	0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
	0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
	0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
	0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
	0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
	0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
	0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
	0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
	0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
	0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
	0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
	0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
	0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
	0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
	0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
	0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
	0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
];

% Putting together the ODE/DAE system

vars = [
	Ca(t); C(t);
	gam_u(t); gam_p(t); zet_u(t); zet_p(t); k10(t);
	B0(t); B1(t); B2(t); B3(t); B4(t); B5(t); B6(t); B7(t); B8(t);
    B9(t); B10(t); B11(t); B12(t); B13(t); chi(t); nu(t); rr(t);
	Sp(t); Su(t); AMPA_bnd(t);
	vPKA_I1(t); vPKA_phos(t); vCaN_I1(t); vCaN_endo(t); vPP1_pase(t); vCK2_exo(t);
	PP1(t); I1P(t);
    mu(t); U(t)                                                          
];

initCond = [
	0.5; 0;
	0; 0; 0; 0; 0;
	0; 0; 0; 0; 0; 0; 0; 0; 0;
    0; 0; 0; 0; 0; 0; 0; 0;
	0; 0; 0;
	0; 0; 0; 0; 0; 0;
	0; 0;
    100; 0
];

initCondFixed = [
	1; 0;
	0; 0; 0; 0; 0;
	0; 0; 0; 0; 0; 0; 0; 0; 0;
    0; 0; 0; 0; 0; 0; 0; 0;
	0; 0; 0;
	0; 0; 0; 0; 0; 0;
	0; 0;
    0; 0
];

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
	N; g0; g1; g2; M
];

paramVals = [
	0.012; 0.2;
	33.3; 10;
	0.1; 0.0001;
	0.1; 0.025; 0.32; 0.40;
	6; 6; 6; 6;
	10; 0.0005;
	0.4; 6000;
	500; 0.1; 1; 0.2;
	0.053; 3; 0.1; 18; 0.1; 18;
	0.1; 18;
	0.11; 8; 0.00359; 100; 0.00359; 100;
	0.0005; 0.1;
	1000; 0.0010; 0.0017; 0.0024; 100
];

odes = [
	ode_Ca;
	ode_B1;
	ode_B2;
	ode_B3;
	ode_B4;
	ode_B5;
	ode_B6;
	ode_B7;
	ode_B8;
	ode_B9;
	ode_B10;
	ode_B11;
	ode_B12;
	ode_B13;
	ode_PP1;
	ode_I1P;
	ode_mu;
	ode_U
       % TOTAL 18 ODEs
];

%eqs = [odes;daes];
%F = daeFunction(eqs, vars, params);
%f = @(t, y, yp)  F(t, y, yp, paramVals);
%t0 = 0;
%tfinal = 300;
%opt = odeset('NonNegative', 1);
%[y0,yp0] = decic(f, t0, initCond, initCondFixed, zeros(37,1), zeros(37,1), opt);
%[t,y] = ode15i(f, [t0, tfinal], y0, yp0);

t0 = 0;
tfinal = 300;
opt = odeset('Mass', M, 'NonNegative', 1);
f = @(t, y, yp)  get_RHS(t, y, yp, paramVals);
[y0,yp0] = decic(f, t0, initCond, initCondFixed, zeros(37,1), zeros(37,1), opt);
[t,y] = ode15s(f, [t0, tfinal], y0, opt);

function ode_RHS = get_RHS(t,y)
	ode_RHS = [
		% ODEs
    	(CaBas-y(1))/tauCa;
    	6*k6*y(3)^2*y(8) - 4*k6*y(3)^2*y(9) - y(22)*y(3)*y(9) + y(23)*(2*(y(10)+y(11)+y(12))-y(9)) - kNMDA_bind*(M-y(36))*y(9);
		k6*y(3)^2*y(9) + y(22)*y(3)*y(9) + y(23)*(3*(y(13)+y(14)+y(15)+y(16))-2*y(10)) - 3*k6*y(3)^2*y(10) - y(22)*y(3)*y(10) - 2*kNMDA_bind*(M-y(36))*y(10);
		2*k6*y(3)^2*y(9) + y(23)*(3*(y(13)+y(14)+y(15)+y(16))-2*y(11)) - 3*k6*y(3)^2*y(11) - y(22)*y(3)*y(11) - 2*kNMDA_bind*(M-y(36))*y(11);
		k6*y(3)^2*y(9) + y(23)*(3*(y(13)+y(14)+y(15)+y(16))-2*y(12)) - 2*k6*y(3)^2*y(12) - 2*y(22)*y(3)*y(12) - 2*kNMDA_bind*(M-y(36))*y(12);
		k6*y(3)^2*(y(10)+y(11)) + y(22)*y(3)*y(10) + y(23)*(4*(y(17)+y(18)+y(19))-3*y(13)) - 2*k6*y(3)^2*y(13) - y(22)*y(3)*y(13) - 3*kNMDA_bind*(M-y(36))*y(13);
		k6*y(3)^2*(y(10)+y(11)) + 2*y(22)*y(3)*y(12) + y(23)*(4*(y(17)+y(18)+y(19))-3*y(14)) - k6*y(3)^2*y(14) - 2*y(22)*y(3)*y(14) - 3*kNMDA_bind*(M-y(36))*y(14);
		k6*y(3)^2*(y(10)+2*y(12)) + y(22)*y(3)*y(11) + y(23)*(4*(y(17)+y(18)+y(19))-3*y(15)) - k6*y(3)^2*y(15) - 2*y(22)*y(3)*y(15) - 3*kNMDA_bind*(M-y(36))*y(15);
		k6*y(3)^2*y(11) + y(23)*(4*(y(17)+y(18)+y(19))-3*y(16)) - 3*y(22)*y(3)*y(16) - 3*kNMDA_bind*(M-y(36))*y(16);
		k6*y(3)^2*y(13) + y(22)*y(3)*(y(14)+y(15)) + y(23)*(5*y(20)-4*y(17)) - k6*y(3)^2*y(17) - y(22)*y(3)*y(17) - 4*kNMDA_bind*(M-y(36))*y(17);
		k6*y(3)^2*(y(13)+y(14)) + y(22)*y(3)*(y(15)+y(16)) + y(23)*(5*y(20)-4*y(18)) - 2*y(22)*y(3)*y(18) - 4*kNMDA_bind*(M-y(36))*y(18);
		k6*y(3)^2*y(15) + y(22)*y(3)*y(14) + y(23)*(5*y(20)-4*y(19)) - 2*y(22)*y(3)*y(19) - 4*kNMDA_bind*(M-y(36))*y(19);
		k6*y(3)^2*y(17) + y(22)*y(3)*(y(17)+2*y(18)+2*y(19)) + y(23)*(6*y(21)-5*y(20)) - y(22)*y(3)*y(20) - 5*kNMDA_bind*(M-y(36))*y(19);
		y(22)*y(3)*y(20) - y(23)*6*y(21) - 6*kNMDA_bind*(M-y(36))*y(21);
		-k11*y(35)*y(34) + km11*(PP10 - y(34));
		-k11*y(35)*y(34) + km11*(PP10 - y(34)) + y(28)*(I10-y(35)) - y(30)*y(35);
		kNMDA_bind*(M-y(36))*y(25);
		kNMDA_bind*(M-y(36))*y(27);
		% DAEs
		(1-y(3)-y(5))*(y(2) - y(3)*y(26) - y(4)*y(25)) - K5*y(3);
		(1-y(4)-y(6))*(y(2) - y(3)*y(26) - y(4)*y(25)) - K9*y(4);
		k18*(1-y(3)-y(5)) - y(7)*y(5)*y(34);
		k17*(1-y(4)-y(6)) - y(7)*y(6)*y(34);
    	-y(24) + y(9) + y(10) + y(11) + y(12) + y(13) + y(14) + y(15) + y(16) + y(17) + y(18) + y(19) + y(20) + y(21);
    	-y(8) + Stot - y(24);
    	-y(25) + y(9) + 2*(y(10) + y(11) + y(12)) + 3*(y(13) + y(14) + y(15) + y(16)) + 4*(y(17) + y(18) + y(19)) + 5*y(20) + 6*y(21);
    	-y(26) + 6*Stot - y(25);
    	-y(7) + k12*y(34)/(KM + (1+y(6))*y(25) + y(5)*y(26));
    	-y(2) + CaM/(1 + L4/y(1) + L3*L4/(y(1)^2) + L2*L3*L4/(y(1)^3) + L1*L2*L3*L4/(y(1)^4));
    	-y(22) + k7*y(4) + k8*(1 - y(4) - y(6)) + k19*y(6);
    	-y(23) + y(7)*(1+y(6));
    	-y(27) + 2*(y(10)+y(11)+y(12)) + 6*(y(13)+y(14)+y(15)+y(16)) + 12*(y(17)+y(18)+y(19)) + 20*y(20) + 30*y(21);
    	-y(28) + kpka0_I1 + kpka_I1/(1 + (Kdpka/y(2))^npka);
    	-y(29) + kpka0_phos + kpka_phos/(1 + (Kdpka/y(2))^npka);
    	-y(30) + kcan0_I1 + kcan_I1/(1 + (Kdcan/y(2))^ncan);
    	-y(31) + kcan0_endo + kcan_endo/(1 + (Kdcan/y(2))^ncan);
    	-y(32) + kPP10_pase + kPP1_pase/(1 + (Kdcan/y(2))^ncan);
    	-y(33) + kCK2_exo*y(25)
	];
end