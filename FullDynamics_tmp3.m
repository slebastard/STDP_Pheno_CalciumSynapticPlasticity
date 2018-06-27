clear all;

NA=6.022140857e23;  

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

% Membrane AMPA dynamics
% ode_n0 = diff(n0) == vCK2_exo*(N-n0-n1) - vCaN_endo*n0 - vPKA_phos*n0*U + vPP1_pase*n1*PP1;
% ode_n1 = iff(n1) == vPKA_phos*n0*U - vPP1_pase*n1*PP1;

% W = g0*n0 + g1*n1;

mass = [
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
	0.012; 0.1;
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
	0.0005; 0; %kNMDA_bind temporarilly set to 0
	1000; 0.0010; 0.0017; 0.0024; 100
];

absTol = [
    1e-8; 1e-8;
    1e-8; 1e-8; 1e-8; 1e-8; 1e-8;
    1e-8; 1e-8; 1e-8; 1e-8; 1e-8; 1e-8; 1e-8; 1e-8; 1e-8;
    1e-8; 1e-8; 1e-8; 1e-8; 1e-8; 1e-8; 1e-8; 1e-8;
    1e-8; 1e-8; 1e-8;
    1e-8; 1e-8; 1e-8; 1e-8; 1e-8; 1e-8;
    1e-8; 1e-8;
    1e-2; 1e-2;
];

t0=0; tfinal=3000;
opt = odeset('Mass', mass, 'AbsTol',5e-4,'RelTol',1e-3);
f = @(t, y)  get_RHS(t, y, paramVals);

y0 = [
	0.5000; 4.2685;
	0.0000; 0.0000; 1.0000; 1.0000; 0.0000;
	33.3; 0; 0; 0; 0; 0; 0; 0; 0;
    0; 0; 0; 0; 0; 6.0000; 0; 0;
    0; 199.8; 0;
    100.0036; 100.0036; 18.1000; 18.1000; 18.1000; 0.0001;
    0.0000; 0;
	0; 0
];

[t,y] = ode15s(f, [t0, tfinal], y0, opt);

plt_h=4; plt_l=4;
for idx = 1:length(vars)
    if mod(idx,plt_h*plt_l)==1
        figure(1 + fix(idx/(plt_h*plt_l)))
        set(gcf, 'Position', get(0, 'Screensize'));
    end
    subplot(plt_h,plt_l,1+mod(idx-1,plt_h*plt_l))
    plot(t(:,1),y(:,idx), 'x')
    title(char(vars(idx)))
end

function ode_RHS = get_RHS(t, y, pars)
	ode_RHS = [
		% ODEs
    	(pars(2)-y(1))/pars(1);
    	6*pars(11)*y(3)^2*y(8) - 4*pars(11)*y(3)^2*y(9) - y(22)*y(3)*y(9) + y(23)*(2*(y(10)+y(11)+y(12))-y(9)) - pars(38)*(pars(43)-y(36))*y(9);
		pars(11)*y(3)^2*y(9) + y(22)*y(3)*y(9) + y(23)*(3*(y(13)+y(14)+y(15)+y(16))-2*y(10)) - 3*pars(11)*y(3)^2*y(10) - y(22)*y(3)*y(10) - 2*pars(38)*(pars(43)-y(36))*y(10);
		2*pars(11)*y(3)^2*y(9) + y(23)*(3*(y(13)+y(14)+y(15)+y(16))-2*y(11)) - 3*pars(11)*y(3)^2*y(11) - y(22)*y(3)*y(11) - 2*pars(38)*(pars(43)-y(36))*y(11);
		pars(11)*y(3)^2*y(9) + y(23)*(3*(y(13)+y(14)+y(15)+y(16))-2*y(12)) - 2*pars(11)*y(3)^2*y(12) - 2*y(22)*y(3)*y(12) - 2*pars(38)*(pars(43)-y(36))*y(12);
		pars(11)*y(3)^2*(y(10)+y(11)) + y(22)*y(3)*y(10) + y(23)*(4*(y(17)+y(18)+y(19))-3*y(13)) - 2*pars(11)*y(3)^2*y(13) - y(22)*y(3)*y(13) - 3*pars(38)*(pars(43)-y(36))*y(13);
		pars(11)*y(3)^2*(y(10)+y(11)) + 2*y(22)*y(3)*y(12) + y(23)*(4*(y(17)+y(18)+y(19))-3*y(14)) - pars(11)*y(3)^2*y(14) - 2*y(22)*y(3)*y(14) - 3*pars(38)*(pars(43)-y(36))*y(14);
		pars(11)*y(3)^2*(y(10)+2*y(12)) + y(22)*y(3)*y(11) + y(23)*(4*(y(17)+y(18)+y(19))-3*y(15)) - pars(11)*y(3)^2*y(15) - 2*y(22)*y(3)*y(15) - 3*pars(38)*(pars(43)-y(36))*y(15);
		pars(11)*y(3)^2*y(11) + y(23)*(4*(y(17)+y(18)+y(19))-3*y(16)) - 3*y(22)*y(3)*y(16) - 3*pars(38)*(pars(43)-y(36))*y(16);
		pars(11)*y(3)^2*y(13) + y(22)*y(3)*(y(14)+y(15)) + y(23)*(5*y(20)-4*y(17)) - pars(11)*y(3)^2*y(17) - y(22)*y(3)*y(17) - 4*pars(38)*(pars(43)-y(36))*y(17);
		pars(11)*y(3)^2*(y(13)+y(14)) + y(22)*y(3)*(y(15)+y(16)) + y(23)*(5*y(20)-4*y(18)) - 2*y(22)*y(3)*y(18) - 4*pars(38)*(pars(43)-y(36))*y(18);
		pars(11)*y(3)^2*y(15) + y(22)*y(3)*y(14) + y(23)*(5*y(20)-4*y(19)) - 2*y(22)*y(3)*y(19) - 4*pars(38)*(pars(43)-y(36))*y(19);
		pars(11)*y(3)^2*y(17) + y(22)*y(3)*(y(17)+2*y(18)+2*y(19)) + y(23)*(6*y(21)-5*y(20)) - y(22)*y(3)*y(20) - 5*pars(38)*(pars(43)-y(36))*y(19);
		y(22)*y(3)*y(20) - y(23)*6*y(21) - 6*pars(38)*(pars(43)-y(36))*y(21);
		-pars(19)*y(35)*y(34) + pars(20)*(pars(22) - y(34));
		-pars(19)*y(35)*y(34) + pars(20)*(pars(22) - y(34)) + y(28)*(pars(21)-y(35)) - y(30)*y(35);
		pars(38)*(pars(43)-y(36))*y(25);
		pars(38)*(pars(43)-y(36))*y(27);
		% DAEs
		(1-y(3)-y(5))*(y(2) - y(3)*y(26) - y(4)*y(25)) - pars(5)*y(3);
		(1-y(4)-y(6))*(y(2) - y(3)*y(26) - y(4)*y(25)) - pars(6)*y(4);
		pars(16)*(1-y(3)-y(5)) - y(7)*y(5)*y(34);
		pars(15)*(1-y(4)-y(6)) - y(7)*y(6)*y(34);
    	-y(24) + y(9) + y(10) + y(11) + y(12) + y(13) + y(14) + y(15) + y(16) + y(17) + y(18) + y(19) + y(20) + y(21);
    	-y(8) + pars(3) - y(24);
    	-y(25) + y(9) + 2*(y(10) + y(11) + y(12)) + 3*(y(13) + y(14) + y(15) + y(16)) + 4*(y(17) + y(18) + y(19)) + 5*y(20) + 6*y(21);
    	-y(26) + 6*pars(3) - y(25);
    	-y(7) + pars(18)*y(34)/(pars(17) + (1+y(6))*y(25) + y(5)*y(26));
    	-y(2) + pars(4)/(1 + pars(10)/y(1) + pars(9)*pars(10)/(y(1)^2) + pars(8)*pars(9)*pars(10)/(y(1)^3) + pars(7)*pars(8)*pars(9)*pars(10)/(y(1)^4));
    	-y(22) + pars(12)*y(4) + pars(13)*(1 - y(4) - y(6)) + pars(14)*y(6);
    	-y(23) + y(7)*(1+y(6));
    	-y(27) + 2*(y(10)+y(11)+y(12)) + 6*(y(13)+y(14)+y(15)+y(16)) + 12*(y(17)+y(18)+y(19)) + 20*y(20) + 30*y(21);
    	-y(28) + pars(33) + pars(34)/(1 + (pars(31)/y(2))^pars(32));
    	-y(29) + pars(35) + pars(36)/(1 + (pars(31)/y(2))^pars(32));
    	-y(30) + pars(25) + pars(26)/(1 + (pars(23)/y(2))^pars(24));
    	-y(31) + pars(27) + pars(28)/(1 + (pars(23)/y(2))^pars(24));
    	-y(32) + pars(29) + pars(30)/(1 + (pars(23)/y(2))^pars(24));
    	-y(33) + pars(37)*y(25)
	];
end