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
% 25	vCaN_I1		Phosphatase rate of I1 by CaN
% 26 	PP1 		Concentration of active PP1
% 27 	I1P 		Concentration of phosphorylated PP inhibitor 1
% 28    gam_u       Rate of subunits that are bound to Ca2+/CaM among
% unphosphorylated subunits. Constant only through assumptions
% 29    gam_p       Rate of subunits that are bound to Ca2+/CaM among
% phosphorylated subunits. Constant only through assumptions
% 30    k10         Rate of dephosphorylation from phosphatase
% 31	mu 			Concentration of CaMKII-NMDA complexes

% Variables for PSD dynamics
% 32	Bp0 		Number of CaMKII-NR2B in state 0 (fully dephosphorylated)
% 33	Bp1 		Number of CaMKII-NR2B in state 1
% 34    Bp2 		Number of CaMKII-NR2B in state 2
% 35    Bp3 		Number of CaMKII-NR2B in state 3
% 36    Bp4 		Number of CaMKII-NR2B in state 4
% 37    Bp5 		Number of CaMKII-NR2B in state 5
% 38    Bp6 		Number of CaMKII-NR2B in state 6
% 39	Bp7 		Number of CaMKII-NR2B in state 7
% 40	Bp8 		Number of CaMKII-NR2B in state 8
% 41	Bp9 		Number of CaMKII-NR2B in state 9
% 42	Bp10 		Number of CaMKII-NR2B in state 10
% 43	Bp11		Number of CaMKII-NR2B in state 11
% 44	Bp12		Number of CaMKII-NR2B in state 12 
% 45	Bp13		Number of CaMKII-NR2B in state 13
% 46    S1          Sum(m_i*Bp_i,i=0..13)
% 47    S2          Sum(m_i^2*Bp_i,i=0..13)
% 48 	U 			Phosphorylation potential from CaMKII-NMDA complexes

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
savePlots = false;

t0=0; tfinal=500;
base_step = 0.5;

% VARIABLES AND PARAMETERS

%Variables
syms C(t)
syms B0(t) B1(t) B2(t) B3(t) B4(t) B5(t) B6(t) B7(t) B8(t) B9(t) B10(t) B11(t) B12(t) B13(t) chi(t) nu(t) S0(t)
syms Sp(t) Su(t) AMPA_bnd(t) Cb(t) 
syms vPKA_I1(t) vCaN_I1(t)
syms PP1(t) I1P(t) mu(t)
syms gam_u(t) gam_p(t) k10(t)
syms Ca(t)
syms S1(t) S2(t) Bpsd0(t) Bpsd1(t) Bpsd2(t) Bpsd3(t) Bpsd4(t)
syms Bpsd5(t) Bpsd6(t) Bpsd7(t) Bpsd8(t) Bpsd9(t) Bpsd10(t) Bpsd11(t) Bpsd12(t) Bpsd13(t)
syms Uvar(t)
%Params
syms tauCa CaBas Stot CaM K5 K9 L1 L2 L3 L4 k6 k7 k8 KM k12 k11 km11 I10 PP10 Kdcan ncan
syms kcan0_I1 kcan_I1 Kdpka npka kpka0_I1 kpka_I1
syms M kbc kbpc kbp kbpp
syms hgt rad_spn rad_psd NA rate_PSD

% QUANTITIES AND EQUATIONS

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
[p0,i0] = getP0(C0, paramVals(20:23), paramVals(24:27), paramVals(16:19));

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
[y0, yp0] = decic(implicitDAE, t0, y0est, y0fix, zeros(31,1), [], opt_ode);

[tauCa, CaBas, Stot, CaM, K5, K9, L1, L2, L3, L4, k6, k7, k8, KM, k12, k11, km11, I10, PP10, Kdcan, ncan, kcan0_I1, kcan_I1, Kdpka, npka, kpka0_I1, kpka_I1, M, kbc, kbpc, kbp, hgt, rad_spn, rad_psd, NA, rate_PSD] = paramVals;

Cas = y0(1);
Cs = y0(2);
B0s = y0(3);
B1s = y0(4);
B2s = y0(5);
B3s = y0(6);
B4s = y0(7);
B5s = y0(8);
B6s = y0(9);
B7s = y0(10);
B8s = y0(11);
B9s = y0(12);
B10s = y0(13);
B11s = y0(14);
B12s = y0(15);
B13s = y0(16);
chis = y0(17);
nus = y0(18);
S0s = y0(19);
Sps = y0(20);
Sus = y0(21);
AMPA_bnds = y0(22);
Cbs = y0(23);
vPKA_I1s = y0(24);
vCaN_I1s = y0(25);
PP1s = y0(26);
I1Ps = y0(27);
gam_us = y0(28);
gam_ps = y0(29);
k10s = y0(30);
mus = y0(31);

ts = 0;
dt = base_step;

while t < tmax
    % ODEs update
    tCa = Ca;
    tB1 = B1;
    tB2 = B2;
    tB3 = B3;
    tB4 = B4;
    tB5 = B5;
    tB6 = B6;
    tB7 = B7;
    tB8 = B8;
    tB9 = B9;
    tB10 = B10;
    tB11 = B11;
    tB12 = B12;
    tB13 = B13;
    tPP1 = PP1;
    tI1P = I1P;
    tmu = mu;

    Ca = Ca + dt*( -1/(tauCa)*(Ca - CatBas));
    B1 = tB1 + dt*( 6*k6*gam_u^2*tB0 - 4*k6*gam_u^2*tB1 - chi*gam_u*tB1 + nu*(2*(tB2+tB3+tB4)-tB1));
    B2 = tB2 + dt*( k6*gam_u^2*tB1 + chi*gam_u*tB1 + nu*(3*(tB5+tB6+tB7+tB8)-2*tB2) - 3*k6*gam_u^2*tB2 - chi*gam_u*tB2);
    B3 = tB3 + dt*( 2*k6*gam_u^2*tB1 + nu*(3*(tB5+tB6+tB7+tB8)-2*tB3) - 3*k6*gam_u^2*tB3 - chi*gam_u*tB3);
    B4 = tB4 + dt*( k6*gam_u^2*tB1 + nu*(3*(tB5+tB6+tB7+tB8)-2*tB4) - 2*k6*gam_u^2*tB4 - 2*chi*gam_u*tB4);
    B5 = tB5 + dt*( k6*gam_u^2*(tB2+tB3) + chi*gam_u*tB2 + nu*(4*(tB9+tB10+tB11)-3*tB5) - 2*k6*gam_u^2*tB5 - chi*gam_u*tB5);
    B6 = tB6 + dt*( k6*gam_u^2*(tB2+tB3) + 2*chi*gam_u*tB4 + nu*(4*(tB9+tB10+tB11)-3*tB6) - k6*gam_u^2*tB6 - 2*chi*gam_u*tB6);
    B7 = tB7 + dt*( k6*gam_u^2*(tB2+2*tB4) + chi*gam_u*tB3 + nu*(4*(tB9+tB10+tB11)-3*tB7) - k6*gam_u^2*tB7 - 2*chi*gam_u*tB7);
    B8 = tB8 + dt*( k6*gam_u^2*tB3 + nu*(4*(tB9+tB10+tB11)-3*tB8) - 3*chi*gam_u*tB8 - (M-mu)*(3*kbc*gam_u + 3*kbpc*gam_p + 3*kbp*(1-gam_p))*tB8);
    B9 = tB9 + dt*( k6*gam_u^2*tB5 + chi*gam_u*(tB6+tB7) + nu*(5*tB12-4*tB9) - k6*gam_u^2*tB9 - chi*gam_u*tB9);
    B10 = tB10 + dt*( k6*gam_u^2*(tB5+tB6) + chi*gam_u*(tB7+tB8) + nu*(5*tB12-4*tB10) - 2*chi*gam_u*tB10);
    B11 = tB11 + dt*( k6*gam_u^2*tB7 + chi*gam_u*tB6 + nu*(5*tB12-4*tB11) - 2*chi*gam_u*tB11);
    B12 = tB12 + dt*( k6*gam_u^2*tB9 + chi*gam_u*(tB9+2*tB10+2*tB11) + nu*(6*tB13-5*tB12) - chi*gam_u*tB12);
    B13 = tB13 + dt*( chi*gam_u*tB12 - nu*6*tB13);
    PP1 = tPP1 + dt*( -k11*tI1P*tPP1 + km11*(PP10 - tPP1));
    I1P = tI1P + dt*( -k11*tI1P*tPP1 + km11*(PP10 - tPP1) + vPKA_I1*(I10-tI1P) - vCaN_I1*tI1P);
    mu = tmu + dt*( 6*(M-tmu)*(kbc*gam_u)*S0 + (M-tmu)*(kbpc*gam_p + kbp*(1-gam_p) - kbc*gam_u)*Sp);

    dCa = bool(Ca) * (Ca - tCa)/tCa;
    dB1 = bool(B1) * (B1 - tB1)/tB1;
    dB2 = bool(B2) * (B2 - tB2)/tB2;
    dB3 = bool(B3) * (B3 - tB3)/tB3;
    dB4 = bool(B4) * (B4 - tB4)/tB4;
    dB5 = bool(B5) * (B5 - tB5)/tB5;
    dB6 = bool(B6) * (B6 - tB6)/tB6;
    dB7 = bool(B7) * (B7 - tB7)/tB7;
    dB8 = bool(B8) * (B8 - tB8)/tB8;
    dB9 = bool(B9) * (B9 - tB9)/tB9;
    dB10 = bool(B10) * (B10 - tB10)/tB10;
    dB11 = bool(B11) * (B11 - tB11)/tB11;
    dB12 = bool(B12) * (B12 - tB12)/tB12;
    dB13 = bool(B13) * (B13 - tB13)/tB13;
    dPP1 = bool(PP1) * (PP1 - tPP1)/tPP1;
    dI1P = bool(I1P) * (I1P - tI1P)/tI1P;
    dmu = bool(mu) * (mu - tmu)/tmu;    

    maxrelvar = max([Ca, B1, B2, B3, B4, B5, B6, B7, B8, B9, B10, B11, B12, B13, PP1, I1P, mu])
    if maxrelvar > 3e-2
        dt = dt/2;
        B1 = tB1;
        B2 = tB2;
        B3 = tB3;
        B4 = tB4;
        B5 = tB5;
        B6 = tB6;
        B7 = tB7;
        B8 = tB8;
        B9 = tB9;
        B10 = tB10;
        B11 = tB11;
        B12 = tB12;
        B13 = tB13;
        PP1 = tPP1;
        I1P = tI1P;
        mu = tmu;
        continue
    else
        B1s = [B1s; B1];
        B2s = [B2s; B2];
        B3s = [B3s; B3];
        B4s = [B4s; B4];
        B5s = [B5s; B5];
        B6s = [B6s; B6];
        B7s = [B7s; B7];
        B8s = [B8s; B8];
        B9s = [B9s; B9];
        B10s = [B10s; B10];
        B11s = [B11s; B11];
        B12s = [B12s; B12];
        B13s = [B13s; B13];
        PP1s = [PP1s; PP1];
        I1Ps = [I1Ps; I1P];
        mus = [mus; mu];

        % Linear entities update
        S0s = [S0s; B1 + B2 + B3 + B4 + B5 + B6 + B7 + B8 + B9 + B10 + B11 + B12 + B13];
        B0s = [B0s; Stot - S0];

        k10s = [k10s; k12/(KM + Sp)];
        
        gam_us = [gam_us; C/(K5 + C)];
        gam_ps = [gam_ps; C/(K9 + C)];
        Sps = [Sps, B1 + 2*(B2 + B3 + B4) + 3*(B5 + B6 + B7 + B8) + 4*(B9 + B10 + B11) + 5*B12 + 6*B13];
        Sus = [Sus, 6*Stot - Sp];
        
        Cs = [Cs, CaM/(1 + L4/Ca + L3*L4/(Ca^2) + L2*L3*L4/(Ca^3) + L1*L2*L3*L4/(Ca^4))];
        Cbs = [Cbs, gam_u*Su + gam_p*Sp];

        chis = [chis, k7*gam_p + k8*(1 - gam_p)];
        nus = [nus, k10*PP1];
        
        AMPA_bnds = [AMPA_bnds; 2*(B2+B3+B4) + 6*(B5+B6+B7+B8) + 12*(B9+B10+B11) + 20*B12 + 30*B13];

        vPKA_I1s = [vPKA_I1s; (heaviside((C-Cb)^2 - 1e-6))*(kpka0_I1 + kpka_I1/(1 + (Kdpka/(C-Cb))^npka))];
        vCaN_I1s = [vCaN_I1s; (heaviside((C-Cb)^2 - 1e-6))*(kcan0_I1 + kcan_I1/(1 + (Kdcan/(C-Cb))^ncan))];


        % Step computation
        dt = base_step;
        t = t + dt;
        ts = [ts; t];
    end
end

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
            0.012; 0.1;
            33.3; CaM;
            0.1; 0.0001;
            0.1; 0.025; 0.32; 0.40;
            6; 6; 6;
            0.4; 6000;
            500; 0.1; 1; 0.2;
            0.053; 3; 0.1; 18;
            0.11; 8; 0.00359; 100;
            1000; 0; 0; 0;
            0.5; 0.2; 0.1; 6.02e23; 0.05
        ];
    else
        paramVals = [
            0.012; 0.1;
            33.3; CaM;
            0.1; 0.0001;
            20; 0.57; 100; 5;
            6; 6; 4.8;
            11; 1.72;
            500; 0.1; 1; 0.2;
            0.053; 3; 0.1; 18;
            0.11; 8; 0.00359; 100;
            1000; 0; 0; 0;
            0.5; 0.2; 0.1; 6.02e23; 0.05
        ];
    end
end


function euler_update(step_size)

end