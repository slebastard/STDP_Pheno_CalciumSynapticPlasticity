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
in_CaM = 0.1;
savePlots = false;

t0=0; tfinal=500;
base_step = 0.5;
sigma = 0.5;
n_threads = 10;

% VARIABLES AND PARAMETERS

tauCa = 0.012;
CaBas = 0.1;
Stot = 33.3;
CaM = in_CaM;
K5 = 0.1;
K9 = 0.0001;
L1 = 0.1;
L2 = 0.025;
L3 = 0.32;
L4 = 0.40;
k6 = 6;
k7 = 6;
k8 = 6;
KM = 0.4;
k12 = 6000;
k11 = 500;
km11 = 0.1;
I10 = 1;
PP10 = 0.2;
Kdcan = 0.053;
ncan = 3;
kcan0_I1 = 0.1;
kcan_I1 = 18;
Kdpka = 0.11;
npka = 8;
kpka0_I1 = 0.00359;
kpka_I1 = 100;
M = 1000;
kbc = 0;
kbpc = 0;
kbp = 0;
hgt = 0.5;
rad_spn = 0.2;
rad_psd = 0.1;
NA = 6.02e23;
rate_PSD = 0.05;

% Sample uniformly and boostrap to get set of natural conditions

in_CaInit = 1;
t_CaInit = 10*rand(n_threads,1);

initConds = Stot.*uniform_simplex(n_samples, 14);

B0 = initConds(1);
B1 = initConds(2);
B2 = initConds(3);
B3 = initConds(4);
B4 = initConds(5);
B5 = initConds(6);
B6 = initConds(7);
B7 = initConds(8);
B8 = initConds(9);
B9 = initConds(10);
B10 = initConds(11);
B11 = initConds(12);
B12 = initConds(13);
B13 = initConds(14);

[y_init,t] = stimulate(in_CaInit, t_stim, B0, B1, B2, B3, B4, B5, B6, B7, B8, B9, B10, B11, B12, B13);

% Use those states as initial conditions and stimulate
B0 = y_init(3);
B1 = y_init(4);
B2 = y_init(5);
B3 = y_init(6);
B4 = y_init(7);
B5 = y_init(8);
B6 = y_init(9);
B7 = y_init(10);
B8 = y_init(11);
B9 = y_init(12);
B10 = y_init(13);
B11 = y_init(14);
B12 = y_init(15);
B13 = y_init(16);
t_stim = 10*rand(n_threads,1);

[y,t] = stimulate(in_CaInit, t_stim, B0, B1, B2, B3, B4, B5, B6, B7, B8, B9, B10, B11, B12, B13);

t_alpha = tstim + tauCa*log(Ca_stim./(in_CaInit-CaBas));
id_alpha = max(find(t > t_alpha), 2);

Sp = y(20);
Sp_stop = Sp(id_alpha,:);
Sp_ifty = Sp(end,:);

% Plotting

figure()
plot(Sp_stop, Sp_ifty);
title('Transfer curve from stopping criterion to stable system');
xlabel('Phosphorylated CaMKII at stop calcium concentration');
ylabel('Phosphorylated CaMKII in stable state');

function [y,t] = stimulate(in_CaInit, tstim, B0, B1, B2, B3, B4, B5, B6, B7, B8, B9, B10, B11, B12, B13)

    C0 = getC(in_CaInit, CaM, L1, L2, L3, L4);
    [p0,i0] = getP0(C0, [Kdcan,ncan,kcan0_I1,kcan_I1], [Kdpka,npka,kpka0_I1,kpka_I1], [k11,km11,I10,PP10]);

    Ca = in_CaInit; C = C0;
    chi = 0.0; nu = 0.0;
    S0 = B0 + B1 + B2 + B3 + B4 + B5 + B6 + B7 + B8 + B9 + B10 + B11 + B12 + B13;
    Sp = B1 + 2*(B2 + B3 + B4) + 3*(B5 + B6 + B7 + B8) + 4*(B9 + B10 + B11) + 5*B12 + 6*B13;
    Su = 6*Stot - Sp; AMPA_bnd = 2*(B2+B3+B4) + 6*(B5+B6+B7+B8) + 12*(B9+B10+B11) + 20*B12 + 30*B13; Cb = gam_u*Su + gam_p*Sp;
    vPKA_I1 = kpka0_I1 + kpka_I1./(1 + (Kdpka./(C-Cb)).^npka; vCaN_I1 = kcan0_I1 + kcan_I1./(1 + (Kdcan./(C-Cb)).^ncan;
    PP1 = p0; I1P = i0;
    gam_u = C0./(K5 + C0); gam_p = C0./(K9 + C0); k10 = k12./(KM + Sp);
    mu = 0.0;

    Cas = Ca; Cs = C;
    B0s = B0; B1s = B1; B2s = B2; B3s = B3; B4s = B4; B5s = B5; B6s = B6; B7s = B7; B8s = B8;
    B9s = B9; B10s = B10; B11s = B11; B12s = B12; B13s = B13; chis = chi; nus = nu; S0s = S0;
    Sps = Sp; Sus = Su; AMPA_bnds = AMPA_bnd; Cbs = Cb;
    vPKA_I1s = vPKA_I1; vCaN_I1s = vCaN_I1;
    PP1s = PP1; I1Ps = I1P;
    gam_us = gam_u; gam_ps = gam_p; k10s = k10; 
    mus = mu;

    t = 0;
    dt = base_step;

    while t < tfinal
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

        Ca = Ca + (t>tstim).*dt.*( -1./(tauCa).*(Ca - CaBas));
        B1 = tB1 + dt.*( 6.*k6.*gam_u.^2.*B0 - 4.*k6.*gam_u.^2.*tB1 - chi.*gam_u.*tB1 + nu.*(2.*(tB2+tB3+tB4)-tB1));
        B2 = tB2 + dt.*( k6.*gam_u.^2.*tB1 + chi.*gam_u.*tB1 + nu.*(3.*(tB5+tB6+tB7+tB8)-2.*tB2) - 3.*k6.*gam_u.^2.*tB2 - chi.*gam_u.*tB2);
        B3 = tB3 + dt.*( 2.*k6.*gam_u.^2.*tB1 + nu.*(3.*(tB5+tB6+tB7+tB8)-2.*tB3) - 3.*k6.*gam_u.^2.*tB3 - chi.*gam_u.*tB3);
        B4 = tB4 + dt.*( k6.*gam_u.^2.*tB1 + nu.*(3.*(tB5+tB6+tB7+tB8)-2.*tB4) - 2.*k6.*gam_u.^2.*tB4 - 2.*chi.*gam_u.*tB4);
        B5 = tB5 + dt.*( k6.*gam_u.^2.*(tB2+tB3) + chi.*gam_u.*tB2 + nu.*(4.*(tB9+tB10+tB11)-3.*tB5) - 2.*k6.*gam_u.^2.*tB5 - chi.*gam_u.*tB5);
        B6 = tB6 + dt.*( k6.*gam_u.^2.*(tB2+tB3) + 2.*chi.*gam_u.*tB4 + nu.*(4.*(tB9+tB10+tB11)-3.*tB6) - k6.*gam_u.^2.*tB6 - 2.*chi.*gam_u.*tB6);
        B7 = tB7 + dt.*( k6.*gam_u.^2.*(tB2+2.*tB4) + chi.*gam_u.*tB3 + nu.*(4.*(tB9+tB10+tB11)-3.*tB7) - k6.*gam_u.^2.*tB7 - 2.*chi.*gam_u.*tB7);
        B8 = tB8 + dt.*( k6.*gam_u.^2.*tB3 + nu.*(4.*(tB9+tB10+tB11)-3.*tB8) - 3.*chi.*gam_u.*tB8 - (M-mu).*(3.*kbc.*gam_u + 3.*kbpc.*gam_p + 3.*kbp.*(1-gam_p)).*tB8);
        B9 = tB9 + dt.*( k6.*gam_u.^2.*tB5 + chi.*gam_u.*(tB6+tB7) + nu.*(5.*tB12-4.*tB9) - k6.*gam_u.^2.*tB9 - chi.*gam_u.*tB9);
        B10 = tB10 + dt.*( k6.*gam_u.^2.*(tB5+tB6) + chi.*gam_u.*(tB7+tB8) + nu.*(5.*tB12-4.*tB10) - 2.*chi.*gam_u.*tB10);
        B11 = tB11 + dt.*( k6.*gam_u.^2.*tB7 + chi.*gam_u.*tB6 + nu.*(5.*tB12-4.*tB11) - 2.*chi.*gam_u.*tB11);
        B12 = tB12 + dt.*( k6.*gam_u.^2.*tB9 + chi.*gam_u.*(tB9+2.*tB10+2.*tB11) + nu.*(6.*tB13-5.*tB12) - chi.*gam_u.*tB12);
        B13 = tB13 + dt.*( chi.*gam_u.*tB12 - nu.*6.*tB13);
        PP1 = tPP1 + dt.*( -k11.*tI1P.*tPP1 + km11.*(PP10 - tPP1));
        I1P = tI1P + dt.*( -k11.*tI1P.*tPP1 + km11.*(PP10 - tPP1) + vPKA_I1.*(I10-tI1P) - vCaN_I1.*tI1P);
        mu = tmu + dt.*( 6.*(M-tmu).*(kbc.*gam_u).*S0 + (M-tmu).*(kbpc.*gam_p + kbp.*(1-gam_p) - kbc.*gam_u).*Sp);

        dCa = (Ca ~= 0) .* (Ca - tCa)./tCa;
        dB1 = (B1 ~= 0) .* (B1 - tB1)./tB1;
        dB2 = (B2 ~= 0) .* (B2 - tB2)./tB2;
        dB3 = (B3 ~= 0) .* (B3 - tB3)./tB3;
        dB4 = (B4 ~= 0) .* (B4 - tB4)./tB4;
        dB5 = (B5 ~= 0) .* (B5 - tB5)./tB5;
        dB6 = (B6 ~= 0) .* (B6 - tB6)./tB6;
        dB7 = (B7 ~= 0) .* (B7 - tB7)./tB7;
        dB8 = (B8 ~= 0) .* (B8 - tB8)./tB8;
        dB9 = (B9 ~= 0) .* (B9 - tB9)./tB9;
        dB10 = (B10 ~= 0) .* (B10 - tB10)./tB10;
        dB11 = (B11 ~= 0) .* (B11 - tB11)./tB11;
        dB12 = (B12 ~= 0) .* (B12 - tB12)./tB12;
        dB13 = (B13 ~= 0) .* (B13 - tB13)./tB13;
        dPP1 = (PP1 ~= 0) .* (PP1 - tPP1)./tPP1;
        dI1P = (I1P ~= 0) .* (I1P - tI1P)./tI1P;
        dmu = (mu ~= 0) * (mu - tmu)/tmu;    

        inspect_id = randi(n_threads);
        maxrelvar = max([dCa, dB1, dB2, dB3, dB4, dB5, dB6, dB7, dB8, dB9, dB10, dB11, dB12, dB13, dPP1, dI1P, dmu]);
        if maxrelvar > 3e-2
            dt = dt/2;
            Ca = tCa;
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
            Cas = [Cas; Ca];
            B1s = [B1s; B1 + sigma.*sqrt(dt.*( 6.*k6.*gam_u.^2.*B0 + 4.*k6.*gam_u.^2.*tB1 + chi.*gam_u.*tB1 + nu.*(2.*(tB2+tB3+tB4)-tB1) )) .* randn(n_threads, 1)];
            B2s = [B2s; B2 + sigma.*sqrt(dt.*( k6.*gam_u.^2.*tB1 + chi.*gam_u.*tB1 + nu.*(3.*(tB5+tB6+tB7+tB8)-2.*tB2) + 3.*k6.*gam_u.^2.*tB2 + chi.*gam_u.*tB2 )) .* randn(n_threads, 1)];
            B3s = [B3s; B3 + sigma.*sqrt(dt.*( 2.*k6.*gam_u.^2.*tB1 + nu.*(3.*(tB5+tB6+tB7+tB8)-2.*tB3) + 3.*k6.*gam_u.^2.*tB3 + chi.*gam_u.*tB3 )) .* randn(n_threads, 1)];
            B4s = [B4s; B4 + sigma.*sqrt(dt.*( k6.*gam_u.^2.*tB1 + nu.*(3.*(tB5+tB6+tB7+tB8)-2.*tB4) + 2.*k6.*gam_u.^2.*tB4 + 2.*chi.*gam_u.*tB4 )) .* randn(n_threads, 1)];
            B5s = [B5s; B5 + sigma.*sqrt(dt.*( k6.*gam_u.^2.*(tB2+tB3) + chi.*gam_u.*tB2 + nu.*(4.*(tB9+tB10+tB11)-3.*tB5) + 2.*k6.*gam_u.^2.*tB5 + chi.*gam_u.*tB5 )) .* randn(n_threads, 1)];
            B6s = [B6s; B6 + sigma.*sqrt(dt.*( k6.*gam_u.^2.*(tB2+tB3) + 2.*chi.*gam_u.*tB4 + nu.*(4.*(tB9+tB10+tB11)-3.*tB6) + k6.*gam_u.^2.*tB6 + 2.*chi.*gam_u.*tB6 )) .* randn(n_threads, 1)];
            B7s = [B7s; B7 + sigma.*sqrt(dt.*( k6.*gam_u.^2.*(tB2+2.*tB4) + chi.*gam_u.*tB3 + nu.*(4.*(tB9+tB10+tB11)-3.*tB7) + k6.*gam_u.^2.*tB7 + 2.*chi.*gam_u.*tB7 )) .* randn(n_threads, 1)];
            B8s = [B8s; B8 + sigma.*sqrt(dt.*( k6.*gam_u.^2.*tB3 + nu.*(4.*(tB9+tB10+tB11)-3.*tB8) + 3.*chi.*gam_u.*tB8 )) .* randn(n_threads, 1)];
            B9s = [B9s; B9 + sigma.*sqrt(dt.*( k6.*gam_u.^2.*tB5 + chi.*gam_u.*(tB6+tB7) + nu.*(5.*tB12-4.*tB9) + k6.*gam_u.^2.*tB9 + chi.*gam_u.*tB9 )) .* randn(n_threads, 1)];
            B10s = [B10s; B10 + sigma.*sqrt(dt.*( k6.*gam_u.^2.*(tB5+tB6) + chi.*gam_u.*(tB7+tB8) + nu.*(5.*tB12-4.*tB10) + 2.*chi.*gam_u.*tB10 )) .* randn(n_threads, 1)];
            B11s = [B11s; B11 + sigma.*sqrt(dt.*( k6.*gam_u.^2.*tB7 + chi.*gam_u.*tB6 + nu.*(5.*tB12-4.*tB11) + 2.*chi.*gam_u.*tB11 )) .* randn(n_threads, 1)];
            B12s = [B12s; B12 + sigma.*sqrt(dt.*( k6.*gam_u.^2.*tB9 + chi.*gam_u.*(tB9+2.*tB10+2.*tB11) + nu.*(6.*tB13-5.*tB12) + chi.*gam_u.*tB12 )) .* randn(n_threads, 1)];
            B13s = [B13s; B13 + sigma.*sqrt(dt.*( chi.*gam_u.*tB12 + nu.*6.*tB13 )) .* randn(n_threads, 1)];
            PP1s = [PP1s; PP1 + sigma.*sqrt(dt.*( k11.*tI1P.*tPP1 + km11.*(PP10-tPP1) )) .* randn(n_threads, 1)];
            I1Ps = [I1Ps; I1P + sigma.*sqrt(dt.*( k11.*tI1P.*tPP1 + km11.*(PP10-tPP1) + vPKA_I1.*(I10-tI1P) + vCaN_I1.*tI1P )) .* randn(n_threads, 1)];
            mus = [mus; mu + sigma.*sqrt(dt.*(6.*(M-tmu).*(kbc.*gam_u).*S0 + (M-tmu).*(kbpc.*gam_p + kbp.*(1-gam_p) + kbc.*gam_u).*Sp )) .* randn(n_threads, 1)];

            % Linear entities update
            S0s = [S0s; B1 + B2 + B3 + B4 + B5 + B6 + B7 + B8 + B9 + B10 + B11 + B12 + B13];
            B0s = [B0s; Stot - S0];

            k10s = [k10s; k12./(KM + Sp)];
            
            gam_us = [gam_us; C./(K5 + C)];
            gam_ps = [gam_ps; C./(K9 + C)];
            Sps = [Sps, B1 + 2.*(B2 + B3 + B4) + 3.*(B5 + B6 + B7 + B8) + 4.*(B9 + B10 + B11) + 5.*B12 + 6.*B13];
            Sus = [Sus, 6.*Stot - Sp];
            
            Cs = [Cs, CaM./(1 + L4./Ca + L3.*L4./(Ca.^2) + L2.*L3.*L4./(Ca.^3) + L1.*L2.*L3.*L4./(Ca.^4))];
            Cbs = [Cbs, gam_u.*Su + gam_p.*Sp];

            chis = [chis, k7.*gam_p + k8.*(1 - gam_p)];
            nus = [nus, k10.*PP1];
            
            AMPA_bnds = [AMPA_bnds; 2.*(B2+B3+B4) + 6.*(B5+B6+B7+B8) + 12.*(B9+B10+B11) + 20.*B12 + 30.*B13];

            vPKA_I1s = [vPKA_I1s; kpka0_I1 + kpka_I1./(1 + (Kdpka./(C-Cb)).^npka)];
            vCaN_I1s = [vCaN_I1s; kcan0_I1 + kcan_I1./(1 + (Kdcan./(C-Cb)).^ncan)];


            % Step computation
            dt = base_step;
            t = t + dt;
            ts = [ts; t];
        end
    end

    y = [Cas, Cs, B0s, B1s, B2s, B3s, B4s, B5s, B6s, B7s, B8s, B9s, B10s, B11s, B12s, B13s, chis, nus, S0s, Sps, Sus, AMPA_bnds, Cbs, vPKA_I1s, vCaN_I1s, PP1s, I1Ps, gam_us, gam_ps, k10s, mus];
end

function K = getC(Ca, CaM, L1, L2, L3, L4)
    K = double(CaM./(1 + L4./Ca + L3.*L4./(Ca.^2) + L2.*L3.*L4./(Ca.^3) + L1.*L2.*L3.*L4./(Ca.^4)));
end


function [p0,i0] = getP0(C, kin_CaN, kin_PKA, par_PP1)
    vCaN = kin_CaN(3) + kin_CaN(4)./(1 + (kin_CaN(1)./C).^kin_CaN(2));
    vPKA = kin_PKA(3) + kin_PKA(4)./(1 + (kin_PKA(1)./C).^kin_PKA(2));

    i0 = min(par_PP1(3).*vPKA./vCaN,par_PP1(3));
    p0 = par_PP1(4)./(1+i0.*par_PP1(1)./par_PP1(2));
end


function samples = uniform_simplex(n_samples, dim)
    a = rand(dim-1,n_samples);
    a = sort(a,1);
    a = [zeros(1,n_samples); a; ones(1,n_samples)];
    samples = a - circshift(a, -1, 1);
    samples = samples(1:end-1,:);
end