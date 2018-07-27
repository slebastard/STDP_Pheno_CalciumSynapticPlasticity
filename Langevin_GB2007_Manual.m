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

varNames = ["Ca", "C", "B0", "B1", "B2", "B3", "B4", "B5", "B6", "B7", "B8", "B9", "B10", "B11", "B12", "B13", "chi", "nu", "S0", "Sp", "Su", "AMPA_bnd", "Cb", "vPKA_I1", "vCaN_I1", "PP1", "I1P", "gam_u", "gam_p", "k10", "mu"];

% INPUTS
in_CaM = 0.1;
savePlots = false;

% VARIABLES AND PARAMETERS
params.tauCa = 0.012;
params.CaBas = 0.4;
params.Stot = 33.3;
params.CaM = in_CaM;
params.K5 = 0.1;
params.K9 = 0.0001;
params.L1 = 0.1;
params.L2 = 0.025;
params.L3 = 0.32;
params.L4 = 0.40;
params.k6 = 6;
params.k7 = 6;
params.k8 = 6;
params.KM = 0.4;
params.k12 = 6000;
params.k11 = 500;
params.km11 = 0.1;
params.I10 = 1;
params.PP10 = 0.2;
params.Kdcan = 0.053;
params.ncan = 3;
params.kcan0_I1 = 0.1;
params.kcan_I1 = 18;
params.Kdpka = 0.11;
params.npka = 8;
params.kpka0_I1 = 0.00359;
params.kpka_I1 = 100;
params.M = 1000;
params.kbc = 0;
params.kbpc = 0;
params.kbp = 0;
params.hgt = 0.5;
params.rad_spn = 0.2;
params.rad_psd = 0.1;
params.NA = 1e-6*6.02e23; % expressed in 1/mumol
params.rate_PSD = 0.05;

V=(1e-15)*params.hgt*pi*params.rad_psd^2;

simu_bootstrap.t0=0; simu_bootstrap.tfinal=40;
simu_bootstrap.base_step = 1;
simu_bootstrap.n_threads = 10;
simu_bootstrap.n_repetes = 1;
simu_bootstrap.sigma = 0; %1/sqrt(V*params.NA);
simu_bootstrap.tol = 1e-2;

simu.t0=0; simu.tfinal=30;
simu.base_step = 1;
simu.n_threads = 10;
simu.n_repetes = 1;
simu.sigma = 0; %1/sqrt(V*params.NA);
simu.tol = 1e-2;

% Sample uniformly and boostrap to get set of natural conditions

Ca_stim = 1*ones(1,simu.n_threads*simu.n_repetes);
Ca_alpha = 0.6*ones(1,simu.n_threads*simu.n_repetes);
t_stim = 2*rand(simu.n_threads,1)';
t_stim = repelem(t_stim, 1, simu.n_repetes);

initConds = params.Stot.*uniform_simplex(simu.n_threads, 14);
initConds = repelem(initConds, 1, simu.n_repetes);

vars.B0 = initConds(1,:);
vars.B1 = initConds(2,:);
vars.B2 = initConds(3,:);
vars.B3 = initConds(4,:);
vars.B4 = initConds(5,:);
vars.B5 = initConds(6,:);
vars.B6 = initConds(7,:);
vars.B7 = initConds(8,:);
vars.B8 = initConds(9,:);
vars.B9 = initConds(10,:);
vars.B10 = initConds(11,:);
vars.B11 = initConds(12,:);
vars.B12 = initConds(13,:);
vars.B13 = initConds(14,:);

[y_init,t_init] = stimulate(Ca_stim, t_stim, params, vars, simu_bootstrap);

% Use those states as initial conditions and stimulate
vars.B0 = y_init.B0(end,:);
vars.B1 = y_init.B1(end,:);
vars.B2 = y_init.B2(end,:);
vars.B3 = y_init.B3(end,:);
vars.B4 = y_init.B4(end,:);
vars.B5 = y_init.B5(end,:);
vars.B6 = y_init.B6(end,:);
vars.B7 = y_init.B7(end,:);
vars.B8 = y_init.B8(end,:);
vars.B9 = y_init.B9(end,:);
vars.B10 = y_init.B10(end,:);
vars.B11 = y_init.B11(end,:);
vars.B12 = y_init.B12(end,:);
vars.B13 = y_init.B13(end,:);

t_stim = 2*rand(simu.n_threads,1)';
t_stim = repelem(t_stim, 1, simu.n_repetes);


[y,t] = stimulate(Ca_stim, t_stim, params, vars, simu);

t_alpha = t_stim + params.tauCa*log(Ca_stim./(Ca_alpha-params.CaBas));

Sp = y.Sp;
id_alpha = sum(t > t_alpha, 1);
Sp_stop = Sp(id_alpha);
Sp_ifty = Sp(end,:);

% What is final strength vs CaMKII level at stop criterion?
% figure(1)
% plot(Sp_stop, Sp_ifty, 'x');
% title('Transfer curve from stopping criterion to stable system');
% xlabel('Phosphorylated CaMKII at stop calcium concentration');
% ylabel('Phosphorylated CaMKII in stable state');

% Plotting randomly selected response
plt_id = 1;
vars = [y.Ca(:, plt_id), ...
y.C(:, plt_id), ...
y.B0(:, plt_id), ...
y.B1(:, plt_id), ...
y.B2(:, plt_id), ...
y.B3(:, plt_id), ...
y.B4(:, plt_id), ...
y.B5(:, plt_id), ...
y.B6(:, plt_id), ...
y.B7(:, plt_id), ...
y.B8(:, plt_id), ...
y.B9(:, plt_id), ...
y.B10(:, plt_id), ...
y.B11(:, plt_id), ...
y.B12(:, plt_id), ...
y.B13(:, plt_id), ...
y.chi(:, plt_id), ...
y.nu(:, plt_id), ...
y.S0(:, plt_id), ...
y.Sp(:, plt_id), ...
y.Su(:, plt_id), ...
y.AMPA_bnd(:, plt_id), ...
y.Cb(:, plt_id), ...
y.vPKA_I1(:, plt_id), ...
y.vCaN_I1(:, plt_id), ...
y.PP1(:, plt_id), ...
y.I1P(:, plt_id), ...
y.gam_u(:, plt_id), ...
y.gam_p(:, plt_id), ...
y.k10(:, plt_id), ...
y.mu(:, plt_id)];

plt_h=4; plt_l=4;
subt = sprintf('Randomly chosen system response');
subplot = @(m,n,p) subtightplot (m, n, p, [0.05 0.035], [0.1 0.01], [0.1 0.01]);

for idx = 1:numel(fieldnames(y))
    if mod(idx,plt_h*plt_l)==1
        ax = subtitle(subt);
        axes(ax);
        fig = figure(2 + fix(idx/(plt_h*plt_l)));
        set(gcf, 'Position', get(0, 'Screensize'));
    end
    h = subplot(plt_h,plt_l,1+mod(idx-1,plt_h*plt_l));
    p = get(h, 'Position');
    set(h,'pos',p+[0 -0.04 0 -0.04]);
    plot(t(:,1),vars(:,idx), 'x')
    title(varNames(idx))
end

ax = subtitle(subt);
axes(ax);

function [y,ts] = stimulate(in_CaInit, tstim, p, v, s)

    B0 = v.B0; B1 = v.B1; B2 = v.B2; B3 = v.B3; B4 = v.B4; B5 = v.B5; B6 = v.B6;
    B7 = v.B7; B8 = v.B8; B9 = v.B9; B10 = v.B10; B11 = v.B11; B12 = v.B12; B13 = v.B13;

    tauCa = p.tauCa; CaBas = p.CaBas; Stot = p.Stot; CaM = p.CaM; K5 = p.K5; K9 = p.K9;
    L1 = p.L1; L2 = p.L2; L3 = p.L3; L4 = p.L4; k6 = p.k6; k7 = p.k7; k8 = p.k8; KM = p.KM; k12 = p.k12;
    k11 = p.k11; km11 = p.km11; I10 = p.I10; PP10 = p.PP10;
    Kdcan = p.Kdcan; ncan = p.ncan; kcan0_I1 = p.kcan0_I1; kcan_I1 = p.kcan_I1; Kdpka = p.Kdpka; npka = p.npka; kpka0_I1 = p.kpka0_I1; kpka_I1 = p.kpka_I1; M = p.M;
    kbc = p.kbc; kbpc = p.kbpc; kbp = p.kbp; hgt = p.hgt; rad_spn = p.rad_spn; rad_psd = p.rad_psd; NA = p.NA; rate_PSD = p.rate_PSD;

    C0 = getC(in_CaInit, CaM, L1, L2, L3, L4);
    [p0,i0] = getP0(C0, [Kdcan,ncan,kcan0_I1,kcan_I1], [Kdpka,npka,kpka0_I1,kpka_I1], [k11,km11,I10,PP10]);

    Ca = in_CaInit; C = C0;
    S0 = B0 + B1 + B2 + B3 + B4 + B5 + B6 + B7 + B8 + B9 + B10 + B11 + B12 + B13;
    gam_u = C0./(K5 + C0); gam_p = C0./(K9 + C0);
    Sp = B1 + 2*(B2 + B3 + B4) + 3*(B5 + B6 + B7 + B8) + 4*(B9 + B10 + B11) + 5*B12 + 6*B13;
    Su = 6*Stot - Sp; AMPA_bnd = 2*(B2+B3+B4) + 6*(B5+B6+B7+B8) + 12*(B9+B10+B11) + 20*B12 + 30*B13; Cb = gam_u.*Su + gam_p.*Sp;
    vPKA_I1 = kpka0_I1 + kpka_I1./(1 + (Kdpka./(C-Cb)).^npka); vCaN_I1 = kcan0_I1 + kcan_I1./(1 + (Kdcan./(C-Cb)).^ncan);
    PP1 = p0; I1P = i0;
    k10 = k12./(KM + Sp);
    chi = (k7*gam_p+k8*(1-gam_p)).*ones(1, s.n_threads*s.n_repetes); nu = k10;
    mu = zeros(1, s.n_threads*s.n_repetes);

    Cas = Ca; Cs = C;
    B0s = B0; B1s = B1; B2s = B2; B3s = B3; B4s = B4; B5s = B5; B6s = B6; B7s = B7; B8s = B8;
    B9s = B9; B10s = B10; B11s = B11; B12s = B12; B13s = B13; chis = chi; nus = nu; S0s = S0;
    Sps = Sp; Sus = Su; AMPA_bnds = AMPA_bnd; Cbs = Cb;
    vPKA_I1s = vPKA_I1; vCaN_I1s = vCaN_I1;
    PP1s = PP1; I1Ps = I1P;
    gam_us = gam_u; gam_ps = gam_p; k10s = k10; 
    mus = mu;

    t = s.t0*ones(1,s.n_threads*s.n_repetes); ts = t;
    dt = s.base_step*ones(1,s.n_threads*s.n_repetes);

    err = 1;
    while any(t < s.tfinal)
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

        stop_dt_eval = 0;    
        while ~stop_dt_eval
            
            dt = 0.9*dt.*min(max((s.tol./abs(err)),0.1),5);

            Ca_0 = tCa + (t>tstim).*dt.*( -1./(tauCa).*(tCa - CaBas));
            B1_0 = tB1 + dt.*( 6.*k6.*gam_u.^2.*B0 - 4.*k6.*gam_u.^2.*tB1 - chi.*gam_u.*tB1 + nu.*(2.*(tB2+tB3+tB4)-tB1));
            B2_0 = tB2 + dt.*( k6.*gam_u.^2.*tB1 + chi.*gam_u.*tB1 + nu.*(3.*(tB5+tB6+tB7+tB8)-2.*tB2) - 3.*k6.*gam_u.^2.*tB2 - chi.*gam_u.*tB2);
            B3_0 = tB3 + dt.*( 2.*k6.*gam_u.^2.*tB1 + nu.*(3.*(tB5+tB6+tB7+tB8)-2.*tB3) - 3.*k6.*gam_u.^2.*tB3 - chi.*gam_u.*tB3);
            B4_0 = tB4 + dt.*( k6.*gam_u.^2.*tB1 + nu.*(3.*(tB5+tB6+tB7+tB8)-2.*tB4) - 2.*k6.*gam_u.^2.*tB4 - 2.*chi.*gam_u.*tB4);
            B5_0 = tB5 + dt.*( k6.*gam_u.^2.*(tB2+tB3) + chi.*gam_u.*tB2 + nu.*(4.*(tB9+tB10+tB11)-3.*tB5) - 2.*k6.*gam_u.^2.*tB5 - chi.*gam_u.*tB5);
            B6_0 = tB6 + dt.*( k6.*gam_u.^2.*(tB2+tB3) + 2.*chi.*gam_u.*tB4 + nu.*(4.*(tB9+tB10+tB11)-3.*tB6) - k6.*gam_u.^2.*tB6 - 2.*chi.*gam_u.*tB6);
            B7_0 = tB7 + dt.*( k6.*gam_u.^2.*(tB2+2.*tB4) + chi.*gam_u.*tB3 + nu.*(4.*(tB9+tB10+tB11)-3.*tB7) - k6.*gam_u.^2.*tB7 - 2.*chi.*gam_u.*tB7);
            B8_0 = tB8 + dt.*( k6.*gam_u.^2.*tB3 + nu.*(4.*(tB9+tB10+tB11)-3.*tB8) - 3.*chi.*gam_u.*tB8 );
            B9_0 = tB9 + dt.*( k6.*gam_u.^2.*tB5 + chi.*gam_u.*(tB6+tB7) + nu.*(5.*tB12-4.*tB9) - k6.*gam_u.^2.*tB9 - chi.*gam_u.*tB9);
            B10_0 = tB10 + dt.*( k6.*gam_u.^2.*(tB5+tB6) + chi.*gam_u.*(tB7+tB8) + nu.*(5.*tB12-4.*tB10) - 2.*chi.*gam_u.*tB10);
            B11_0 = tB11 + dt.*( k6.*gam_u.^2.*tB7 + chi.*gam_u.*tB6 + nu.*(5.*tB12-4.*tB11) - 2.*chi.*gam_u.*tB11);
            B12_0 = tB12 + dt.*( k6.*gam_u.^2.*tB9 + chi.*gam_u.*(tB9+2.*tB10+2.*tB11) + nu.*(6.*tB13-5.*tB12) - chi.*gam_u.*tB12);
            B13_0 = tB13 + dt.*( chi.*gam_u.*tB12 - nu.*6.*tB13);
            PP1_0 = tPP1 + dt.*( -k11.*tI1P.*tPP1 + km11.*(PP10 - tPP1));
            I1P_0 = tI1P + dt.*( -k11.*tI1P.*tPP1 + km11.*(PP10 - tPP1) + vPKA_I1.*(I10-tI1P) - vCaN_I1.*tI1P);
            mu_0 = tmu + dt.*( 6.*(M-tmu).*(kbc.*gam_u).*S0 + (M-tmu).*(kbpc.*gam_p + kbp.*(1-gam_p) - kbc.*gam_u).*Sp);
    
            C_0 = CaM./(1 + L4./Ca_0 + L3.*L4./(Ca_0.^2) + L2.*L3.*L4./(Ca_0.^3) + L1.*L2.*L3.*L4./(Ca_0.^4));
            gam_u_0 = C_0./(K5 + C_0);
            gam_p_0 = C_0./(K9 + C_0); 
            S0_0 = B1_0 + B2_0 + B3_0 + B4_0 + B5_0 + B6_0 + B7_0 + B8_0 + B9_0 + B10_0 + B11_0 + B12_0 + B13_0;
            B0_0 = Stot - S0_0;
            Sp_0 = B1_0 + 2.*(B2_0 + B3_0 + B4_0) + 3.*(B5_0 + B6_0 + B7_0 + B8_0) + 4.*(B9_0 + B10_0 + B11_0) + 5.*B12_0 + 6.*B13_0;
            Su_0 = 6.*Stot - Sp_0; 
            k10_0 = k12./(KM + Sp_0);   
            Cb_0 = gam_u_0.*Su_0 + gam_p_0.*Sp_0;
            chi_0 = k7.*gam_p_0 + k8.*(1 - gam_p_0);
            nu_0 = k10_0.*PP1_0;
            AMPA_bnd_0 = 2.*(B2_0+B3_0+B4_0) + 6.*(B5_0+B6_0+B7_0+B8_0) + 12.*(B9_0+B10_0+B11_0) + 20.*B12_0 + 30.*B13_0;
            vPKA_I1_0 = kpka0_I1 + kpka_I1./(1 + (Kdpka./(C_0-Cb_0)).^npka);
            vCaN_I1_0 = kcan0_I1 + kcan_I1./(1 + (Kdcan./(C_0-Cb_0)).^ncan);

            Ca_demi = Ca + 0.5*(t>tstim).*dt.*( -1./(tauCa).*(tCa - CaBas));
            B1_demi = tB1 + 0.5*dt.*( 6.*k6.*gam_u.^2.*B0 - 4.*k6.*gam_u.^2.*tB1 - chi.*gam_u.*tB1 + nu.*(2.*(tB2+tB3+tB4)-tB1));
            B2_demi = tB2 + 0.5*dt.*( k6.*gam_u.^2.*tB1 + chi.*gam_u.*tB1 + nu.*(3.*(tB5+tB6+tB7+tB8)-2.*tB2) - 3.*k6.*gam_u.^2.*tB2 - chi.*gam_u.*tB2);
            B3_demi = tB3 + 0.5*dt.*( 2.*k6.*gam_u.^2.*tB1 + nu.*(3.*(tB5+tB6+tB7+tB8)-2.*tB3) - 3.*k6.*gam_u.^2.*tB3 - chi.*gam_u.*tB3);
            B4_demi = tB4 + 0.5*dt.*( k6.*gam_u.^2.*tB1 + nu.*(3.*(tB5+tB6+tB7+tB8)-2.*tB4) - 2.*k6.*gam_u.^2.*tB4 - 2.*chi.*gam_u.*tB4);
            B5_demi = tB5 + 0.5*dt.*( k6.*gam_u.^2.*(tB2+tB3) + chi.*gam_u.*tB2 + nu.*(4.*(tB9+tB10+tB11)-3.*tB5) - 2.*k6.*gam_u.^2.*tB5 - chi.*gam_u.*tB5);
            B6_demi = tB6 + 0.5*dt.*( k6.*gam_u.^2.*(tB2+tB3) + 2.*chi.*gam_u.*tB4 + nu.*(4.*(tB9+tB10+tB11)-3.*tB6) - k6.*gam_u.^2.*tB6 - 2.*chi.*gam_u.*tB6);
            B7_demi = tB7 + 0.5*dt.*( k6.*gam_u.^2.*(tB2+2.*tB4) + chi.*gam_u.*tB3 + nu.*(4.*(tB9+tB10+tB11)-3.*tB7) - k6.*gam_u.^2.*tB7 - 2.*chi.*gam_u.*tB7);
            B8_demi = tB8 + 0.5*dt.*( k6.*gam_u.^2.*tB3 + nu.*(4.*(tB9+tB10+tB11)-3.*tB8) - 3.*chi.*gam_u.*tB8 );
            B9_demi = tB9 + 0.5*dt.*( k6.*gam_u.^2.*tB5 + chi.*gam_u.*(tB6+tB7) + nu.*(5.*tB12-4.*tB9) - k6.*gam_u.^2.*tB9 - chi.*gam_u.*tB9);
            B10_demi = tB10 + 0.5*dt.*( k6.*gam_u.^2.*(tB5+tB6) + chi.*gam_u.*(tB7+tB8) + nu.*(5.*tB12-4.*tB10) - 2.*chi.*gam_u.*tB10);
            B11_demi = tB11 + 0.5*dt.*( k6.*gam_u.^2.*tB7 + chi.*gam_u.*tB6 + nu.*(5.*tB12-4.*tB11) - 2.*chi.*gam_u.*tB11);
            B12_demi = tB12 + 0.5*dt.*( k6.*gam_u.^2.*tB9 + chi.*gam_u.*(tB9+2.*tB10+2.*tB11) + nu.*(6.*tB13-5.*tB12) - chi.*gam_u.*tB12);
            B13_demi = tB13 + 0.5*dt.*( chi.*gam_u.*tB12 - nu.*6.*tB13);
            PP1_demi = tPP1 + 0.5*dt.*( -k11.*tI1P.*tPP1 + km11.*(PP10 - tPP1));
            I1P_demi = tI1P + 0.5*dt.*( -k11.*tI1P.*tPP1 + km11.*(PP10 - tPP1) + vPKA_I1.*(I10-tI1P) - vCaN_I1.*tI1P);
            mu_demi = tmu + 0.5*dt.*( 6.*(M-tmu).*(kbc.*gam_u).*S0 + (M-tmu).*(kbpc.*gam_p + kbp.*(1-gam_p) - kbc.*gam_u).*Sp); 
    
            C_demi = CaM./(1 + L4./Ca_demi + L3.*L4./(Ca_demi.^2) + L2.*L3.*L4./(Ca_demi.^3) + L1.*L2.*L3.*L4./(Ca_demi.^4));
            gam_u_demi = C_demi./(K5 + C_demi);
            gam_p_demi = C_demi./(K9 + C_demi); 
            S0_demi = B1_demi + B2_demi + B3_demi + B4_demi + B5_demi + B6_demi + B7_demi + B8_demi + B9_demi + B10_demi + B11_demi + B12_demi + B13_demi;
            B0_demi = Stot - S0_demi;
            Sp_demi = B1_demi + 2.*(B2_demi + B3_demi + B4_demi) + 3.*(B5_demi + B6_demi + B7_demi + B8_demi) + 4.*(B9_demi + B10_demi + B11_demi) + 5.*B12_demi + 6.*B13_demi;
            Su_demi = 6.*Stot - Sp_demi; 
            k10_demi = k12./(KM + Sp_demi);   
            Cb_demi = gam_u_demi.*Su_demi + gam_p_demi.*Sp_demi;
            chi_demi = k7.*gam_p_demi + k8.*(1 - gam_p_demi);
            nu_demi = k10_demi.*PP1_demi;
            AMPA_bnd_demi = 2.*(B2_demi+B3_demi+B4_demi) + 6.*(B5_demi+B6_demi+B7_demi+B8_demi) + 12.*(B9_demi+B10_demi+B11_demi) + 20.*B12_demi + 30.*B13_demi;
            vPKA_I1_demi = kpka0_I1 + kpka_I1./(1 + (Kdpka./(C_demi-Cb_demi)).^npka);
            vCaN_I1_demi = kcan0_I1 + kcan_I1./(1 + (Kdcan./(C_demi-Cb_demi)).^ncan);
    
            Ca_1 = Ca_demi + 0.5*(t>tstim).*dt.*( -1./(tauCa).*(Ca_demi - CaBas));
            B1_1 = B1_demi + 0.5*dt.*( 6.*k6.*gam_u_demi.^2.*B0_demi - 4.*k6.*gam_u_demi.^2.*B1_demi - chi_demi.*gam_u_demi.*B1_demi + nu_demi.*(2.*(B2_demi+B3_demi+B4_demi)-B1_demi));
            B2_1 = B2_demi + 0.5*dt.*( k6.*gam_u_demi.^2.*B1_demi + chi_demi.*gam_u_demi.*B1_demi + nu_demi.*(3.*(B5_demi+B6_demi+B7_demi+B8_demi)-2.*B2_demi) - 3.*k6.*gam_u_demi.^2.*B2_demi - chi_demi.*gam_u_demi.*B2_demi);
            B3_1 = B3_demi + 0.5*dt.*( 2.*k6.*gam_u_demi.^2.*B1_demi + nu_demi.*(3.*(B5_demi+B6_demi+B7_demi+B8_demi)-2.*B3_demi) - 3.*k6.*gam_u_demi.^2.*B3_demi - chi_demi.*gam_u_demi.*B3_demi);
            B4_1 = B4_demi + 0.5*dt.*( k6.*gam_u_demi.^2.*B1_demi + nu_demi.*(3.*(B5_demi+B6_demi+B7_demi+B8_demi)-2.*B4_demi) - 2.*k6.*gam_u_demi.^2.*B4_demi - 2.*chi_demi.*gam_u_demi.*B4_demi);
            B5_1 = B5_demi + 0.5*dt.*( k6.*gam_u_demi.^2.*(B2_demi+B3_demi) + chi_demi.*gam_u_demi.*B2_demi + nu_demi.*(4.*(B9_demi+B10_demi+B11_demi)-3.*B5_demi) - 2.*k6.*gam_u_demi.^2.*B5_demi - chi_demi.*gam_u_demi.*B5_demi);
            B6_1 = B6_demi + 0.5*dt.*( k6.*gam_u_demi.^2.*(B2_demi+B3_demi) + 2.*chi_demi.*gam_u_demi.*B4_demi + nu_demi.*(4.*(B9_demi+B10_demi+B11_demi)-3.*B6_demi) - k6.*gam_u_demi.^2.*B6_demi - 2.*chi_demi.*gam_u_demi.*B6_demi);
            B7_1 = B7_demi + 0.5*dt.*( k6.*gam_u_demi.^2.*(B2_demi+2.*B4_demi) + chi_demi.*gam_u_demi.*B3_demi + nu_demi.*(4.*(B9_demi+B10_demi+B11_demi)-3.*B7_demi) - k6.*gam_u_demi.^2.*B7_demi - 2.*chi_demi.*gam_u_demi.*B7_demi);
            B8_1 = B8_demi + 0.5*dt.*( k6.*gam_u_demi.^2.*B3_demi + nu_demi.*(4.*(B9_demi+B10_demi+B11_demi)-3.*B8_demi) - 3.*chi_demi.*gam_u_demi.*B8_demi );
            B9_1 = B9_demi + 0.5*dt.*( k6.*gam_u_demi.^2.*B5_demi + chi_demi.*gam_u_demi.*(B6_demi+B7_demi) + nu_demi.*(5.*B12_demi-4.*B9_demi) - k6.*gam_u_demi.^2.*B9_demi - chi_demi.*gam_u_demi.*B9_demi);
            B10_1 = B10_demi + 0.5*dt.*( k6.*gam_u_demi.^2.*(B5_demi+B6_demi) + chi_demi.*gam_u_demi.*(B7_demi+B8_demi) + nu_demi.*(5.*B12_demi-4.*B10_demi) - 2.*chi_demi.*gam_u_demi.*B10_demi);
            B11_1 = B11_demi + 0.5*dt.*( k6.*gam_u_demi.^2.*B7_demi + chi_demi.*gam_u_demi.*B6_demi + nu_demi.*(5.*B12_demi-4.*B11_demi) - 2.*chi_demi.*gam_u_demi.*B11_demi);
            B12_1 = B12_demi + 0.5*dt.*( k6.*gam_u_demi.^2.*B9_demi + chi_demi.*gam_u_demi.*(B9_demi+2.*B10_demi+2.*B11_demi) + nu_demi.*(6.*B13_demi-5.*B12_demi) - chi_demi.*gam_u_demi.*B12_demi);
            B13_1 = B13_demi + 0.5*dt.*( chi_demi.*gam_u_demi.*B12_demi - nu_demi.*6.*B13_demi);
            PP1_1 = PP1_demi + 0.5*dt.*( -k11.*I1P_demi.*PP1_demi + km11.*(PP10 - PP1_demi));
            I1P_1 = I1P_demi + 0.5*dt.*( -k11.*I1P_demi.*PP1_demi + km11.*(PP10 - PP1_demi) + vPKA_I1_demi.*(I10-I1P_demi) - vCaN_I1_demi.*I1P_demi);
            mu_1 = mu_demi + 0.5*dt.*( 6.*(M-mu_demi).*(kbc.*gam_u_demi).*S0_demi + (M-mu_demi).*(kbpc.*gam_p_demi + kbp.*(1-gam_p_demi) - kbc.*gam_u_demi).*Sp_demi);
    
            err_Ca = abs((Ca>0).*(Ca_1 - Ca_0)./Ca);
            err_B1 = abs((B1>0).*(B1_1 - B1_0)./B1);
            err_B2 = abs((B2>0).*(B2_1 - B2_0)./B2);
            err_B3 = abs((B3>0).*(B3_1 - B3_0)./B3);
            err_B4 = abs((B4>0).*(B4_1 - B4_0)./B4);
            err_B5 = abs((B5>0).*(B5_1 - B5_0)./B5);
            err_B6 = abs((B6>0).*(B6_1 - B6_0)./B6);
            err_B7 = abs((B7>0).*(B7_1 - B7_0)./B7);
            err_B8 = abs((B8>0).*(B8_1 - B8_0)./B8);
            err_B9 = abs((B9>0).*(B9_1 - B9_0)./B9);
            err_B10 = abs((B10>0).*(B10_1 - B10_0)./B10);
            err_B11 = abs((B11>0).*(B11_1 - B11_0)./B11);
            err_B12 = abs((B12>0).*(B12_1 - B12_0)./B12);
            err_B13 = abs((B13>0).*(B13_1 - B13_0)./B13);
            err_PP1 = abs((PP1>0).*(PP1_1 - PP1_0)./PP1);
            err_I1P = abs((I1P>0).*(I1P_1 - I1P_0)./I1P);
            err_mu = zeros(1,s.n_threads*s.n_repetes);
    
            err = max([err_Ca; err_B1; err_B2; err_B3; err_B4; err_B5; err_B6; err_B7; err_B8; err_B9; err_B10; err_B11; err_B12; err_B13; err_PP1; err_I1P; err_mu]);
            stop_dt_eval = ~any(err > s.tol);
        end

        Ca = Ca_1;
        B1 = B1_1 + s.sigma.*sqrt(dt.*( 6.*k6.*gam_u.^2.*B0 + 4.*k6.*gam_u.^2.*tB1 + chi.*gam_u.*tB1 + nu.*(2.*(tB2+tB3+tB4)-tB1))) .* randn(1, s.n_threads*s.n_repetes);
        B2 = B2_1 + s.sigma.*sqrt(dt.*( k6.*gam_u.^2.*tB1 + chi.*gam_u.*tB1 + nu.*(3.*(tB5+tB6+tB7+tB8)-2.*tB2) + 3.*k6.*gam_u.^2.*tB2 + chi.*gam_u.*tB2)) .* randn(1, s.n_threads*s.n_repetes);
        B3 = B3_1 + s.sigma.*sqrt(dt.*( 2.*k6.*gam_u.^2.*tB1 + nu.*(3.*(tB5+tB6+tB7+tB8)-2.*tB3) + 3.*k6.*gam_u.^2.*tB3 + chi.*gam_u.*tB3)) .* randn(1, s.n_threads*s.n_repetes);
        B4 = B4_1 + s.sigma.*sqrt(dt.*( k6.*gam_u.^2.*tB1 + nu.*(3.*(tB5+tB6+tB7+tB8)-2.*tB4) + 2.*k6.*gam_u.^2.*tB4 + 2.*chi.*gam_u.*tB4)) .* randn(1, s.n_threads*s.n_repetes);
        B5 = B5_1 + s.sigma.*sqrt(dt.*( k6.*gam_u.^2.*(tB2+tB3) + chi.*gam_u.*tB2 + nu.*(4.*(tB9+tB10+tB11)-3.*tB5) + 2.*k6.*gam_u.^2.*tB5 + chi.*gam_u.*tB5)) .* randn(1, s.n_threads*s.n_repetes);
        B6 = B6_1 + s.sigma.*sqrt(dt.*( k6.*gam_u.^2.*(tB2+tB3) + 2.*chi.*gam_u.*tB4 + nu.*(4.*(tB9+tB10+tB11)-3.*tB6) + k6.*gam_u.^2.*tB6 + 2.*chi.*gam_u.*tB6)) .* randn(1, s.n_threads*s.n_repetes);
        B7 = B7_1 + s.sigma.*sqrt(dt.*( k6.*gam_u.^2.*(tB2+2.*tB4) + chi.*gam_u.*tB3 + nu.*(4.*(tB9+tB10+tB11)-3.*tB7) + k6.*gam_u.^2.*tB7 + 2.*chi.*gam_u.*tB7)) .* randn(1, s.n_threads*s.n_repetes);
        B8 = B8_1 + s.sigma.*sqrt(dt.*( k6.*gam_u.^2.*tB3 + nu.*(4.*(tB9+tB10+tB11)-3.*tB8) + 3.*chi.*gam_u.*tB8 )) .* randn(1, s.n_threads*s.n_repetes);
        B9 = B9_1 + s.sigma.*sqrt(dt.*( k6.*gam_u.^2.*tB5 + chi.*gam_u.*(tB6+tB7) + nu.*(5.*tB12-4.*tB9) + k6.*gam_u.^2.*tB9 - chi.*gam_u.*tB9)) .* randn(1, s.n_threads*s.n_repetes);
        B10 = B10_1 + s.sigma.*sqrt(dt.*( k6.*gam_u.^2.*(tB5+tB6) + chi.*gam_u.*(tB7+tB8) + nu.*(5.*tB12-4.*tB10) + 2.*chi.*gam_u.*tB10)) .* randn(1, s.n_threads*s.n_repetes);
        B11 = B11_1 + s.sigma.*sqrt(dt.*( k6.*gam_u.^2.*tB7 + chi.*gam_u.*tB6 + nu.*(5.*tB12-4.*tB11) + 2.*chi.*gam_u.*tB11)) .* randn(1, s.n_threads*s.n_repetes);
        B12 = B12_1 + s.sigma.*sqrt(dt.*( k6.*gam_u.^2.*tB9 + chi.*gam_u.*(tB9+2.*tB10+2.*tB11) + nu.*(6.*tB13-5.*tB12) + chi.*gam_u.*tB12)) .* randn(1, s.n_threads*s.n_repetes);
        B13 = B13_1 + s.sigma.*sqrt(dt.*( chi.*gam_u.*tB12 + nu.*6.*tB13)) .* randn(1, s.n_threads*s.n_repetes);
        PP1 = PP1_1 + s.sigma.*sqrt(dt.*( k11.*tI1P.*tPP1 + km11.*(PP10-tPP1))) .* randn(1, s.n_threads*s.n_repetes);
        I1P = I1P_1 + s.sigma.*sqrt(dt.*( k11.*tI1P.*tPP1 + km11.*(PP10-tPP1) + vPKA_I1.*(I10-tI1P) + vCaN_I1.*tI1P)) .* randn(1, s.n_threads*s.n_repetes);
        mu = mu_1 + s.sigma.*sqrt(dt.*( 6.*(M-tmu).*(kbc.*gam_u).*S0 + (M-tmu).*(kbpc.*gam_p + kbp.*(1-gam_p) + kbc.*gam_u).*Sp)) .* randn(1, s.n_threads*s.n_repetes);
        
        Cas = [Cas; Ca];
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

        C = CaM./(1 + L4./Ca + L3.*L4./(Ca.^2) + L2.*L3.*L4./(Ca.^3) + L1.*L2.*L3.*L4./(Ca.^4));
        gam_u = C./(K5 + C);
        gam_p = C./(K9 + C); 
        S0 = B1 + B2 + B3 + B4 + B5 + B6 + B7 + B8 + B9 + B10 + B11 + B12 + B13;
        B0 = Stot - S0;
        Sp = B1 + 2.*(B2 + B3 + B4) + 3.*(B5 + B6 + B7 + B8) + 4.*(B9 + B10 + B11) + 5.*B12 + 6.*B13;
        Su = 6.*Stot - Sp; 
        k10 = k12./(KM + Sp);   
        Cb = gam_u.*Su + gam_p.*Sp;
        chi = k7.*gam_p + k8.*(1 - gam_p);
        nu = k10.*PP1;
        AMPA_bnd = 2.*(B2+B3+B4) + 6.*(B5+B6+B7+B8) + 12.*(B9+B10+B11) + 20.*B12 + 30.*B13;
        vPKA_I1 = kpka0_I1 + kpka_I1./(1 + (Kdpka./(C-Cb)).^npka);
        vCaN_I1 = kcan0_I1 + kcan_I1./(1 + (Kdcan./(C-Cb)).^ncan);

        Cs = [Cs; C];
        gam_us = [gam_us; gam_u];
        gam_ps = [gam_ps; gam_p];
        S0s = [S0s; S0];
        B0s = [B0s; B0];
        Sps = [Sps; Sp];
        Sus = [Sus; Su];
        k10s = [k10s; k10];
        Cbs = [Cbs; Cb];
        chis = [chis; chi];
        nus = [nus; nu];
        AMPA_bnds = [AMPA_bnds; AMPA_bnd];
        vPKA_I1s = [vPKA_I1s; vPKA_I1];
        vCaN_I1s = [vCaN_I1s; vCaN_I1];


        % Step computation
        t = t + dt;
        ts = [ts; t];
    end

    y.Ca = Cas;
    y.C = Cs;
    y.B0 = B0s;
    y.B1 = B1s;
    y.B2 = B2s;
    y.B3 = B3s;
    y.B4 = B4s;
    y.B5 = B5s;
    y.B6 = B6s;
    y.B7 = B7s;
    y.B8 = B8s;
    y.B9 = B9s;
    y.B10 = B10s;
    y.B11 = B11s;
    y.B12 = B12s;
    y.B13 = B13s;
    y.chi = chis;
    y.nu = nus;
    y.S0 = S0s;
    y.Sp = Sps;
    y.Su = Sus;
    y.AMPA_bnd = AMPA_bnds;
    y.Cb = Cbs;
    y.vPKA_I1 = vPKA_I1s;
    y.vCaN_I1 = vCaN_I1s;
    y.PP1 = PP1s;
    y.I1P = I1Ps;
    y.gam_u = gam_us;
    y.gam_p = gam_ps;
    y.k10 = k10s;
    y.mu = mus;
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
    samples = circshift(a, -1, 1) - a;
    samples = samples(1:end-1,:);
end