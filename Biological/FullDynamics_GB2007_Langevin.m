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
% ID    Name        Description
% 1     Ca          Free calcium concentration
% 2     B0          Concentration of olygomer 0 (fully dephosphorylated)
% 3     B1          Concentration of olygomer 1
% 4     B2          Concentration of olygomer 2
% 5     B3          Concentration of olygomer 3
% 6     B4          Concentration of olygomer 4
% 7     B5          Concentration of olygomer 5
% 8     B6          Concentration of olygomer 6
% 9     B7          Concentration of olygomer 7
% 10    B8          Concentration of olygomer 8
% 11    B9          Concentration of olygomer 9
% 12    B10         Concentration of olygomer 10
% 13    B11         Concentration of olygomer 11
% 14    B12         Concentration of olygomer 12 
% 15    B13         Concentration of olygomer 13
% 16    PP1         Concentration of active PP1
% 17    I1P         Concentration of phosphorylated PP inhibitor 1

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

% QUANTITIES AND EQUATIONS

function [eqs_det, eqs_stoch] = getEqs(vars, params)  
    C = @(vars, params) (params(4)/(1 + params(10)/vars(1) + params(9)*params(10)/(vars(1)^2) + params(8)*params(9)*params(10)/(vars(1)^3) + params(7)*params(8)*params(9)*params(10)/(vars(1)^4)));
    S0 = @(vars, params) (vars(3) + vars(4) + vars(5) + vars(6) + vars(7) + vars(8) + vars(9) + vars(10 + vars(11 + vars(12 + vars(13 + vars(14 + vars(15));
    S1 = @(vars, params) ()
    S2 = @(vars, params) ()
    B0 = @(vars, params) (params(3) - S0);
    k10 = @(vars, params) (params(15)/(params(14) + Sp));
    gam_u = @(vars, params) (C/(params(5) + C));
    gam_p = @(vars, params) (C/(params(6) + C));
    Sp = @(vars, params) (vars(3) + 2*(vars(4) + vars(5) + vars(6)) + 3*(vars(7) + vars(8) + vars(9) + vars(10)) + 4*(vars(11) + vars(12) + vars(13)) + 5*vars(14) + 6*vars(15));
    Su = @(vars, params) (6*params(3) - Sp);
    Cb = @(vars, params) (gam_u*Su + gam_p*Sp);
    chi = @(vars, params) (params(12)*gam_p + params(13)*(1 - gam_p));
    nu = @(vars, params) (k10*vars(16));
    AMPA_bnd = @(vars, params) (2*(vars(4) + vars(5) + vars(6)) ...
    + 6*(vars(7) + vars(8) + vars(9) + vars(10)) ...
    + 12*(vars(11) + vars(12) + vars(13)) ...
    + 20*vars(14) + 30*vars(15));
    vPKA_I1 = @(vars, params) (params(26) + params(27)/(1 + (params(24)/(C-Cb))^params(25)));
    vCaN_I1 = @(vars, params) (params(22) + params(23)/(1 + (params(20)/(C-Cb))^params(21)));

    eqs_det = @(vars, params) [
        -1/(params(1))*(vars(1) - params(2));
        6*params(11)*gam_u^2*vars(2) - 4*params(11)*gam_u^2*vars(3) - chi*gam_u*vars(3) + nu*(2*(vars(4)+vars(5)+vars(6))-vars(3));
        params(11)*gam_u^2*vars(3) + chi*gam_u*vars(3) + nu*(3*(vars(7)+vars(8)+vars(9)+vars(10))-2*vars(4)) - 3*params(11)*gam_u^2*vars(4) - chi*gam_u*vars(4);
        2*params(11)*gam_u^2*vars(3) + nu*(3*(vars(7)+vars(8)+vars(9)+vars(10))-2*vars(5)) - 3*params(11)*gam_u^2*vars(5) - chi*gam_u*vars(5);
        params(11)*gam_u^2*vars(3) + nu*(3*(vars(7)+vars(8)+vars(9)+vars(10))-2*vars(6)) - 2*params(11)*gam_u^2*vars(6) - 2*chi*gam_u*vars(6);
        params(11)*gam_u^2*(vars(4)+vars(5)) + chi*gam_u*vars(4) + nu*(4*(vars(11)+vars(12)+vars(13))-3*vars(7)) - 2*params(11)*gam_u^2*vars(7) - chi*gam_u*vars(7);
        params(11)*gam_u^2*(vars(4)+vars(5)) + 2*chi*gam_u*vars(6) + nu*(4*(vars(11)+vars(12)+vars(13))-3*vars(8)) - params(11)*gam_u^2*vars(8) - 2*chi*gam_u*vars(8);
        params(11)*gam_u^2*(vars(4)+2*vars(6)) + chi*gam_u*vars(5) + nu*(4*(vars(11)+vars(12)+vars(13))-3*vars(9)) - params(11)*gam_u^2*vars(9) - 2*chi*gam_u*vars(9);
        params(11)*gam_u^2*vars(5) + nu*(4*(vars(11)+vars(12)+vars(13))-3*vars(10)) - 3*chi*gam_u*vars(10) - (params(28)-mu)*(3*params(29)*gam_u + 3*params(30)*gam_p + 3*params(31)*(1-gam_p))*vars(10);
        params(11)*gam_u^2*vars(7) + chi*gam_u*(vars(8)+vars(9)) + nu*(5*vars(14)-4*vars(11)) - params(11)*gam_u^2*vars(11) - chi*gam_u*vars(11);
        params(11)*gam_u^2*(vars(7)+vars(8)) + chi*gam_u*(vars(9)+vars(10)) + nu*(5*vars(14)-4*vars(12)) - 2*chi*gam_u*vars(12);
        params(11)*gam_u^2*vars(9) + chi*gam_u*vars(8) + nu*(5*vars(14)-4*vars(13)) - 2*chi*gam_u*vars(13);
        params(11)*gam_u^2*vars(11) + chi*gam_u*(vars(11)+2*vars(12)+2*vars(13)) + nu*(6*vars(15)-5*vars(14)) - chi*gam_u*vars(14);
        chi*gam_u*vars(14) - nu*6*vars(15);
        -k11*vars(17)*vars(16) + km11*(vars(16)0 - vars(16));
        -k11*vars(17)*vars(16) + km11*(vars(16)0 - vars(16)) + vPKA_I1*(I10-vars(17)) - vCaN_I1*vars(17);
        6*(params(28)-mu)*(params(29)*gam_u)*S0 + (params(28)-mu)*(params(30)*gam_p + params(31)*(1-gam_p) - params(29)*gam_u)*S1;
    ];

    eqs_stoch = @(vars, params)[
        0;
        ;
        sqrt(6*params(11)*gam_u^2*vars(2) + 4*params(11)*gam_u^2*vars(3) + chi*gam_u*vars(3) + nu*(2*(vars(4)+vars(5)+vars(6))+vars(3)));
        sqrt(params(11)*gam_u^2*vars(3) + chi*gam_u*vars(3) + nu*(3*(vars(7)+vars(8)+vars(9)+vars(10))+2*vars(4)) + 3*params(11)*gam_u^2*vars(4) + chi*gam_u*vars(4));
        sqrt(2*params(11)*gam_u^2*vars(3) + nu*(3*(vars(7)+vars(8)+vars(9)+vars(10))+2*vars(5)) + 3*params(11)*gam_u^2*vars(5) + chi*gam_u*vars(5));
        sqrt(params(11)*gam_u^2*vars(3) + nu*(3*(vars(7)+vars(8)+vars(9)+vars(10))+2*vars(6)) + 2*params(11)*gam_u^2*vars(6) + 2*chi*gam_u*vars(6));
        sqrt(params(11)*gam_u^2*(vars(4)+vars(5)) + chi*gam_u*vars(4) + nu*(4*(vars(11)+vars(12)+vars(13))+3*vars(7)) + 2*params(11)*gam_u^2*vars(7) + chi*gam_u*vars(7));
        sqrt(params(11)*gam_u^2*(vars(4)+vars(5)) + 2*chi*gam_u*vars(6) + nu*(4*(vars(11)+vars(12)+vars(13))+3*vars(8)) + params(11)*gam_u^2*vars(8) + 2*chi*gam_u*vars(8));
        sqrt(params(11)*gam_u^2*(vars(4)+2*vars(6)) + chi*gam_u*vars(5) + nu*(4*(vars(11)+vars(12)+vars(13))+3*vars(9)) + params(11)*gam_u^2*vars(9) + 2*chi*gam_u*vars(9));
        sqrt(params(11)*gam_u^2*vars(5) + nu*(4*(vars(11)+vars(12)+vars(13))+3*vars(10)) + 3*chi*gam_u*vars(10) + (params(28)+mu)*(3*params(29)*gam_u + 3*params(30)*gam_p + 3*params(31)*(1+gam_p))*vars(10));
        sqrt(params(11)*gam_u^2*vars(7) + chi*gam_u*(vars(8)+vars(9)) + nu*(5*vars(14)+4*vars(11)) + params(11)*gam_u^2*vars(11) + chi*gam_u*vars(11));
        sqrt(params(11)*gam_u^2*(vars(7)+vars(8)) + chi*gam_u*(vars(9)+vars(10)) + nu*(5*vars(14)+4*vars(12)) + 2*chi*gam_u*vars(12));
        sqrt(params(11)*gam_u^2*vars(9) + chi*gam_u*vars(8) + nu*(5*vars(14)+4*vars(13)) + 2*chi*gam_u*vars(13));
        sqrt(params(11)*gam_u^2*vars(11) + chi*gam_u*(vars(11)+2*vars(12)+2*vars(13)) + nu*(6*vars(15)+5*vars(14)) + chi*gam_u*vars(14));
        sqrt(chi*gam_u*vars(14) + nu*6*vars(15));
        ;
        ;

    ];
end

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

y0 = [
    in_CaInit; C0;
    33.3; 0; 0; 0; 0; 0; 0; 0; 0;
    0; 0; 0; 0; 0; 0; 0; 0;
    0; 199.8; 0; 0;
    0; 0;
    p0; i0;
    0.5; 0.5; 10000;
    0
];

[eqs_det, eqs_stoch] = getEqs(vars, params);
[eqs_det, eqs_stoch] = @(vars) [eqs_det(vars, paramVals), eqs_stoch(vars, paramVals)];
 
opts = sdeset('NonNegative',17,...
              'RandSeed',1,...
              'SDEType','Ito');

[tSimu,yCaMKII] = sde_euler(eqs_det, eqs_stoch, [t0, tfinal], y0', opts);

% Compute NMDA variables

Bp0 = floor(max(0,1e-21* paramVals(36) * paramVals(35) * paramVals(32)*pi*paramVals(34)^2 * yCaMKII(:,3)));
Bp1 = floor(max(0,1e-21* paramVals(36) * paramVals(35) * paramVals(32)*pi*paramVals(34)^2 * yCaMKII(:,4)));
Bp2 = floor(max(0,1e-21* paramVals(36) * paramVals(35) * paramVals(32)*pi*paramVals(34)^2 * yCaMKII(:,5)));
Bp3 = floor(max(0,1e-21* paramVals(36) * paramVals(35) * paramVals(32)*pi*paramVals(34)^2 * yCaMKII(:,6)));
Bp4 = floor(max(0,1e-21* paramVals(36) * paramVals(35) * paramVals(32)*pi*paramVals(34)^2 * yCaMKII(:,7)));
Bp5 = floor(max(0,1e-21* paramVals(36) * paramVals(35) * paramVals(32)*pi*paramVals(34)^2 * yCaMKII(:,8)));
Bp6 = floor(max(0,1e-21* paramVals(36) * paramVals(35) * paramVals(32)*pi*paramVals(34)^2 * yCaMKII(:,9)));
Bp7 = floor(max(0,1e-21* paramVals(36) * paramVals(35) * paramVals(32)*pi*paramVals(34)^2 * yCaMKII(:,10)));
Bp8 = floor(max(0,1e-21* paramVals(36) * paramVals(35) * paramVals(32)*pi*paramVals(34)^2 * yCaMKII(:,11)));
Bp9 = floor(max(0,1e-21* paramVals(36) * paramVals(35) * paramVals(32)*pi*paramVals(34)^2 * yCaMKII(:,12)));
Bp10 = floor(max(0,1e-21* paramVals(36) * paramVals(35) * paramVals(32)*pi*paramVals(34)^2 * yCaMKII(:,13)));
Bp11 = floor(max(0,1e-21* paramVals(36) * paramVals(35) * paramVals(32)*pi*paramVals(34)^2 * yCaMKII(:,14)));
Bp12 = floor(max(0,1e-21* paramVals(36) * paramVals(35) * paramVals(32)*pi*paramVals(34)^2 * yCaMKII(:,15)));
Bp13 = floor(max(0,1e-21* paramVals(36) * paramVals(35) * paramVals(32)*pi*paramVals(34)^2 * yCaMKII(:,16)));
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
        U; U(end) + (tSimu(stp+1) - tSimu(stp))*((paramVals(28)-y(end,31))*(6*paramVals(29)*y(end,28) - paramVals(30)*y(end,29) - paramVals(31)*(1-y(end,29)))*y(end,46) + (paramVals(28)-y(end,31))*(paramVals(30)*y(end,29) + paramVals(31)*(1-y(end,29)) - paramVals(29)*y(end,28))*y(end,47))
    ];
end

y = [y, U];

%%
plt_h=4; plt_l=4;
subt = sprintf('Response for impulse CaInit=%0.1f for total CaM=%0.2f', in_CaInit, in_CaM);
figName = sprintf('OutputGB_%s_CaInit%0.1f_CaM%0.2f', paramSetName, in_CaInit, in_CaM);
subplot = @(m,n,p) subtightplot (m, n, p, [0.05 0.035], [0.1 0.01], [0.1 0.01]);

for idx = 1:length(vars)
    if mod(idx,plt_h*plt_l)=1
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