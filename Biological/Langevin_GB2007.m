function [Sp_alpha, Sp_infty] = Langevin_GB2007(t0, tfinal, dt, Ca_stim, tstim)

% Naming conventions
% k --- reaction constant
% r --- reaction rate ------------------------------- r_{11} = k_{11}*[C_p]*[PP1_{free}]
% K --- dissociation constant ----------------------- K_5 = k_{-5} / k_{5}
% KM -- Michaelis constant for catalytic reaction --- KM_Dpu = (k_{-11} + k_{12})/k_{11}

% All time units in s
% All concentrations in microMolars = 1mumol/L

% Variables for cytosolic CaMKII pathways
% ID    Name        Description
% 1     Ca          Free calcium concentration
% 2     B1          Concentration of olygomer 1
% 3     B2          Concentration of olygomer 2
% 4     B3          Concentration of olygomer 3
% 5     B4          Concentration of olygomer 4
% 6     B5          Concentration of olygomer 5
% 7     B6          Concentration of olygomer 6
% 8     B7          Concentration of olygomer 7
% 9     B8          Concentration of olygomer 8
% 10    B9          Concentration of olygomer 9
% 11    B10         Concentration of olygomer 10
% 12    B11         Concentration of olygomer 11
% 13    B12         Concentration of olygomer 12 
% 14    B13         Concentration of olygomer 13
% 15    PP1         Concentration of active PP1
% 16    I1P         Concentration of phosphorylated PP inhibitor 1

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

% --> 37 t_stim     Duration of stimulation at [Ca2+] = 1muM

% INPUTS (default below)
% t0=0; tfinal=100; dt=0.01; in_CaInit = 0.1;

Ca_alpha = 0.219;
paramSetName = 'Graupner';
in_CaM = 0.1;
%Ca_stim = 1;

showPlots = true;

% GET PARAMETERS AND INITIAL VALUES
paramVals = getParams(paramSetName, in_CaM, tstim);

C0 = getC(Ca_stim,paramVals(4),paramVals(7),paramVals(8),paramVals(9),paramVals(10));
[p0,i0] = getP0(C0, paramVals(20:23), paramVals(24:27), paramVals(16:19));

y0 = [
    Ca_stim;
    0; 0; 0; 0; 0; 0; 0; 0;
    0; 0; 0; 0; 0
    p0; i0
];

% Integrate equations
[eqs_det1, eqs_det2, eqs_stoch, ~] = getEqs(paramVals);
 
opts = sdeset('NonNegative',1:16,...
              'RandSeed',1,...
              'SDEType','Ito');

% y1 = sde_euler(eqs_det1, eqs_stoch, t0:dt:tstim, y0', opts);
% y2 = sde_euler(eqs_det2, eqs_stoch, tstim+dt:dt:tfinal, y1(end,:), opts);
y1 = ode45(eqs_det1, [t0 tstim], y0');
y2 = ode45(eqs_det2, [tstim tfinal], y1.y(:,end));
y = [y1.y';y2.y'];

% t = t0:dt:tfinal;
t = [y1.x,y2.x];

% Finding Sp when calcium goes beyond given threshold
t_alpha = tstim + paramVals(1)*log(Ca_stim/(Ca_alpha-paramVals(2)));
Sp = y(:,2) + 2*(y(:,3) + y(:,4) + y(:,5)) + 3*(y(:,6) + y(:,7) + y(:,8) + y(:,9)) + 4*(y(:,10) + y(:,11) + y(:,12)) + 5*y(:,13) + 6*y(:,14);
S0 = y(:,2) + y(:,3) + y(:,4) + y(:,5) + y(:,6) + y(:,7) + y(:,8) + y(:,9) + y(:,10) + y(:,11) + y(:,12) + y(:,13) + y(:,14);
B0 = paramVals(3) - S0;
id_alpha = floor(t_alpha/dt);
% Sp_alpha = Sp(id_alpha,1);
Sp_infty = Sp(end,1);

% [Optionnal, disabled by default] Plot and save variables
if showPlots
    figure(1);
    plot(t,y(:,1),'+');
    title('Ca2+ as a function of time');
    
    figure(2);
    plot(t,Sp,'+');
    title('Phosphorylated CaMKII as a function of time');
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

function paramVals = getParams(author, CaM, t_stim)
    if strcmp(author,'Graupner')
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
            0.5; 0.2; 0.1; 6.02e23; 0.05;
            
            t_stim
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
            0.5; 0.2; 0.1; 6.02e23; 0.05;
            
            t_stim
        ];
    end
end

function [eqs_det1, eqs_det2, eqs_stoch, aux] = getEqs(params)  
    C = @(t, vars) (params(4)/(1 + params(10)/vars(1) + params(9)*params(10)/(vars(1)^2) + params(8)*params(9)*params(10)/(vars(1)^3) + params(7)*params(8)*params(9)*params(10)/(vars(1)^4)));
    S0 = @(t, vars) (vars(2) + vars(3) + vars(4) + vars(5) + vars(6) + vars(7) + vars(8) + vars(9) + vars(10) + vars(11) + vars(12) + vars(13) + vars(14));
    B0 = @(t, vars) (params(3) - S0(t, vars));
    gam_u = @(t, vars) (C(t, vars)/(params(5) + C(t, vars)));
    gam_p = @(t, vars) (C(t, vars)/(params(6) + C(t, vars)));
    Sp = @(t, vars) (vars(2) + 2*(vars(3) + vars(4) + vars(5)) + 3*(vars(6) + vars(7) + vars(8) + vars(9)) + 4*(vars(10) + vars(11) + vars(12)) + 5*vars(13) + 6*vars(14));
    Su = @(t, vars) (6*params(3) - Sp(t, vars));
    Cb = @(t, vars) (gam_u(t, vars)*Su(t, vars) + gam_p(t, vars)*Sp(t, vars));
    k10 = @(t, vars) (params(15)/(params(14) + Sp(t, vars)));
    chi = @(t, vars) (params(12)*gam_p(t, vars) + params(13)*(1 - gam_p(t, vars)));
    nu = @(t, vars) (k10(t, vars)*vars(15));
    AMPA_bnd = @(t, vars) (2*(vars(3) + vars(4) + vars(5)) ...
    + 6*(vars(6) + vars(7) + vars(8) + vars(9)) ...
    + 12*(vars(10) + vars(11) + vars(12)) ...
    + 20*vars(13) + 30*vars(14));
    vPKA_I1 = @(t, vars) (params(26) + params(27)/(1 + (params(24)/(C(t, vars)))^params(25)));
    vCaN_I1 = @(t, vars) (params(22) + params(23)/(1 + (params(20)/(C(t, vars)))^params(21)));
    
    eqs_det1 = @(t, vars) [
        0;
        6*params(11)*gam_u(t, vars)^2*B0(t, vars) - 4*params(11)*gam_u(t, vars)^2*vars(2) - chi(t, vars)*gam_u(t, vars)*vars(2) + nu(t, vars)*(2*(vars(3)+vars(4)+vars(5))-vars(2));
        params(11)*gam_u(t, vars)^2*vars(2) + chi(t, vars)*gam_u(t, vars)*vars(2) + nu(t, vars)*(3*(vars(6)+vars(7)+vars(8)+vars(9))-2*vars(3)) - 3*params(11)*gam_u(t, vars)^2*vars(3) - chi(t, vars)*gam_u(t, vars)*vars(3);
        2*params(11)*gam_u(t, vars)^2*vars(2) + nu(t, vars)*(3*(vars(6)+vars(7)+vars(8)+vars(9))-2*vars(4)) - 3*params(11)*gam_u(t, vars)^2*vars(4) - chi(t, vars)*gam_u(t, vars)*vars(4);
        params(11)*gam_u(t, vars)^2*vars(2) + nu(t, vars)*(3*(vars(6)+vars(7)+vars(8)+vars(9))-2*vars(5)) - 2*params(11)*gam_u(t, vars)^2*vars(5) - 2*chi(t, vars)*gam_u(t, vars)*vars(5);
        params(11)*gam_u(t, vars)^2*(vars(3)+vars(4)) + chi(t, vars)*gam_u(t, vars)*vars(3) + nu(t, vars)*(4*(vars(10)+vars(11)+vars(12))-3*vars(6)) - 2*params(11)*gam_u(t, vars)^2*vars(6) - chi(t, vars)*gam_u(t, vars)*vars(6);
        params(11)*gam_u(t, vars)^2*(vars(3)+vars(4)) + 2*chi(t, vars)*gam_u(t, vars)*vars(5) + nu(t, vars)*(4*(vars(10)+vars(11)+vars(12))-3*vars(7)) - params(11)*gam_u(t, vars)^2*vars(7) - 2*chi(t, vars)*gam_u(t, vars)*vars(7);
        params(11)*gam_u(t, vars)^2*(vars(3)+2*vars(5)) + chi(t, vars)*gam_u(t, vars)*vars(4) + nu(t, vars)*(4*(vars(10)+vars(11)+vars(12))-3*vars(8)) - params(11)*gam_u(t, vars)^2*vars(8) - 2*chi(t, vars)*gam_u(t, vars)*vars(8);
        params(11)*gam_u(t, vars)^2*vars(4) + nu(t, vars)*(4*(vars(10)+vars(11)+vars(12))-3*vars(9)) - 3*chi(t, vars)*gam_u(t, vars)*vars(9);
        params(11)*gam_u(t, vars)^2*vars(6) + chi(t, vars)*gam_u(t, vars)*(vars(7)+vars(8)) + nu(t, vars)*(5*vars(13)-4*vars(10)) - params(11)*gam_u(t, vars)^2*vars(10) - chi(t, vars)*gam_u(t, vars)*vars(10);
        params(11)*gam_u(t, vars)^2*(vars(6)+vars(7)) + chi(t, vars)*gam_u(t, vars)*(vars(8)+vars(9)) + nu(t, vars)*(5*vars(13)-4*vars(11)) - 2*chi(t, vars)*gam_u(t, vars)*vars(11);
        params(11)*gam_u(t, vars)^2*vars(8) + chi(t, vars)*gam_u(t, vars)*vars(7) + nu(t, vars)*(5*vars(13)-4*vars(12)) - 2*chi(t, vars)*gam_u(t, vars)*vars(12);
        params(11)*gam_u(t, vars)^2*vars(10) + chi(t, vars)*gam_u(t, vars)*(vars(10)+2*vars(11)+2*vars(12)) + nu(t, vars)*(6*vars(14)-5*vars(13)) - chi(t, vars)*gam_u(t, vars)*vars(13);
        chi(t, vars)*gam_u(t, vars)*vars(13) - nu(t, vars)*6*vars(14);
        -params(16)*vars(16)*vars(15) + params(17)*(params(19) - vars(15));
        -params(16)*vars(16)*vars(15) + params(17)*(params(19) - vars(15)) + vPKA_I1(t, vars)*(params(18)-vars(16)) - vCaN_I1(t, vars)*vars(16);
    ];

    eqs_det2 = @(t, vars) [
        -(vars(1) - params(2))/(params(1));
        6*params(11)*gam_u(t, vars)^2*B0(t, vars) - 4*params(11)*gam_u(t, vars)^2*vars(2) - chi(t, vars)*gam_u(t, vars)*vars(2) + nu(t, vars)*(2*(vars(3)+vars(4)+vars(5))-vars(2));
        params(11)*gam_u(t, vars)^2*vars(2) + chi(t, vars)*gam_u(t, vars)*vars(2) + nu(t, vars)*(3*(vars(6)+vars(7)+vars(8)+vars(9))-2*vars(3)) - 3*params(11)*gam_u(t, vars)^2*vars(3) - chi(t, vars)*gam_u(t, vars)*vars(3);
        2*params(11)*gam_u(t, vars)^2*vars(2) + nu(t, vars)*(3*(vars(6)+vars(7)+vars(8)+vars(9))-2*vars(4)) - 3*params(11)*gam_u(t, vars)^2*vars(4) - chi(t, vars)*gam_u(t, vars)*vars(4);
        params(11)*gam_u(t, vars)^2*vars(2) + nu(t, vars)*(3*(vars(6)+vars(7)+vars(8)+vars(9))-2*vars(5)) - 2*params(11)*gam_u(t, vars)^2*vars(5) - 2*chi(t, vars)*gam_u(t, vars)*vars(5);
        params(11)*gam_u(t, vars)^2*(vars(3)+vars(4)) + chi(t, vars)*gam_u(t, vars)*vars(3) + nu(t, vars)*(4*(vars(10)+vars(11)+vars(12))-3*vars(6)) - 2*params(11)*gam_u(t, vars)^2*vars(6) - chi(t, vars)*gam_u(t, vars)*vars(6);
        params(11)*gam_u(t, vars)^2*(vars(3)+vars(4)) + 2*chi(t, vars)*gam_u(t, vars)*vars(5) + nu(t, vars)*(4*(vars(10)+vars(11)+vars(12))-3*vars(7)) - params(11)*gam_u(t, vars)^2*vars(7) - 2*chi(t, vars)*gam_u(t, vars)*vars(7);
        params(11)*gam_u(t, vars)^2*(vars(3)+2*vars(5)) + chi(t, vars)*gam_u(t, vars)*vars(4) + nu(t, vars)*(4*(vars(10)+vars(11)+vars(12))-3*vars(8)) - params(11)*gam_u(t, vars)^2*vars(8) - 2*chi(t, vars)*gam_u(t, vars)*vars(8);
        params(11)*gam_u(t, vars)^2*vars(4) + nu(t, vars)*(4*(vars(10)+vars(11)+vars(12))-3*vars(9)) - 3*chi(t, vars)*gam_u(t, vars)*vars(9);
        params(11)*gam_u(t, vars)^2*vars(6) + chi(t, vars)*gam_u(t, vars)*(vars(7)+vars(8)) + nu(t, vars)*(5*vars(13)-4*vars(10)) - params(11)*gam_u(t, vars)^2*vars(10) - chi(t, vars)*gam_u(t, vars)*vars(10);
        params(11)*gam_u(t, vars)^2*(vars(6)+vars(7)) + chi(t, vars)*gam_u(t, vars)*(vars(8)+vars(9)) + nu(t, vars)*(5*vars(13)-4*vars(11)) - 2*chi(t, vars)*gam_u(t, vars)*vars(11);
        params(11)*gam_u(t, vars)^2*vars(8) + chi(t, vars)*gam_u(t, vars)*vars(7) + nu(t, vars)*(5*vars(13)-4*vars(12)) - 2*chi(t, vars)*gam_u(t, vars)*vars(12);
        params(11)*gam_u(t, vars)^2*vars(10) + chi(t, vars)*gam_u(t, vars)*(vars(10)+2*vars(11)+2*vars(12)) + nu(t, vars)*(6*vars(14)-5*vars(13)) - chi(t, vars)*gam_u(t, vars)*vars(13);
        chi(t, vars)*gam_u(t, vars)*vars(13) - nu(t, vars)*6*vars(14);
        -params(16)*vars(16)*vars(15) + params(17)*(params(19) - vars(15));
        -params(16)*vars(16)*vars(15) + params(17)*(params(19) - vars(15)) + vPKA_I1(t, vars)*(params(18)-vars(16)) - vCaN_I1(t, vars)*vars(16);
    ];
    
    lvl=0.001;
    eqs_stoch = @(t, vars) [
%         0;
%         lvl*(1/sqrt(params(35)*params(32)*pi*params(33)^2))*sqrt(6*params(11)*gam_u(t, vars)^2*S0(t, vars) + 4*params(11)*gam_u(t, vars)^2*vars(2) + chi(t, vars)*gam_u(t, vars)*vars(2) + nu(t, vars)*(2*(vars(3)+vars(4)+vars(5))+vars(2)));
%         lvl*(1/sqrt(params(35)*params(32)*pi*params(33)^2))*sqrt(params(11)*gam_u(t, vars)^2*vars(2) + chi(t, vars)*gam_u(t, vars)*vars(2) + nu(t, vars)*(3*(vars(6)+vars(7)+vars(8)+vars(9))+2*vars(3)) + 3*params(11)*gam_u(t, vars)^2*vars(3) + chi(t, vars)*gam_u(t, vars)*vars(3));
%         lvl*(1/sqrt(params(35)*params(32)*pi*params(33)^2))*sqrt(2*params(11)*gam_u(t, vars)^2*vars(2) + nu(t, vars)*(3*(vars(6)+vars(7)+vars(8)+vars(9))+2*vars(4)) + 3*params(11)*gam_u(t, vars)^2*vars(4) + chi(t, vars)*gam_u(t, vars)*vars(4));
%         lvl*(1/sqrt(params(35)*params(32)*pi*params(33)^2))*sqrt(params(11)*gam_u(t, vars)^2*vars(2) + nu(t, vars)*(3*(vars(6)+vars(7)+vars(8)+vars(9))+2*vars(5)) + 2*params(11)*gam_u(t, vars)^2*vars(5) + 2*chi(t, vars)*gam_u(t, vars)*vars(5));
%         lvl*(1/sqrt(params(35)*params(32)*pi*params(33)^2))*sqrt(params(11)*gam_u(t, vars)^2*(vars(3)+vars(4)) + chi(t, vars)*gam_u(t, vars)*vars(3) + nu(t, vars)*(4*(vars(10)+vars(11)+vars(12))+3*vars(6)) + 2*params(11)*gam_u(t, vars)^2*vars(6) + chi(t, vars)*gam_u(t, vars)*vars(6));
%         lvl*(1/sqrt(params(35)*params(32)*pi*params(33)^2))*sqrt(params(11)*gam_u(t, vars)^2*(vars(3)+vars(4)) + 2*chi(t, vars)*gam_u(t, vars)*vars(5) + nu(t, vars)*(4*(vars(10)+vars(11)+vars(12))+3*vars(7)) + params(11)*gam_u(t, vars)^2*vars(7) + 2*chi(t, vars)*gam_u(t, vars)*vars(7));
%         lvl*(1/sqrt(params(35)*params(32)*pi*params(33)^2))*sqrt(params(11)*gam_u(t, vars)^2*(vars(3)+2*vars(5)) + chi(t, vars)*gam_u(t, vars)*vars(4) + nu(t, vars)*(4*(vars(10)+vars(11)+vars(12))+3*vars(8)) + params(11)*gam_u(t, vars)^2*vars(8) + 2*chi(t, vars)*gam_u(t, vars)*vars(8));
%         lvl*(1/sqrt(params(35)*params(32)*pi*params(33)^2))*sqrt(params(11)*gam_u(t, vars)^2*vars(4) + nu(t, vars)*(4*(vars(10)+vars(11)+vars(12))+3*vars(9)) + 3*chi(t, vars)*gam_u(t, vars)*vars(9));
%         lvl*(1/sqrt(params(35)*params(32)*pi*params(33)^2))*sqrt(params(11)*gam_u(t, vars)^2*vars(6) + chi(t, vars)*gam_u(t, vars)*(vars(7)+vars(8)) + nu(t, vars)*(5*vars(13)+4*vars(10)) + params(11)*gam_u(t, vars)^2*vars(10) + chi(t, vars)*gam_u(t, vars)*vars(10));
%         lvl*(1/sqrt(params(35)*params(32)*pi*params(33)^2))*sqrt(params(11)*gam_u(t, vars)^2*(vars(6)+vars(7)) + chi(t, vars)*gam_u(t, vars)*(vars(8)+vars(9)) + nu(t, vars)*(5*vars(13)+4*vars(11)) + 2*chi(t, vars)*gam_u(t, vars)*vars(11));
%         lvl*(1/sqrt(params(35)*params(32)*pi*params(33)^2))*sqrt(params(11)*gam_u(t, vars)^2*vars(8) + chi(t, vars)*gam_u(t, vars)*vars(7) + nu(t, vars)*(5*vars(13)+4*vars(12)) + 2*chi(t, vars)*gam_u(t, vars)*vars(12));
%         lvl*(1/sqrt(params(35)*params(32)*pi*params(33)^2))*sqrt(params(11)*gam_u(t, vars)^2*vars(10) + chi(t, vars)*gam_u(t, vars)*(vars(10)+2*vars(11)+2*vars(12)) + nu(t, vars)*(6*vars(14)+5*vars(13)) + chi(t, vars)*gam_u(t, vars)*vars(13));
%         lvl*(1/sqrt(params(35)*params(32)*pi*params(33)^2))*sqrt(chi(t, vars)*gam_u(t, vars)*vars(13) + nu(t, vars)*6*vars(14));
%         lvl*(1/sqrt(params(35)*params(32)*pi*params(33)^2))*sqrt(+params(16)*vars(16)*vars(15) + params(17)*(params(19) + vars(15)));
%         lvl*(1/sqrt(params(35)*params(32)*pi*params(33)^2))*sqrt(+params(16)*vars(16)*vars(15) + params(17)*(params(19) + vars(15)) + vPKA_I1(t, vars)*(params(18)+vars(16)) + vCaN_I1(t, vars)*vars(16));
        0;
        0;
        0;
        0;
        0;
        0;
        0;
        0;
        0;
        0;
        0;
        0;
        0;
        0;
        0;
        0;
    ];

    aux = {C; S0; B0; k10; gam_u; gam_p; Sp; Su; Cb; chi; nu; AMPA_bnd; vPKA_I1; vCaN_I1};
end

end