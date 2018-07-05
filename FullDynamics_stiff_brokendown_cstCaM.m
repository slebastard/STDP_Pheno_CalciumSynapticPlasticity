clear all;

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
syms C(t)
syms B0(t) B1(t) B2(t) B3(t) B4(t) B5(t) B6(t) B7(t) B8(t) B9(t) B10(t) B11(t) B12(t) B13(t) chi(t) nu(t) rr(t)
syms Sp(t) Su(t) AMPA_bnd(t) Cb(t) 
syms vPKA_I1(t) vPKA_phos(t) vCaN_I1(t) vCaN_endo(t) vPP1_pase(t) vCK2_exo(t)
syms PP1(t) I1P(t) mu(t) U(t)
syms k10(t)
%Params
syms tauCa CaBas Stot CaM K5 K9 L1 L2 L3 L4 k6 k7 k8 k19 k17 k18 KM k12 k11 km11 I10 PP10 Kdcan ncan
syms kcan0_I1 kcan_I1 kcan0_endo kcan_endo kPP10_pase kPP1_pase Kdpka npka kpka0_I1 kpka_I1 kpka0_phos kpka_phos
syms kCK2_exo kNMDA_bind N g0 g1 g2 M
syms gam_u gam_p zet_u zet_p

% QUANTITIES AND EQUATIONS

% Membrane AMPA dynamics
% ode_n0 = diff(n0) == vCK2_exo*(N-n0-n1) - vCaN_endo*n0 - vPKA_phos*n0*U + vPP1_pase*n1*PP1;
% ode_n1 = iff(n1) == vPKA_phos*n0*U - vPP1_pase*n1*PP1;

% W = g0*n0 + g1*n1;

daes = [
    rr(t) == B1(t) + B2(t) + B3(t) + B4(t) + B5(t) + B6(t) + B7(t) + B8(t) + B9(t) + B10(t) + B11(t) + B12(t) + B13(t);
    B0(t) == Stot - rr(t);

    Sp(t) == B1(t) + 2*(B2(t) + B3(t) + B4(t)) + 3*(B5(t) + B6(t) + B7(t) + B8(t)) + 4*(B9(t) + B10(t) + B11(t)) + 5*B12(t) + 6*B13(t);
    Su(t) == 6*Stot - Sp(t);
    k10(t) == k12*PP1(t)/(KM + (1+zet_p)*Sp(t) + zet_u*Su(t));
    
    C(t) == CaM/(1 + L4/CaBas + L3*L4/(CaBas^2) + L2*L3*L4/(CaBas^3) + L1*L2*L3*L4/(CaBas^4));
    Cb(t) == gam_u*Su(t) + gam_p*Sp(t);

    chi(t) == k7*gam_p + k8*(1 - gam_p - zet_p) + k19*zet_p;
    nu(t) == k10*PP1(t)*(1+zet_p);
    
    AMPA_bnd(t) == 2*(B2(t)+B3(t)+B4(t)) + 6*(B5(t)+B6(t)+B7(t)+B8(t)) + 12*(B9(t)+B10(t)+B11(t)) + 20*B12(t) + 30*B13(t);

    vPKA_I1(t) == kpka0_I1 + kpka_I1/(1 + (Kdpka/C(t))^npka);
    vPKA_phos(t) == kpka0_phos + kpka_phos/(1 + (Kdpka/C(t))^npka);

    vCaN_I1(t) == kcan0_I1 + kcan_I1/(1 + (Kdcan/C(t))^ncan);
    vCaN_endo(t) == kcan0_endo + kcan_endo/(1 + (Kdcan/C(t))^ncan);

    vPP1_pase(t) == kPP10_pase + kPP1_pase/(1 + (Kdcan/C(t))^ncan);

    vCK2_exo(t) == kCK2_exo*Sp(t);
];

odes = [
	diff(B1(t), t) == 6*k6*gam_u^2*B0(t) - 4*k6*gam_u^2*B1(t) - chi(t)*gam_u*B1(t) + nu(t)*(2*(B2(t)+B3(t)+B4(t))-B1(t)) - kNMDA_bind*(M-mu(t))*B1(t);
	diff(B2(t), t) == k6*gam_u^2*B1(t) + chi(t)*gam_u*B1(t) + nu(t)*(3*(B5(t)+B6(t)+B7(t)+B8(t))-2*B2(t)) - 3*k6*gam_u^2*B2(t) - chi(t)*gam_u*B2(t) - 2*kNMDA_bind*(M-mu(t))*B2(t);
	diff(B3(t), t) == 2*k6*gam_u^2*B1(t) + nu(t)*(3*(B5(t)+B6(t)+B7(t)+B8(t))-2*B3(t)) - 3*k6*gam_u^2*B3(t) - chi(t)*gam_u*B3(t) - 2*kNMDA_bind*(M-mu(t))*B3(t);
	diff(B4(t), t) == k6*gam_u^2*B1(t) + nu(t)*(3*(B5(t)+B6(t)+B7(t)+B8(t))-2*B4(t)) - 2*k6*gam_u^2*B4(t) - 2*chi(t)*gam_u*B4(t) - 2*kNMDA_bind*(M-mu(t))*B4(t);
	diff(B5(t), t) == k6*gam_u^2*(B2(t)+B3(t)) + chi(t)*gam_u*B2(t) + nu(t)*(4*(B9(t)+B10(t)+B11(t))-3*B5(t)) - 2*k6*gam_u^2*B5(t) - chi(t)*gam_u*B5(t) - 3*kNMDA_bind*(M-mu(t))*B5(t);
	diff(B6(t), t) == k6*gam_u^2*(B2(t)+B3(t)) + 2*chi(t)*gam_u*B4(t) + nu(t)*(4*(B9(t)+B10(t)+B11(t))-3*B6(t)) - k6*gam_u^2*B6(t) - 2*chi(t)*gam_u*B6(t) - 3*kNMDA_bind*(M-mu(t))*B6(t);
	diff(B7(t), t) == k6*gam_u^2*(B2(t)+2*B4(t)) + chi(t)*gam_u*B3(t) + nu(t)*(4*(B9(t)+B10(t)+B11(t))-3*B7(t)) - k6*gam_u^2*B7(t) - 2*chi(t)*gam_u*B7(t) - 3*kNMDA_bind*(M-mu(t))*B7(t);
	diff(B8(t), t) == k6*gam_u^2*B3(t) + nu(t)*(4*(B9(t)+B10(t)+B11(t))-3*B8(t)) - 3*chi(t)*gam_u*B8(t) - 3*kNMDA_bind*(M-mu(t))*B8(t);
	diff(B9(t), t) == k6*gam_u^2*B5(t) + chi(t)*gam_u*(B6(t)+B7(t)) + nu(t)*(5*B12(t)-4*B9(t)) - k6*gam_u^2*B9(t) - chi(t)*gam_u*B9(t) - 4*kNMDA_bind*(M-mu(t))*B9(t);
	diff(B10(t), t) == k6*gam_u^2*(B5(t)+B6(t)) + chi(t)*gam_u*(B7(t)+B8(t)) + nu(t)*(5*B12(t)-4*B10(t)) - 2*chi(t)*gam_u*B10(t) - 4*kNMDA_bind*(M-mu(t))*B10(t);
	diff(B11(t), t) == k6*gam_u^2*B7(t) + chi(t)*gam_u*B6(t) + nu(t)*(5*B12(t)-4*B11(t)) - 2*chi(t)*gam_u*B11(t) - 4*kNMDA_bind*(M-mu(t))*B11(t);
	diff(B12(t), t) == k6*gam_u^2*B9(t) + chi(t)*gam_u*(B9(t)+2*B10(t)+2*B11(t)) + nu(t)*(6*B13(t)-5*B12(t)) - chi(t)*gam_u*B12(t) - 5*kNMDA_bind*(M-mu(t))*B11(t);
	diff(B13(t), t) == chi(t)*gam_u*B12(t) - nu(t)*6*B13(t) - 6*kNMDA_bind*(M-mu(t))*B13(t);
	diff(PP1(t), t) == -k11*I1P(t)*PP1(t) + km11*(PP10 - PP1(t));
	diff(I1P(t), t) == -k11*I1P(t)*PP1(t) + km11*(PP10 - PP1(t)) + vPKA_I1(t)*(I10-I1P(t)) - vCaN_I1(t)*I1P(t);
	diff(mu(t), t) == kNMDA_bind*(M-mu(t))*Sp(t);
	diff(U(t), t) == kNMDA_bind*(M-mu(t))*AMPA_bnd(t);
       % TOTAL 18 ODEs
];

eqs = [odes;daes];

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
    gam_u; gam_p; zet_u; zet_p
];

paramVals = [
	0.012; 0.5;
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
	1000; 0.0010; 0.0017; 0.0024; 100
];

C0 = getC(paramVals(4),paramVals(2),paramVals(7),paramVals(8),paramVals(9),paramVals(10));
[p0,i0] = getP0(C0, paramVals(23:26), paramVals(31:34), paramVals(19:22));

[gu,gp,zu,zp] = getRates(C0, 199.8, 0, paramVals(5), paramVals(6), paramVals(15), paramVals(16), paramVals(18), paramVals(17), p0);

paramVals = [
    paramVals;
    gu; gp; zu; zp
];

vars = [
	C(t);
	B0(t); B1(t); B2(t); B3(t); B4(t); B5(t); B6(t); B7(t); B8(t);
    B9(t); B10(t); B11(t); B12(t); B13(t); chi(t); nu(t); rr(t);
	Sp(t); Su(t); AMPA_bnd(t); Cb(t);
	vPKA_I1(t); vPKA_phos(t); vCaN_I1(t); vCaN_endo(t); vPP1_pase(t); vCK2_exo(t);
	PP1(t); I1P(t);
    mu(t); U(t);
    k10(t)
];

y0est = [
	C0;
	33.3; 0; 0; 0; 0; 0; 0; 0; 0;
    0; 0; 0; 0; 0; 0; 0; 0;
	0; 199.8; 0; 0;
	0; 0; 0; 0; 0; 0;
	p0; i0;
    0; 0;
    0
];

y0fix = [
	0;
	0; 0; 0; 0; 0; 0; 1; 1; 1;
    1; 1; 1; 1; 1; 0; 0; 0;
	0; 0; 0; 0;
	0; 0; 0; 0; 0; 0;
	1; 1;
    0; 0;
    0
];

[mass,F] = massMatrixForm(eqs, vars);
mass = odeFunction(mass, vars);
F = odeFunction(F, vars, params);
f = @(t, y) F(t, y, paramVals);

t0=0; tfinal=10; step=0.001;
nsteps = (tfinal - t0)/step;

opt_ode = odeset('Mass', mass,...
'AbsTol',1e-6,...
'RelTol',1e-3);

opt_fsolve = optimoptions('fsolve','Display','off');

implicitDAE = @(t,y,yp) mass(t,y)*yp - f(t,y);
[y0, yp0] = decic(implicitDAE, t0, y0est, y0fix, zeros(33,1), [], opt_ode);

yi = y0';
t = []; y = [];
ratesHist=[];
gamuRoots = [];

eps = (paramVals(15)*paramVals(6))/(paramVals(16)*paramVals(5));

for nstep=1:nsteps
    lbd = paramVals(6)*(1 - paramVals(15)/paramVals(16))/yi(end,1) + 1 - eps;
    
%     a = vpasolve(...
%         paramVals(16)*paramVals(5)*gam_u*((paramVals(17)+2*yi(end,19)+yi(end,20))*yi(end,1)*(lbd*gam_u+eps)...
%         - yi(end,20)*(paramVals(5)+yi(end,1))*gam_u*(lbd*gam_u+eps) - yi(end,19)*(paramVals(6)+yi(end,1))*gam_u)...
%         - paramVals(18)*yi(end,29)*yi(end,1)*(yi(end,1)*(lbd*gam_u+eps)-(paramVals(6)+yi(end,1))*gam_u) == 0 ...
%         , gam_u, paramVals(44)...
%     );

    fun =  @(y) root_gu(y, yi, paramVals, eps, lbd);
    a = fsolve(fun,paramVals(44), opt_fsolve);

    gamuRoots = [gamuRoots; a];
    
    % appl_gu = a(a >=0 & a <= 1 & isreal(a));
    appl_gu = min(max(0, a(imag(a) == 0)), 1);
    dist_gu = abs(appl_gu - paramVals(44));
    [dist_gu_sorted, dist_gu_order] = sort(dist_gu);
    sorted_gu = appl_gu(dist_gu_order,:);

    pvals = paramVals(44:47);
    try
        paramVals(44) = sorted_gu(1);
        if paramVals(44) == 0 || paramVals(44) == 1
            warning('Warning: gam_u saturated')
        end
        
        paramVals(45) = min(max(0, paramVals(44)/(lbd*paramVals(44) + eps)), 1);
        if paramVals(45) == 0 || paramVals(45) == 1
            warning('Warning: gam_p saturated')
        end
        
        paramVals(46) = min(max(0, 1 - paramVals(44)*(paramVals(5)+yi(end,1))/yi(end,1)), 1);
        if paramVals(46) == 0 || paramVals(46) == 1
            warning('Warning: zet_u saturated')
        end
        
        paramVals(47) = min(max(0, 1 - paramVals(45)*(paramVals(6)+yi(end,1))/yi(end,1)), 1);
        if paramVals(47) == 0 || paramVals(47) == 1
            warning('Warning: zet_p saturated')
        end
    catch
        warning('Warning: solution to system not found. Using previous values of gamma and zeta, continuing to integrate equation');
        paramVals(44:47) = pvals(1:4);
    end
    f = @(t, y) F(t, y, paramVals);

    tstep = t0 + (nstep-1)*step;
    [ti,yi] = ode15s(f, [tstep, tstep+step], yi(end,:), opt_ode);
    t = [t;ti];
    y = [y;yi];
    
    ratesHist = [ratesHist; repmat(paramVals(44:47).',length(ti),1)];
end

%%
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

idx0=idx;
for idx=1:4
    subplot(plt_h,plt_l,1+mod(idx0+idx-1,plt_h*plt_l))
    plot(t(:,1),ratesHist(:,idx), 'x')
    title(char(params(43+idx)))
end

function K = getC(CaM,CaBas, L1, L2, L3, L4)
    K = double(CaM/(1 + L4/CaBas + L3*L4/(CaBas^2) + L2*L3*L4/(CaBas^3) + L1*L2*L3*L4/(CaBas^4)));
end

function [p0,i0] = getP0(C, kin_CaN, kin_PKA, par_PP1)
    syms loc_p0 loc_i0

    vCaN = kin_CaN(3) + kin_CaN(4)/(1 + (kin_CaN(1)/C)^kin_CaN(2));
    vPKA = kin_PKA(3) + kin_PKA(4)/(1 + (kin_PKA(1)/C)^kin_PKA(2));

    i0 = min(par_PP1(3)*vPKA/vCaN,par_PP1(3));
    p0 = par_PP1(4)/(1+i0*par_PP1(1)/par_PP1(2));
end


function [gu, gp, zu, zp] = getRates(C, Su, Sp, K5, K9, k17, k18, k12, KM, PP1)
    syms loc_gu
    
    eps = (k17*K9)/(k18*K5);
    lbd = K9*(1 - k17/k18)/C + 1 - eps;
    
    a = vpasolve(...
        k18*K5*loc_gu*((KM+2*Sp+Su)*C*(lbd*loc_gu+eps)...
        - Su*(K5+C)*loc_gu*(lbd*loc_gu+eps) - Sp*(K9+C)*loc_gu)...
        - k12*PP1*C*(C*(lbd*loc_gu+eps)-(K9+C)*loc_gu) == 0 ...
        ,loc_gu,[0 1]...
    );

    gu = min(max(0.0, double(a)), 1.0);
    gp = min(max(0.0, gu/(lbd*gu + eps)), 1.0);
    zu = min(max(0.0, 1 - gu*(K5+C)/C), 1.0);
    zp = min(max(0.0, 1 - gp*(K9+C)/C), 1.0);

end

function r = root_gu(gam_u, yi, paramVals, eps, lbd)
        r(1) = paramVals(16)*paramVals(5)*gam_u*((paramVals(17)+2*yi(end,19)+yi(end,20))*yi(end,1)*(lbd*gam_u+eps)...
        - yi(end,20)*(paramVals(5)+yi(end,1))*gam_u*(lbd*gam_u+eps) - yi(end,19)*(paramVals(6)+yi(end,1))*gam_u)...
        - paramVals(18)*yi(end,29)*yi(end,1)*(yi(end,1)*(lbd*gam_u+eps)-(paramVals(6)+yi(end,1))*gam_u);
end