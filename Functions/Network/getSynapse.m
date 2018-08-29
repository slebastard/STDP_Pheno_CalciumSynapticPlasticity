function syn = getSynapse( )
%GET_SYNAPSE Summary of this function goes here
%   Detailed explanation goes here
    syn.Cpre = 0.4;
    syn.Cpost = 0.84;
    syn.delay = 5e-3;
    
    syn.tDep = 1.05;
    syn.gDep = 200;
    
    syn.tPot = 1.08;
    syn.gPot = 380;
    
    syn.tauCa = 3e-2;
    syn.tauRho = 100;
    
    syn.rho0 = 25;
    syn.rhoMax = 200;
    syn.sAttr = 100;
    syn.sigma = 25;
    
    syn.tauDamp = 0.1; % from 2e-2 (pure GluR4) to 1e-1 (pure GluR1)
    syn.dampFactor = 0.3;
    
end

