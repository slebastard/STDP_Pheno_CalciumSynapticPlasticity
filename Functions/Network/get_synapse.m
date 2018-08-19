function syn = get_synapse( )
%GET_SYNAPSE Summary of this function goes here
%   Detailed explanation goes here
    syn.Cpre = 0.4;
    syn.Cpost = 0.84;
    syn.delay = 5e-3;
    
    syn.tDep = 1;
    syn.gDep = 200;
    
    syn.tPot = 1.08;
    syn.gPot = 120;
    
    syn.tauCa = 3e-2;
    syn.tauRho = 100;
    
    syn.rho0 = 25;
    syn.rhoMax = 200;
    syn.sAttr = 40;
    syn.sigma = 25;
    
end

