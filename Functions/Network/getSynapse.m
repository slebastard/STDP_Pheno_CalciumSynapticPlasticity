function syn = getSynapse( )
%GET_SYNAPSE Summary of this function goes here
%   Detailed explanation goes here
    syn.C_pre = 0.4;
    syn.C_post = 0.84;
    syn.delay_pre = -15e-3;
    
    syn.theta_dep = 1;
    syn.gamma_dep = 200;
    
    syn.theta_pot = 1.2;
    syn.gamma_pot = 70;
    
    syn.tau_Ca = 80e-3;
    syn.tau_rho = 100;
    
    syn.rho_0 = 25;
    syn.rho_max = 200;
    syn.S_attr = 40;
    syn.noise_lvl = 25;
    
    syn.tau_x = 0.1; % from 2e-2 (pure GluR4) to 1e-1 (pure GluR1)
    syn.dampFactor = 0.3;
    syn.w_0 = transfer_ind(syn.rho_0, syn);
    
end

