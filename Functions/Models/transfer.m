function wf = transfer(rho, protocParams)
%ZETA Summary of this function goes here
%   Detailed explanation goes here

S_attr = protocParams.S_attr;
sigma = protocParams.noise_lvl;
gamma_pot = protocParams.gamma_pot;
gamma_dep = protocParams.gamma_dep;
tau = protocParams.tau_rho;
n_iter = protocParams.n_iter;
freq = protocParams.frequency; % in ms^(-1)

alpha_pot = protocParams.alpha_pot;
alpha_dep = protocParams.alpha_dep;

G_pot = gamma_pot*alpha_pot;
G_dep = gamma_dep*alpha_dep;

% Time-dependent noise (actual formula)
sigma_eff = sigma.*sqrt((alpha_pot + alpha_dep)./(G_pot + G_dep));
sigma_eff(isnan(sigma_eff)) = 0;

tau_eff = (alpha_pot + alpha_dep > 0).*tau./(G_pot + G_dep);
tau_eff(isnan(tau_eff)) = 0;

rho_eq = G_pot./(G_pot + G_dep);
rho_eq(isnan(rho_eq)) = 0;

wf = 0.5*erfc((S_attr-real(rho))./(sigma_eff.*sqrt(1-exp(-2*1000*n_iter./(freq*tau_eff)))));
wf(isnan(wf)) = 0.5;

end

