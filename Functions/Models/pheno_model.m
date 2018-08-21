function [rho_hist, w_hist, c_hist] = pheno_model( pre_spikes_hist, post_spikes_hist, params, simu)
%NAIVE_MODEL Simulates the behavior of a synapse whose behavior follows the
%naive Calcium_based dynamics
%   Detailed explanation goes here

%% Default parameter values + unpacking params
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch nargin
    case 0
        pre_spikes_hist = [];
        post_spikes_hist = [];
        params = default_params();
        simu = default_simu();
    case 1
        post_spikes_hist = [];
        params = default_params();
        simu = default_simu();
    case 2
        params = default_params();
        simu = default_simu();
    case 3
        simu = default_simu();
    case 4
    otherwise
        error('4 inputs max are accepted')
end

%%%%%%%%%%%%%%%%%%%%
% Unpacking params %
%%%%%%%%%%%%%%%%%%%%

T = simu.T;
step = simu.int_step;
scheme = simu.int_scheme;
n_steps = T / step;

rho_0 = params.rho_0;
rho_max = params.rho_max;
C_pre = params.C_pre;
C_post = params.C_post;
tau_Ca = params.tau_Ca;
delay_pre = params.delay_pre;
theta_dep = params.theta_dep;
gamma_dep = params.gamma_dep;
theta_pot = params.theta_pot;
gamma_pot = params.gamma_pot;
tau_rho = params.tau_rho;
sigma = params.noise_lvl;
S_attr = params.S_attr;
tau_w = params.tau_w;
theta_act = params.theta_act;
w_0 = params.w_0;

prot = params;
prot.n_iter = 1;

prot.alpha_pot = 0;
prot.alpha_dep = 0;
prot.frequency = 1000/step;

%% Building events list based on calcium hypothesis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

len_pre = length(pre_spikes_hist);
len_post = length(post_spikes_hist);

evts = [];

if len_pre > 0
    for pre_evtID=1:len_pre
        evts = cat(1, evts, [pre_spikes_hist(pre_evtID) + delay_pre, C_pre]);
    end
end

if len_post > 0
    for post_evtID=1:len_post
        evts = cat(1, evts, [post_spikes_hist(post_evtID),C_post]);
    end
end

evts = sortrows(evts, 1);

%% Simulating process
%%%%%%%%%%%%%%%%%%%%%

if strcmp(scheme, 'euler_expl')
    % Initiate simulation %
    %%%%%%%%%%%%%%%%%%%%%%%
    rho = rho_0;
    w = w_0;
    t = 0;
    c = 0;
    
    rho_hist = rho;
    w_hist = w;
    c_hist = 0;

    % Check whether simulation is trivial %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if isempty(evts)
        rho_hist = rho_0 * ones(n_steps);
        w_f = transfer(rho, prot);
        w_hist = (w_0-w_f)*exp(-step*(0:n_steps-1)/tau_w);
    else
        % Scheme propagation %
        %%%%%%%%%%%%%%%%%%%%%%
        while t < T
            c_hist = [c_hist, c];
            evt = find(and(t - double(evts(:,1)) < step, t >= double(evts(:,1))));
            if evt
                C_bump = sum(evts(evt,2));
                c = c + C_bump; 
            end
            
            rho = rho + step/tau_rho * (gamma_pot*(rho_max-rho)*(c > theta_pot) - gamma_dep*rho*(c > theta_dep)) * (c > theta_act);
            
            w_f = transfer(rho, prot);
            w = w + step/tau_w * (w_f-w);
            
            rho_hist = [rho_hist, rho];
            w_hist = [w_hist; w];
            c = c * exp(-step/tau_Ca);
            t = t + step;
            prot.alpha_pot = (1+step/t)*(prot.alpha_pot+(step/t)*(c > theta_pot));
            prot.alpha_dep = (1+step/t)*(prot.alpha_dep+(step/t)*(c > theta_dep));
            prot.frequency = 1000/t;
        end
    end     
        
else
	error('This integration scheme is not currently handled')
end

end
