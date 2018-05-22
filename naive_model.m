function [rho_hist, c_hist] = naive_model( pre_spikes_hist, post_spikes_hist, params, int_scheme, int_step, noise_lvl)
%NAIVE_MODEL Simulates the behavior of a synapse whose behavior follows the
%naive Calcium_based dynamics
%   Detailed explanation goes here

%% Default parameter values + unpacking params
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

def_params = [...
    1000 ...        % T         total simu time     (ms)
    .3 ...          % rho_0     init syn strength
    1 ...           % C_pre
    2 ...           % C_post
    20 ...          % tau_Ca
    3 ...           % delay_pre
    1 ...           % theta_dep
    200 ...         % gamma_dep
    1.3 ...         % theta_pot
    321 ...         % gamma_pot
    150 ...         % tau       syn plast time cst  (ms)
    ];

switch nargin
    case 0
        pre_spikes_hist = [];
        post_spikes_hist = [];
        params = def_params;
        int_scheme = 'euler_expl';
        int_step = 0.1;
        noise_lvl = 0;
    case 1
        post_spikes_hist = [];
        params = def_params;
        int_scheme = 'euler_expl';
        int_step = 0.1;
        noise_lvl = 0;
    case 2
        params = def_params;
        int_scheme = 'euler_expl';
        int_step = 0.1;
        noise_lvl = 0;
    case 3
        int_scheme = 'euler_expl';
        int_step = 0.1;
        noise_lvl = 0;
    case 4
        int_step = 0.1;
        noise_lvl = 0;
    case 5
        noise_lvl = 0;
    case 6
    otherwise
        error('6 inputs max are accepted')
end

%%%%%%%%%%%%%%%%%%%%
% Unpacking params %
%%%%%%%%%%%%%%%%%%%%

T = params(1);
rho_0 = params(2);

step = int_step;
n_steps = T / step;

C_pre = params(3);
C_post = params(4);
tau_Ca = params(5);
delay_pre = params(6);

theta_dep = params(7);
gamma_dep = params(8);

theta_pot = params(9);
gamma_pot = params(10);

tau = params(11);

%% Building events list based on calcium hypothesis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

len_pre = length(pre_spikes_hist);
len_post = length(post_spikes_hist);

evts = [];
simult_thr = 1e-6;

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

if strcmp(int_scheme, 'euler_expl')
    % Initiate simulation %
    %%%%%%%%%%%%%%%%%%%%%%%
    rho = params(2); 
    rho_hist = [];
    c_hist = [];
    t = 0;
    c = 0;

    % Check whether simulation is trivial %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if isempty(evts)
        rho_hist = rho_0 * ones(n_steps);
    else
        % Scheme propagation %
        %%%%%%%%%%%%%%%%%%%%%%
        while t < T
            c_hist = [c_hist, c];
            evt = find(abs(t - double(evts(:,1))) < simult_thr);
            if evt
                C_bump = sum(evts(evt,2));
                c = c + C_bump; 
            end

            if c > theta_pot
                rho = rho + step/tau * (gamma_pot*(1-rho) - gamma_dep*rho + noise_lvl*sqrt(tau)*randn);
            elseif c > theta_dep
                rho = rho - step/tau * gamma_dep * rho;
            end
            rho_hist = [rho_hist, rho];
            c = c * exp(-step/tau_Ca);
            t = t + step;
        end
    end     
        
else
	error('This integration scheme is not currently handled')
end

end

