function [ ...
    dw_freq_prepost, ...
    dw_freq_postpre, ...
    dw_npairs_prepost, ...
    dw_npairs_postpre ...
    ] = get_freq_an(model, params, int_scheme, int_step)

% STDP - IMPACT OF FREQUENCY & NB PAIRING
% - Runs a battery of model simulation with Calcium bumps
% relfecting different temporal differences. Uses those simulation to build
% a dw = f(delta_t) curve
% - All time params are in ms, all frequencies are in Hz

%% Default parameter values + unpacking params
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
    2.85 ...        % sigma     noise level
    -75 ...         % dt
    1 ...           % freq_def                      (Hz)
    75 ...          % n_iter_max
    3 ...           % n_iter_def
    60 ...          % freq_min  min freq for anls   (Hz)
    1 ...           % freq_max  max freq for anls   (Hz)
    ];

switch nargin
    case 0
        model = 'naive';
        params = def_params;
        int_scheme = 'euler_expl';
        int_step = 0.1;
    case 1
        params = def_params;
        int_scheme = 'euler_expl';
        int_step = 0.1;
    case 2
        int_scheme = 'euler_expl';
        int_step = 0.1;
    case 3
        int_step = 0.1;
    case 4
    otherwise
        error('4 inputs max are accepted')
end

%%%%%%%%%%%%%%%%%%%%
% Unpacking params %
%%%%%%%%%%%%%%%%%%%%

T = params(1);
step = int_step;
n_steps = T / step;

rho_0 = params(2);
C_pre = params(3);
C_post = params(4);
tau_Ca = params(5);
delay_pre = params(6);

theta_dep = params(7);
gamma_dep = params(8);

theta_pot = params(9);
gamma_pot = params(10);

tau = params(11);
sigma = params(12);

dt = params(13);

freq_def = params(14);
n_iter_max= params(15);

n_iter_def = params(16);
freq_min = params(17);
freq_max = params(18);

%% Frequency analysis
%%%%%%%%%%%%%%%%%%%%%

dw_freq_prepost = [];
dw_freq_postpre = [];

for freq = linspace(freq_min, freq_max, 50)
    % 0) Compute the simulation parameters
    T = 1000*(n_iter_def-1)/freq + 5*tau;
    naive_params = [T, rho_0, C_pre, C_post, tau_Ca, delay_pre, theta_dep, gamma_dep, theta_pot, gamma_pot, tau, sigma];
    
    % 1) Simulation for the pre-post scenario
    pre_spikes_hist = linspace(0, 1000*(n_iter_def-1)/freq, n_iter_def);
    post_spikes_hist = pre_spikes_hist + dt;
    if strcmp(model, 'naive')
        rho_hist = naive_model(pre_spikes_hist, post_spikes_hist, naive_params, int_scheme, int_step);
        q_rho = rho_hist(end) / rho_hist(1);
        dw_freq_prepost = cat(1, dw_freq_prepost, [freq, q_rho]);
    end
    
    % 2) Simulation for the post-pre scenario
    post_spikes_hist = linspace(0, 1000*(n_iter_def-1)/freq, n_iter_def);
    pre_spikes_hist = post_spikes_hist + dt;
    if strcmp(model, 'naive')
        rho_hist = naive_model(pre_spikes_hist, post_spikes_hist, naive_params, int_scheme, int_step);
        q_rho = rho_hist(end) / rho_hist(1);
        dw_freq_postpre= cat(1, dw_freq_postpre, [freq, q_rho]);
    end
end

%% Spike pairs number analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dw_npairs_prepost = [];
dw_npairs_postpre = [];

for n_iter = linspace(1, n_iter_max, n_iter_max / 2)
    % 0) Compute the simulation parameters
    T = 1000*(n_iter+1)/freq_def + 5*tau;
    naive_params = [T, rho_0, C_pre, C_post, tau_Ca, delay_pre, theta_dep, gamma_dep, theta_pot, gamma_pot, tau, sigma];
    
    % 1) Simulation for the pre-post scenario
    pre_spikes_hist = linspace(0, 1000*(n_iter-1)/freq_def, n_iter);
    post_spikes_hist = pre_spikes_hist + dt;
    if strcmp(model, 'naive')
        rho_hist = naive_model(pre_spikes_hist, post_spikes_hist, naive_params, int_scheme, int_step);
        q_rho = rho_hist(end)/rho_hist(1);
        dw_npairs_prepost = cat(1, dw_npairs_prepost, [n_iter, q_rho]);
    end
    
    % 2) Simulation for the post-pre scenario
    post_spikes_hist = linspace(0, 1000*(n_iter_def-1)/freq, n_iter_def);
    pre_spikes_hist = post_spikes_hist + dt;
    if strcmp(model, 'naive')
        rho_hist = naive_model(pre_spikes_hist, post_spikes_hist, naive_params, int_scheme, int_step);
        q_rho = rho_hist(end)/rho_hist(1);
        dw_npairs_postpre= cat(1, dw_npairs_postpre, [n_iter, q_rho]);
    end
end

end