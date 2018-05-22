function [STDP] = get_STDP(model, params, int_scheme, int_step)
% STDP EXPERIMENT
% - Runs a battery of model simulation with Calcium bumps
% relfecting different temporal differences. Uses those simulation to build
% a dw = f(delta_t) curve
% - All time params are in ms, all frequencies are in Hz
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
    -75 ...         % t_min
    75 ...          % t_max
    3 ...           % dt        step of d_t grid
    60 ...          % n_iter    nb of spike pairs
    1 ...           % freq                          (Hz)
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

t_min = params(12);
t_max = params(13);
dt = params(14);
n_iter = params(15);
freq = params(16);

%% Running simulations, returning STDP curve
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n_points = (t_max - t_min) / dt;
STDP = [];

for d_t = linspace(t_min, t_max, n_points)
    % Define the calcium bumps history
    if d_t > 0
        pre_spikes_hist = linspace(0, n_iter/(1000*freq), n_iter);
        post_spikes_hist = pre_spikes_hist + d_t;
    else
        post_spikes_hist = linspace(0, n_iter/(1000*freq), n_iter);
        pre_spikes_hist = post_spikes_hist + d_t;
    end
    
    % Simulate the evolution of synaptic strength through model
    if strcmp(model, 'naive')
        rho_hist = naive_model(pre_spikes_hist, post_spikes_hist, params(1:11), int_scheme, int_step);
        d_rho = rho_hist(end) - rho_hist(1);
        q_rho = rho_hist(end)/rho_hist(1);
        STDP = cat(1, STDP, [d_t, q_rho]);
    end
end

end