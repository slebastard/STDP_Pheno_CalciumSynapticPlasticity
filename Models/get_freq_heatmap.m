function STDP = get_freq_heatmap(model, mode, params, int_scheme, dt_params, freq_params)
% STDP EXPERIMENT
% - Runs a battery of model simulation with Calcium bumps
% reflecting different temporal differences. Uses those simulation to build
% a dw = f(delta_t) curve
% - All time params are in ms, all frequencies are in Hz
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Default parameter values + unpacking params
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

def_params = [...
    1000 ...        % T         total simu time     (ms)
    .3 ...          % rho_0     init syn strength
    1 ...           % rho_max
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
    10 ...           % n_iter
    ];

switch nargin
    case 0
        model = 'naive';
        mode = 'rel';
        params = def_params;
        int_scheme = 'euler_expl';
        dt_params = [-50, 50, 2];
        freq_params = [100, 2];
    case 1
        mode = 'rel';
        params = def_params;
        int_scheme = 'euler_expl';
        dt_params = [-50, 50, 2];
        freq_params = [100, 2];
    case 2
        params = def_params;
        int_scheme = 'euler_expl';
        dt_params = [-50, 50, 2];
        freq_params = [100, 2];
    case 3
        int_scheme = 'euler_expl';
        dt_params = [-50, 50, 2];
        freq_params = [100, 2];
    case 4
        dt_params = [-50, 50, 2];
        freq_params = [100, 2];
    case 5
        freq_params = [100, 2];
    case 6
    otherwise
        error('6 inputs max are accepted. Please provide freq and dt parameters as arrays')
end

%%%%%%%%%%%%%%%%%%%%
% Unpacking params %
%%%%%%%%%%%%%%%%%%%%

if strcmp(model, 'naive')
    T = params(1);
    rho_0 = params(2);
    rho_max = params(3);
    C_pre = params(4);
    C_post = params(5);
    tau_Ca = params(6);
    delay_pre = params(7);

    theta_dep = params(8);
    gamma_dep = params(9);

    theta_pot = params(10);
    gamma_pot = params(11);

    tau_rho = params(12);
    sigma = params(13);

    n_iter = params(14);
elseif strcmp(model, 'pheno')
    T = params(1);
    rho_0 = params(2);
    rho_max = params(3);
    C_pre = params(4);
    C_post = params(5);
    tau_Ca = params(6);
    delay_pre = params(7);

    theta_dep = params(8);
    gamma_dep = params(9);

    theta_pot = params(10);
    gamma_pot = params(11);

    tau_rho = params(12);
    sigma = params(13);

    tau_w = params(14);
    theta_act = params(15);
    
    n_iter = params(16);
end

dt_min = dt_params(1);
dt_max = dt_params(2);
step_dt = dt_params(3);

freq_max = freq_params(1);
step_freq = freq_params(2);

int_step = 0.5;
S_attr = 40;
w_0 = transfer(rho_0, S_attr, sigma);

%% Running simulations, returning STDP curve
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n_points_freq = 1+floor((freq_max-1)/step_freq);
freqs = linspace(1, freq_max, n_points_freq);

STDP = [];

for freq_id = 1:n_points_freq
    frq = freqs(1,freq_id);
    stdp_params = [params(1:15), dt_params, params(16), frq];
    std = get_STDP(model, 'rel', stdp_params, int_scheme, 0.5);
    std = cat(2, frq*ones(size(std,1),1), std);
    STDP = cat(1, STDP, std);
end

end