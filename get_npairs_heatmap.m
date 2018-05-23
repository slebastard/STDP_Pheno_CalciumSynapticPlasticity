function STDP = get_npairs_heatmap(model, model_params, heatmap_params, int_scheme, int_step)

% STDP - HEATMAP OF SYNAPTIC PLASTICITY AS A FUNCTION OF dt AND n_pairs


%% Default parameter values + unpacking model_params
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

def_model_params = [...
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

def_heatmap_params = [...
    -75 ...         % dt_min
    75 ...          % dt_max
    3               % dt_step
    100 ...         % n_iter_max
    2 ...           % n_iter_step
    60 ...          % freq
    ];

switch nargin
    case 0
        model = 'naive';
        model_params = def_model_params;
        heatmap_params = def_heatmap_params;
        int_scheme = 'euler_expl';
        int_step = 0.1;
        fprintf('All arguments set to default values')
    case 1
        model_params = def_model_params;
        heatmap_params = def_heatmap_params;
        int_scheme = 'euler_expl';
        int_step = 0.1;
    case 2
        heatmap_params = def_heatmap_params;
        int_scheme = 'euler_expl';
        int_step = 0.1;
    case 3
        int_scheme = 'euler_expl';
        int_step = 0.1;
    case 4
        int_step = 0.1;
    case 5
    otherwise
        error('5 input groups max are accepted')
end

%%%%%%%%%%%%%%%%%%%%
% Unpacking model_params %
%%%%%%%%%%%%%%%%%%%%

% model_params
T = model_params(1);
step = int_step;
n_steps = T / step;

rho_0 = model_params(2);
C_pre = model_params(3);
C_post = model_params(4);
tau_Ca = model_params(5);
delay_pre = model_params(6);

theta_dep = model_params(7);
gamma_dep = model_params(8);

theta_pot = model_params(9);
gamma_pot = model_params(10);

tau = model_params(11);

dt = model_params(12);

dt_min = heatmap_params(1);
dt_max = heatmap_params(2);
dt_step = heatmap_params(3);
n_iter_max = heatmap_params(4);
n_iter_step = heatmap_params(5);
freq = heatmap_params(6);

n_dt = (t_max - t_min) / dt;
n_niter = (n_iter_max - 1) / n_iter_step;

%% Running all simulations
%%%%%%%%%%%%%%%%%%%%%

STDP = [];

for dt = linspace(t_min, t_max, n_dt)
    for n_iter = linspace(1, n_iter_max, n_niter)
        % 0) Compute the simulation parameters
        T = n_iter/(1000*freq) + 200;
        naive_params = [T, rho_0, C_pre, C_post, tau_Ca, delay_pre, theta_dep, gamma_dep, theta_pot, gamma_pot, tau];

        % 1) Simulation
        pre_spikes_hist = linspace(0, n_iter/(1000*freq), n_iter);
        post_spikes_hist = pre_spikes_hist + dt;
        if strcmp(model, 'naive')
            rho_hist = naive_model(pre_spikes_hist, post_spikes_hist, naive_params, int_scheme, int_step);
            q_rho = rho_hist(end)/rho_hist(1);
            STDP = cat(1, STDP, [dt, n_iter, q_rho]);
        end
    end
end

end