function STDP = get_freq_heatmap(model, model_params, heatmap_params, int_scheme, int_step)

% STDP - HEATMAP OF SYNAPTIC PLASTICITY AS A FUNCTION OF dt AND frequency


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
    0 ...           % sigma
    ];

def_heatmap_params = [...
    -75 ...         % t_min
    75 ...          % t_max
    3               % dt
    1 ...           % freq_min
    200 ...         % freq_max
    5 ...           % freq_step
    60 ...          % n_iter
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
sigma = model_params(12);

t_min = heatmap_params(1);
t_max = heatmap_params(2);
dt = heatmap_params(3);
freq_min = heatmap_params(4);
freq_max = heatmap_params(5);
freq_step = heatmap_params(6);
n_iter = heatmap_params(6);

n_dt = (t_max - t_min) / dt;
n_freq = (freq_max - 1) / freq_step;

%% Running all simulations
%%%%%%%%%%%%%%%%%%%%%

STDP = [];

for freq = linspace(1, n_freq, n_freq)
    % 0) Compute the simulation parameters
    T = n_iter/(1000*freq) + 200;

    % 1) Simulation
    
    if strcmp(model, 'naive')
        stdp_params = [model_params, dt_min, dt_max, dt, n_iter, freq];
        STDP_dt = cat(2, get_STDP(model, stdp_params, int_scheme, int_step), freq*ones(n_dt));
        STDP = cat(1, STDP, STDP_dt);
    end
end

end