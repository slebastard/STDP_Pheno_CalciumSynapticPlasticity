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

    w_0 = params(14);
    tau_w = params(15);
    theta_act = params(16);
    
    n_iter = params(17);
end

dt_min = dt_params(1);
dt_max = dt_params(2);
step_dt = dt_params(3);

freq_max = freq_params(1);
step_freq = freq_params(2);

int_step = 0.5;
S_attr = 40;

%% Running simulations, returning STDP curve
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n_points_freq = floor(freq_max/step_freq)-1;
n_points_dt = floor((dt_max-dt_min)/step_dt);

dts = linspace(dt_min, dt_max, n_points_dt);
freqs = linspace(1, freq_max, n_points_freq);

STDP = [];

for freq_id = 1:n_points_freq
    frq = freqs(1,freq_id);
    for dt_id = 1:n_points_dt
        dt = dts(dt_id);
        if dt >= 0
            pre_spikes_hist = linspace(0, 1000*(n_iter-1)/frq, n_iter);
            post_spikes_hist = pre_spikes_hist + dt;
            T = max(1000*(n_iter-1)./frq + abs(dt) + 10*tau_w);
            params(1) = T;
            if strcmp(model, 'naive')
                [rho_hist, ~] = naive_model(pre_spikes_hist, post_spikes_hist, params(1:13), int_scheme, int_step);
                q_rho = rho_hist(end)/rho_0;

                if strcmp(mode, 'rel')
                    STDP = cat(1, STDP, [frq, dt, q_rho]);
                elseif strcmp(mode, 'abs')
                    STDP = cat(1, STDP, [frq, dt, rho_hist(end)]);
                elseif strcmp(mode, 'lim')
                    error('Limit mode not supported for transient mode of activity. Please lower frequency')
                else
                    error('Unknown mode')
                end
            elseif strcmp(model, 'pheno')
                [~, w_hist, ~] = pheno_model(pre_spikes_hist, post_spikes_hist, params(1:16), int_scheme, int_step);
                q_w = w_hist(end)/w_0;

                if strcmp(mode, 'rel')
                    STDP = cat(1, STDP, [frq, dt, q_w]);
                elseif strcmp(mode, 'abs')
                    STDP = cat(1, STDP, [frq, dt, w_hist(end)]);
                elseif strcmp(mode, 'lim')
                    error('Limit mode not supported for transient mode of activity. Please lower frequency')
                else
                    error('Unknown mode')
                end 
            end

        else
            post_spikes_hist = linspace(0, 1000*(n_iter-1)/frq, n_iter);
            pre_spikes_hist = post_spikes_hist + dt;
            T = max(1000*(n_iter-1)./frq + abs(dt) + 10*tau_w);
            params(1) = T;
            if strcmp(model, 'naive')
                [rho_hist, ~] = naive_model(pre_spikes_hist, post_spikes_hist, params(1:13), int_scheme, int_step);
                q_rho = rho_hist(end)/rho_0;

                if strcmp(mode, 'rel')
                    STDP = cat(1, STDP, [frq, dt, q_rho]);
                elseif strcmp(mode, 'abs')
                    STDP = cat(1, STDP, [frq, dt, rho_hist(end)]);
                elseif strcmp(mode, 'lim')
                    error('Limit mode not supported for transient mode of activity. Please lower frequency')
                else
                    error('Unknown mode')
                end
            elseif strcmp(model, 'pheno')
                [~, w_hist, ~] = pheno_model(pre_spikes_hist, post_spikes_hist, params(1:16), int_scheme, int_step);
                q_w = w_hist(end)/w_0;

                if strcmp(mode, 'rel')
                    STDP = cat(1, STDP, [frq, dt, q_w]);
                elseif strcmp(mode, 'abs')
                    STDP = cat(1, STDP, [frq, dt, w_hist(end)]);
                elseif strcmp(mode, 'lim')
                    error('Limit mode not supported for transient mode of activity. Please lower frequency')
                else
                    error('Unknown mode')
                end 
            end
        end
    end
end

end