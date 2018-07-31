function STDP = get_npairs_heatmap(model, mode, params, int_scheme, dt_params, pairs_params)

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

switch nargin
    case 0
        model = 'naive';
        mode = 'rel';
        params = def_params;
        int_scheme = 'euler_expl';
        dt_params = [-50, 50, 2];
        pairs_params = [100, 1];
    case 1
        mode = 'rel';
        params = def_params;
        int_scheme = 'euler_expl';
        dt_params = [-50, 50, 2];
        pairs_params = [100, 1];
    case 2
        params = def_params;
        int_scheme = 'euler_expl';
        dt_params = [-50, 50, 2];
        pairs_params = [100, 1];
    case 3
        int_scheme = 'euler_expl';
        dt_params = [-50, 50, 2];
        pairs_params = [100, 1];
    case 4
        dt_params = [-50, 50, 2];
        pairs_params = [100, 1];
    case 5
        pairs_params = [100, 1];
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

    freq = params(14);
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
    
    freq = params(17);
end

dt_min = dt_params(1);
dt_max = dt_params(2);
step_dt = dt_params(3);

pairs_max = pairs_params(1);
step_pairs = pairs_params(2);

int_step = 0.5;
S_attr = 40;

%% Running simulations, returning STDP curve
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n_points_pairs = floor(pairs_max/step_pairs)-1;
n_points_dt = floor((dt_max-dt_min)/step_dt);

dts = linspace(dt_min, dt_max, n_points_dt);
pairs = linspace(1, pairs_max, n_points_pairs);

STDP = [];

for npairs_id = 1:n_points_pairs
    n_iter = pairs(1,npairs_id);
    for dt_id = 1:n_points_dt
        dt = dts(dt_id);
        T = max(1000*(n_iter-1)./freq + abs(dt) + 10*tau_Ca);
        params(1) = T;
        if dt >= 0
            pre_spikes_hist = linspace(0, 1000*(n_iter-1)/freq, n_iter);
            post_spikes_hist = pre_spikes_hist + dt;

            if strcmp(model, 'naive')
                [rho_hist, ~] = naive_model(pre_spikes_hist, post_spikes_hist, params(1:13), int_scheme, int_step);
                q_rho = rho_hist(end)/rho_0;

                if strcmp(mode, 'rel')
                    STDP = cat(1, STDP, [n_iter, dt, q_rho]);
                elseif strcmp(mode, 'abs')
                    STDP = cat(1, STDP, [n_iter, dt, rho_hist(end)]);
                elseif strcmp(mode, 'lim')
                    error('Limit mode not supported for transient mode of activity. Please lower frequency')
                else
                    error('Unknown mode')
                end
            elseif strcmp(model, 'pheno')
                [~, w_end, ~] = pheno_model_efficient(pre_spikes_hist, post_spikes_hist, params(1:16), int_scheme, int_step);
                q_w = w_end/w_0;

                if strcmp(mode, 'rel')
                    STDP = cat(1, STDP, [n_iter, dt, q_w]);
                elseif strcmp(mode, 'abs')
                    STDP = cat(1, STDP, [n_iter, dt, w_end]);
                elseif strcmp(mode, 'lim')
                    error('Limit mode not supported for transient mode of activity. Please lower frequency')
                else
                    error('Unknown mode')
                end 
            end

        else
            post_spikes_hist = linspace(0, 1000*(n_iter-1)/freq, n_iter);
            pre_spikes_hist = post_spikes_hist - dt;
            
            if strcmp(model, 'naive')
                [rho_hist, ~] = naive_model(pre_spikes_hist, post_spikes_hist, params(1:13), int_scheme, int_step);
                q_rho = rho_hist(end)/rho_0;

                if strcmp(mode, 'rel')
                    STDP = cat(1, STDP, [n_iter, dt, q_rho]);
                elseif strcmp(mode, 'abs')
                    STDP = cat(1, STDP, [n_iter, dt, rho_hist(end)]);
                elseif strcmp(mode, 'lim')
                    error('Limit mode not supported for transient mode of activity. Please lower frequency')
                else
                    error('Unknown mode')
                end
            elseif strcmp(model, 'pheno')
                [~, w_end, ~] = pheno_model_efficient(pre_spikes_hist, post_spikes_hist, params(1:16), int_scheme, int_step);
                q_w = w_end/w_0;

                if strcmp(mode, 'rel')
                    STDP = cat(1, STDP, [n_iter, dt, q_w]);
                elseif strcmp(mode, 'abs')
                    STDP = cat(1, STDP, [n_iter, dt, w_end]);
                elseif strcmp(mode, 'lim')
                    error('Limit mode not supported for transient mode of activity. Please lower frequency')
                else
                    error('Unknown mode')
                end 
            end
        end
    end

end