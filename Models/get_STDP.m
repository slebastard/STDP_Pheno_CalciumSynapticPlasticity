function avgSTDP = get_STDP(model, mode, params, int_scheme, int_step)
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
    -75 ...         % t_min
    75 ...          % t_max
    3 ...           % dt        step of d_t grid
    60 ...          % n_iter    nb of spike pairs
    1 ...           % freq                          (Hz)
    ];

switch nargin
    case 0
        model = 'naive';
        mode = 'rel';
        params = def_params;
        int_scheme = 'euler_expl';
        int_step = 0.1;
    case 1
        mode = 'rel';
        params = def_params;
        int_scheme = 'euler_expl';
        int_step = 0.1;
    case 2
        params = def_params;
        int_scheme = 'euler_expl';
        int_step = 0.1;
    case 3
        int_scheme = 'euler_expl';
        int_step = 0.1;
    case 4
        int_step = 0.1;
    case 5
    otherwise
        error('5 inputs max are accepted')
end

%%%%%%%%%%%%%%%%%%%%
% Unpacking params %
%%%%%%%%%%%%%%%%%%%%

if strcmp(model, 'naive')
    T = params(1);
    step = int_step;
    n_steps = T / step;

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

    t_min = params(14);
    t_max = params(15);
    dt = params(16);
    n_iter = params(17);
    freq = params(18);
elseif strcmp(model, 'pheno')
    T = params(1);
    step = int_step;
    n_steps = T / step;

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
    
    t_min = params(16);
    t_max = params(17);
    dt = params(18);
    n_iter = params(19);
    freq = params(20);
end

S_attr = 40;
rho0_step = 5;
muW = 0.5;

%% Running simulations, returning STDP curve
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n_points = 1 + floor((t_max - t_min)/dt);
STDP = [];

perm_regime = (freq/1000 < 1/(t_max + 10*tau_Ca));

function r = Ca_topTheta_rate(theta, dt)

    if C_pre<=theta && C_post<=theta
        dt_crit_low = log((theta-C_pre)/C_post);
        dt_crit_high = log(C_pre/(theta-C_post));
        r = tau_Ca * (freq/1000) * (...
            log((C_pre+C_post*exp(dt/tau_Ca))/theta) .* (dt/tau_Ca  > dt_crit_low) .* (dt/tau_Ca  <= 0) ...
            + log((C_post+C_pre*exp(-dt/tau_Ca))/theta) .* (dt/tau_Ca  > 0) .* (dt/tau_Ca  <= dt_crit_high) ...
        );

    elseif theta<=C_pre && theta<=C_post
        dt_crit_low = log(theta/C_post);
        dt_crit_high = log(C_pre/theta);

        r = tau_Ca * (freq/1000) * (...
            ( log(C_pre/theta) + log((C_pre*exp(-dt/tau_Ca) + C_post)/theta) ) .* (dt/tau_Ca > dt_crit_high) ...
            +  (log((C_pre+C_post*exp(dt/tau_Ca))/theta)-(dt/tau_Ca)) .* (dt/tau_Ca  > dt_crit_low) .* (dt/tau_Ca  <= dt_crit_high) ...
            + ( log(C_post/theta) + log((C_post*exp(dt/tau_Ca) + C_pre)/theta) ) .* (dt/tau_Ca  <= dt_crit_low) ...
            );

    elseif C_pre<theta
        dt_crit_low = log((theta - C_pre)/C_post);
        dt_crit_high = log(theta/C_post);

        r = tau_Ca * (freq/1000) * (...
            log((C_post * exp(dt/tau_Ca) + C_pre)./(theta*exp(dt/tau_Ca))) .* (dt/tau_Ca > dt_crit_high) ...
            + (log(C_post/theta) + log((C_post*exp(dt/tau_Ca) + C_pre)/theta)) .* (dt/tau_Ca  > dt_crit_low) .* (dt/tau_Ca  <= dt_crit_high) ...
            + log(C_post/theta) .* (dt/tau_Ca  <= dt_crit_low) ...
            );

    else
        dt_crit_low = log(C_pre/theta);
        dt_crit_high = log(C_pre/(theta-C_post));

        r = tau_Ca * (freq/1000) * (...
            log((C_pre * exp(-dt/tau_Ca) + C_post)./(theta*exp(-dt/tau_Ca))) .* (dt/tau_Ca  <= dt_crit_low) ...
            + (log(C_pre/theta) + log((C_pre*exp(-dt/tau_Ca) + C_post)/theta)) .* (dt/tau_Ca  > dt_crit_low) .* (dt/tau_Ca  <= dt_crit_high) ...
            + log(C_pre/theta) .* (dt/tau_Ca > dt_crit_high) ...
            );

    end
end

STDP = ones(1+floor((t_max-t_min)/dt), 1+floor(rho_max/rho0_step));
muW = 0.5;

% Can we compute the STDP curve explicitely?
if perm_regime
        
    % If so, compute rate of time spent above thresholds...
    t = linspace(t_min, t_max, n_points);
    r_pot = Ca_topTheta_rate(theta_pot, t-delay_pre);
    r_dep = Ca_topTheta_rate(theta_dep, t-delay_pre) - r_pot;
    % ...then get the analytic STDP curve
    a = exp(-(r_dep*gamma_dep + r_pot*(gamma_dep+gamma_pot))/((freq/1000)*tau_rho));
    b = rho_max*(gamma_pot/(gamma_pot + gamma_dep)) * exp(-(r_dep*gamma_dep)/(tau_rho*(freq/1000))) .* (1 - exp(-(r_pot*(gamma_pot+gamma_dep))/(tau_rho*(freq/1000))));
    c = sigma * sqrt((r_pot + r_dep)./(tau_rho*freq));

    rho_lim = b ./ (1-a);
        
    for rho_0 = 0:rho0_step:rho_max
        rho0_id = 1 + floor(rho_0/rho0_step);
        w_0 = transfer(rho_0, S_attr, sigma);
        rho_lim(isnan(rho_lim)) = rho_0;
        rho = (rho_0 - rho_lim).*(a.^n_iter)...
            + rho_lim...
            + c .* sqrt((1 - a.^(2*n_iter))./(1 - a.^2)) .* randn(1,n_points); % final EPSP amplitude

        rho(isnan(rho)) = rho_0;

        if strcmp(model, 'naive')
            if strcmp(mode, 'rel')
                STDP(:,rho0_id) = rho/rho_0;
            elseif strcmp(mode, 'abs')
                STDP(:,rho0_id) = rho;
            elseif strcmp(mode, 'lim')
                STDP(:,rho0_id) = rho_lim;
            else
                error('Unknown mode')
            end
        elseif strcmp(model, 'pheno')
            if strcmp(mode, 'rel')
                STDP(:,rho0_id) = (1/w_0).*transfer(rho, S_attr, sigma)';
            elseif strcmp(mode, 'abs')
                STDP(:,rho0_id) = transfer(rho, S_attr, sigma)';
            elseif strcmp(mode, 'lim')
                STDP(:,rho0_id) = transfer(rho_lim, S_attr, sigma)';
            else
                error('Unknown mode')
            end        
        end
    end
 
else
    
    % When we cannot consider that pairs of spikes are independent, we
    % compute the curve by individual simulations...
    for t = linspace(t_min, t_max, n_points)
        t_id = 1 + floor(1e-8 + (t-t_min)/dt);
        % Define the calcium bumps history
        if t > 0
            pre_spikes_hist = linspace(0, 1000*(n_iter-1)/freq, n_iter);
            post_spikes_hist = pre_spikes_hist + t;
        else
            post_spikes_hist = linspace(0, 1000*(n_iter-1)/freq, n_iter);
            pre_spikes_hist = post_spikes_hist - t;
        end

        % Simulate the evolution of synaptic strength through model -
        % COMPARE TO ANALYTIC
        params(1) = 1000*(n_iter-1)/freq + 10*tau_Ca;
        for rho_0 = 0:rho0_step:rho_max
            rho0_id = 1 + floor(rho_0/rho0_step);
            w_0 = transfer(rho_0, S_attr, sigma);
            params(2) = rho_0;
            [~, w_end, ~] = pheno_model_efficient(pre_spikes_hist, post_spikes_hist, params(1:16), int_scheme, int_step);
            q_w = w_end/w_0;
            STDP(t_id, rho0_id) = q_w;
        end
    end
    
end

tmp = transfer(0:rho0_step:rho_max, S_attr, sigma);
rhoDistr = (1/(sqrt(2*pi*muW*(1-muW)))) .* exp(-(tmp - muW*ones(1, 1+floor(rho_max/rho0_step))).^2 ./(2*muW*(1-muW))); 
rhoDistr = (1/sum(rhoDistr)).*rhoDistr;
avgSTDP = STDP * rhoDistr';
avgSTDP = cat(2,avgSTDP,(t_min:dt:t_max)');
avgSTDP(:,[1 2]) = avgSTDP(:, [2 1]);

end