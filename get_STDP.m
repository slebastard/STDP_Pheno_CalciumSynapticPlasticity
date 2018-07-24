function STDP = get_STDP(model, mode, params, int_scheme, int_step)
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

t_min = params(13);
t_max = params(14);
dt = params(15);
n_iter = params(16);
freq = params(17);

%% Running simulations, returning STDP curve
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rho_max = 199.8;

n_points = 1 + (t_max - t_min)/dt;
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

% Can we compute the STDP curve explicitely?
if perm_regime
    
    % If so, compute rate of time spent above thresholds...
    dt = linspace(t_min, t_max, n_points);
    r_pot = Ca_topTheta_rate(theta_pot, dt-delay_pre);
    r_dep = Ca_topTheta_rate(theta_dep, dt-delay_pre) - r_pot;
    % ...then get the analytic STDP curve
    a = exp(-(r_dep*gamma_dep + r_pot*(gamma_dep+gamma_pot))/((freq/1000)*tau));
    b = rho_max*(gamma_pot/(gamma_pot + gamma_dep)) * exp(-(r_dep*gamma_dep)/(tau*(freq/1000))) .* (1 - exp(-(r_pot*(gamma_pot+gamma_dep))/(tau*(freq/1000))));
    c = sigma * sqrt((r_pot + r_dep)./(tau*freq));
    
    rho_lim = b ./ (1-a);
    rho_lim(isnan(rho_lim)) = rho_0;
    rho = (rho_0 - rho_lim).*(a.^n_iter)...
        + rho_lim...
        + c .* sqrt((1 - a.^(2*n_iter))./(1 - a.^2)) .* randn(1,n_points); % final EPSP amplitude

    rho(isnan(rho)) = rho_0;
    
    if strcmp(mode, 'rel')
        STDP = transpose(cat(1, dt, rho/rho_0));
    elseif strcmp(mode, 'abs')
        STDP = transpose(cat(1, dt, rho));
    elseif strcmp(mode, 'lim')
        STDP = transpose(cat(1, dt, rho_lim));
    else
        error('Unknown mode')
    end
 
else
    
    % When we cannot consider that pairs of spikes are independent, we
    % compute the curve by individual simulations...
    for dt = linspace(t_min, t_max, n_points)
        % Define the calcium bumps history
        if dt > 0
            pre_spikes_hist = linspace(0, 1000*(n_iter-1)/freq, n_iter);
            post_spikes_hist = pre_spikes_hist + dt;
        else
            post_spikes_hist = linspace(0, 1000*(n_iter-1)/freq, n_iter);
            pre_spikes_hist = post_spikes_hist - dt;
        end

        % Simulate the evolution of synaptic strength through model -
        % COMPARE TO ANALYTIC
        if strcmp(model, 'naive')
            params(1) = 1000*(n_iter-1)/freq + 10*tau_Ca;
            [rho_hist, ~] = naive_model(pre_spikes_hist, post_spikes_hist, params(1:12), int_scheme, int_step);
            q_rho = rho_hist(end)/rho_hist(1);
            
            if strcmp(mode, 'rel')
                STDP = cat(1, STDP, [dt, q_rho]);
            elseif strcmp(mode, 'abs')
                STDP = cat(1, STDP, [dt, rho_hist(end)]);
            elseif strcmp(mode, 'lim')
                error('Limit mode not supported for transient mode of activity. Please lower frequency')
            else
                error('Unknown mode')
            end
           
        end
    end
    
end

end