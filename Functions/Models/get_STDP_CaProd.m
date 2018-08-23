function avgSTDP = get_STDP(STDP, params)
% STDP EXPERIMENT
% - Runs a battery of model simulation with Calcium bumps
% reflecting different temporal differences. Uses those simulation to build
% a dw = f(delta_t) curve
% - All time params are in ms, all frequencies are in Hz
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Default parameter values + unpacking params
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch nargin
    case 0
        error('1 inputs min are accepted')
    case 1
        params = default_params();
    case 2
    otherwise
        error('2 inputs max are accepted')
end

%%%%%%%%%%%%%%%%%%%%
% Unpacking params %
%%%%%%%%%%%%%%%%%%%%

model = STDP.model;
mode = STDP.mode;
T = STDP.T;
step = STDP.int_step;
scheme = STDP.int_scheme;
n_steps = T / step;

rho_max = params.rho_max;
C_pre = params.C_pre;
C_post = params.C_post;
tau_Ca = params.tau_Ca;
delay_pre = params.delay_pre;
theta_dep = params.theta_dep;
gamma_dep = params.gamma_dep;
theta_pot = params.theta_pot;
gamma_pot = params.gamma_pot;
tau_rho = params.tau_rho;
sigma = params.noise_lvl;
S_attr = params.S_attr;
tau_w = params.tau_w;
theta_act = params.theta_act;
S_attr = params.S_attr;

t_min = STDP.dt.min;
t_max = STDP.dt.max;
dt = STDP.dt.step;
n_iter = STDP.n_iter;
freq = STDP.frequency;

rho0_step = 5;

prot = params;
prot.n_iter = STDP.n_iter;
prot.frequency = STDP.frequency;

%% Running simulations, returning STDP curve
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n_points = 1 + floor((t_max - t_min)/dt);

% perm_regime = (freq/1000 < 1/(t_max + 10*tau_Ca));
perm_regime = 0;

function r = Ca_topTheta_rate(theta, dt)

% PRODUCT MODELS
%     if C_pre*C_post >= theta
%         dt_crit = log(C_pre*C_post/theta);
%         r = tau_Ca * (freq/1000) * max(0, dt_crit - abs(dt));
%     else
%         r = 0;
%     end

% ADDITIVE MODELS
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

indivSTDP = ones(1+floor((t_max-t_min)/dt), 1+floor(rho_max/rho0_step));
muW = 0.5;

% Can we compute the STDP curve explicitely?
if perm_regime
        
    % If so, compute rate of time spent above thresholds...
    t = linspace(t_min, t_max, n_points);
    r_pot = Ca_topTheta_rate(theta_pot, t-delay_pre);
    prot.alpha_pot = r_pot;
    r_dep = Ca_topTheta_rate(theta_dep, t-delay_pre) - r_pot;
    prot.alpha_dep = r_dep;
    % ...then get the analytic STDP curve
    a = exp(-(r_dep*gamma_dep + r_pot*(gamma_dep+gamma_pot))/((freq/1000)*tau_rho));
    b = rho_max*(gamma_pot/(gamma_pot + gamma_dep)) * exp(-(r_dep*gamma_dep)/(tau_rho*(freq/1000))) .* (1 - exp(-(r_pot*(gamma_pot+gamma_dep))/(tau_rho*(freq/1000))));
    c = sigma * sqrt((r_pot + r_dep)./(tau_rho*freq));

    rho_lim = b ./ (1-a);
        
    for rho_0 = 0:rho0_step:rho_max
        rho0_id = 1 + floor(rho_0/rho0_step);
        w_0 = transfer(rho_0, prot);
        rho_lim(isnan(rho_lim)) = rho_0;
        rho = (rho_0 - rho_lim).*(a.^n_iter)...
            + rho_lim...
            + c .* sqrt((1 - a.^(2*n_iter))./(1 - a.^2)) .* randn(1,n_points); % final EPSP amplitude

       rho(isnan(rho)) = (rho_0 - rho_lim(isnan(rho))).*(a(isnan(rho)).^n_iter) + rho_lim(isnan(rho));
        if strcmp(mode, 'rel')
            indivSTDP(:,rho0_id) = ((1./w_0).*transfer(rho, prot))';
        elseif strcmp(mode, 'abs')
            indivSTDP(:,rho0_id) = transfer(rho, prot)';
        elseif strcmp(mode, 'lim')
            indivSTDP(:,rho0_id) = transfer(rho_lim, prot)';
        else
            error('Unknown mode')
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
        STDP.T = 1000*(n_iter-1)/freq + 10*tau_Ca;
        for rho_0 = 0:rho0_step:rho_max
            rho0_id = 1 + floor(rho_0/rho0_step);
            w_0 = transfer(rho_0, prot);
            params.rho_0 = rho_0;
            params.w_0 = w_0;
            [~, w_end, ~] = caProd_model_efficient(pre_spikes_hist, post_spikes_hist, params, STDP);
            q_w = w_end/w_0;
            indivSTDP(t_id, rho0_id) = q_w;
        end
    end
    
end

tmp = transfer_ind(0:rho0_step:rho_max, prot);
rhoDistr = (1/(sqrt(2*pi*muW*(1-muW)))) .* exp(-(tmp - muW*ones(1, 1+floor(rho_max/rho0_step))).^2 ./(2*muW*(1-muW))); 
rhoDistr = (1/sum(rhoDistr)).*rhoDistr;
avgSTDP = indivSTDP * rhoDistr';
avgSTDP = cat(2,avgSTDP,(t_min:dt:t_max)');
avgSTDP(:,[1 2]) = avgSTDP(:, [2 1]);

end