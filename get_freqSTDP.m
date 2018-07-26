function [STDP_prepost, STDP_postpre] = get_freqSTDP(model, mode, params, int_scheme, dt, freq_max, step)
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
        dt = 10;
        freq_max = 100;
        step = 0.5;
    case 1
        mode = 'rel';
        params = def_params;
        int_scheme = 'euler_expl';
        dt = 10;
        freq_max = 100;
        step = 0.5;
    case 2
        params = def_params;
        int_scheme = 'euler_expl';
        dt = 10;
        freq_max = 100;
        step = 0.5;
    case 3
        int_scheme = 'euler_expl';
        dt = 10;
        freq_max = 100;
        step = 0.5;
    case 4
        dt = 10;
        freq_max = 100;
        step = 0.5;
    case 5
        freq_max = 100;
        step = 0.5;
    case 6
        step = 0.5;
    case 7
    otherwise
        error('7 inputs max are accepted')
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

int_step = 0.5;
S_attr = 40;

%% Running simulations, returning STDP curve
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n_points = floor(freq_max/step)-1;
STDP_prepost = [];
STDP_postpre = [];

freq = linspace(1, freq_max, n_points);
perm_regime = (freq./1000 < 1./(dt + 10*tau_Ca));

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

% First let's deal with low frequencies, that result in no accumulation.
% We have an analytic formula for these
r_pot = Ca_topTheta_rate(theta_pot, dt-delay_pre);
r_dep = Ca_topTheta_rate(theta_dep, dt-delay_pre) - r_pot;

a = exp(-(r_dep*gamma_dep + r_pot*(gamma_dep+gamma_pot))/((freq/1000)*tau_rho));
b = rho_max*(gamma_pot/(gamma_pot + gamma_dep)) * exp(-(r_dep*gamma_dep)/(tau_rho*(freq/1000))) .* (1 - exp(-(r_pot*(gamma_pot+gamma_dep))/(tau_rho*(freq/1000))));
c = sigma * sqrt((r_pot + r_dep)./(tau_rho*freq));

rho_prepost_lim = b ./ (1-a);
rho_prepost_lim(isnan(rho_prepost_lim)) = rho_0;
rho_prepost = (rho_0 - rho_prepost_lim).*(a.^n_iter)...
    + rho_prepost_lim...
    + c .* sqrt((1 - a.^(2.*n_iter))./(1 - a.^2)) .* randn(1,n_points); % final EPSP amplitude

rho_prepost(isnan(rho_prepost)) = rho_0;
STDP_prepost_stat = transfer(rho_prepost(perm_regime), S_attr, sigma)./w_0;
STDP_prepost_stat = [freq(perm_regime)', STDP_prepost_stat'];


r_pot = Ca_topTheta_rate(theta_pot, -dt-delay_pre);
r_dep = Ca_topTheta_rate(theta_dep, -dt-delay_pre) - r_pot;

a = exp(-(r_dep*gamma_dep + r_pot*(gamma_dep+gamma_pot))/((freq/1000)*tau_rho));
b = rho_max*(gamma_pot/(gamma_pot + gamma_dep)) * exp(-(r_dep*gamma_dep)/(tau_rho*(freq/1000))) .* (1 - exp(-(r_pot*(gamma_pot+gamma_dep))/(tau_rho*(freq/1000))));
c = sigma * sqrt((r_pot + r_dep)./(tau_rho*freq));

rho_postpre_lim = b ./ (1-a);
rho_postpre_lim(isnan(rho_postpre_lim)) = rho_0;
rho_postpre = (rho_0 - rho_postpre_lim).*(a.^n_iter)...
    + rho_postpre_lim...
    + c .* sqrt((1 - a.^(2.*n_iter))./(1 - a.^2)) .* randn(1,n_points); % final EPSP amplitude

rho_postpre(isnan(rho_postpre)) = rho_0;
STDP_postpre_stat = transfer(rho_postpre(perm_regime), S_attr, sigma)./w_0;
STDP_postpre_stat = [freq(perm_regime)', STDP_postpre_stat'];

% Now let's deal with frequencies that lead to accumulation
id_nonstat = find(~perm_regime);
freq_nonstat = freq(id_nonstat);
STDP_prepost_nonstat = [];
STDP_postpre_nonstat = [];

for id = linspace(1, length(freq_nonstat), length(freq_nonstat))
    frq = freq_nonstat(1,id);

    pre_spikes_hist = linspace(0, 1000*(n_iter-1)/frq, n_iter);
    post_spikes_hist = pre_spikes_hist + dt;

    if strcmp(model, 'naive')
        [rho_hist, ~] = naive_model(pre_spikes_hist, post_spikes_hist, params(1:13), int_scheme, int_step);
        q_rho = rho_hist(end)/rho_0;

        if strcmp(mode, 'rel')
            STDP_prepost_nonstat = cat(1, STDP_prepost_nonstat, [frq, q_rho]);
        elseif strcmp(mode, 'abs')
            STDP_prepost_nonstat = cat(1, STDP_prepost_nonstat, [frq, rho_hist(end)]);
        elseif strcmp(mode, 'lim')
            error('Limit mode not supported for transient mode of activity. Please lower frequency')
        else
            error('Unknown mode')
        end
    elseif strcmp(model, 'pheno')
        [~, w_hist, ~] = pheno_model(pre_spikes_hist, post_spikes_hist, params(1:16), int_scheme, int_step);
        q_w = w_hist(end)/w_0;

        if strcmp(mode, 'rel')
            STDP_prepost_nonstat = cat(1, STDP_prepost_nonstat, [frq, q_w]);
        elseif strcmp(mode, 'abs')
            STDP_prepost_nonstat = cat(1, STDP_prepost_nonstat, [frq, w_hist(end)]);
        elseif strcmp(mode, 'lim')
            error('Limit mode not supported for transient mode of activity. Please lower frequency')
        else
            error('Unknown mode')
        end 
    end


    post_spikes_hist = linspace(0, 1000*(n_iter-1)/frq, n_iter);
    pre_spikes_hist = post_spikes_hist + dt;

    if strcmp(model, 'naive')
        [rho_hist, ~] = naive_model(pre_spikes_hist, post_spikes_hist, params(1:13), int_scheme, int_step);
        q_rho = rho_hist(end)/rho_0;

        if strcmp(mode, 'rel')
            STDP_postpre_nonstat = cat(1, STDP_postpre_nonstat, [frq, q_rho]);
        elseif strcmp(mode, 'abs')
            STDP_postpre_nonstat = cat(1, STDP_postpre_nonstat, [frq, rho_hist(end)]);
        elseif strcmp(mode, 'lim')
            error('Limit mode not supported for transient mode of activity. Please lower frequency')
        else
            error('Unknown mode')
        end
    elseif strcmp(model, 'pheno')
        [~, w_hist, ~] = pheno_model(pre_spikes_hist, post_spikes_hist, params(1:16), int_scheme, int_step);
        q_w = w_hist(end)/w_0;

        if strcmp(mode, 'rel')
            STDP_postpre_nonstat = cat(1, STDP_postpre_nonstat, [frq, q_w]);
        elseif strcmp(mode, 'abs')
            STDP_postpre_nonstat = cat(1, STDP_postpre_nonstat, [frq, w_hist(end)]);
        elseif strcmp(mode, 'lim')
            error('Limit mode not supported for transient mode of activity. Please lower frequency')
        else
            error('Unknown mode')
        end 
    end

end

STDP_prepost = [STDP_prepost_stat; STDP_prepost_nonstat];
STDP_postpre = [STDP_postpre_stat; STDP_postpre_nonstat];

end