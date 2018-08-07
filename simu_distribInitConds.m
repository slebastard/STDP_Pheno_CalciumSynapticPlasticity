%% SIMULATION STARTER INSTRUCTIONS %%
% 1) Pick a mode among the following:
% -------------------------------------------------------------------------
% Simulation modes
% - single    Full evolution of syn plast on a single simulation
% - STDP      Provides STDP curve for provided frequency and # of pairings
% - freq3     STDP = f(freq,dt) 3D plot
% - freq      STDP = f(freq) 3D plot for two opposite timings
% - pairs3    STDP = f(n_iter,dt) 3D plot
% - pairs     STDP = f(n_iter) 3D plot for two opposite timings
%
% - dataFit
% -------------------------------------------------------------------------
%
% 2) Choose the model to use for simulation:
% -------------------------------------------------------------------------
% Models
% - naive   Graupner & Brunel 2012 with plasticity break (stopped below
%           calcium threshold, and noise
% - pheno   Graupner & Brunel 2012 with plasticity brake, including a
%           mapping of noise from  the analysis of G&B 2007 realistic model
% -------------------------------------------------------------------------
%
% 3) Depending on which modes you use, you will want to change the default
% dt, number of pairings and frequency in the environment section below.
% Mode-specific parameters are availabel in pertaining section
%
% You should be all set!

mode = 'STDP';
model = 'pheno';

% Parameters controlling excitation history
d_t = 30;
n_iter = 100;
frequency = 1;

%% Environment definition

T = 10;
rho_max = 200;
S_attr = 40;

C_pre = 0.4;
C_post = 0.8;
tau_Ca = 80;
delay_pre = -15;
Ca_params = [C_pre, C_post, tau_Ca, delay_pre];

theta_dep = 1;
gamma_dep = 200;
dep_params = [theta_dep, gamma_dep];

theta_pot = 1.05;
% gamma_pot = balance_coefs(C_pre, C_post, theta_dep, theta_pot, gamma_dep);
gamma_pot = 120;
pot_params = [theta_pot, gamma_pot];

N_A = 6.02e17; %mumol^(-1)
V = 2.5e-16; %L
theta_act = theta_dep;

tau_rho = 50000;
tau_w = 500000;

noise_lvl = 25; %1/sqrt(N_A*V);

% Initialization of variables
rho_0 = 25; % must be between 0 and rho_max

if strcmp(model, 'naive')
    model_params = [T, rho_0, rho_max, Ca_params, dep_params, pot_params, tau_rho, noise_lvl];
elseif strcmp(model, 'pheno')
    model_params = [T, rho_0, rho_max, Ca_params, dep_params, pot_params, tau_rho, noise_lvl, tau_w, theta_act];
end

int_scheme = 'euler_expl';
scheme_step = 0.5;

% Defining default excitation timeline
pre_spikes_hist = linspace(0, 1000*(n_iter-1)./frequency, n_iter);
post_spikes_hist = pre_spikes_hist + d_t;
T = max(1000*(n_iter-1)./frequency + abs(d_t) + 10*tau_Ca);
model_params(1) = T;

%% mode='single') Full evolution of syn plast on a single simulation
if strcmp(mode, 'single') || strcmp(mode, 'all')
    if strcmp(model, 'naive')
        [rho_hist, c_hist] = naive_model(pre_spikes_hist, post_spikes_hist, model_params, int_scheme, scheme_step);
    elseif strcmp(model, 'pheno')
        [rho_hist, w_hist, c_hist] = pheno_model(pre_spikes_hist, post_spikes_hist, model_params, int_scheme, scheme_step);
    end

    % Plotting rho as a function of time
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    t = linspace(0, T, T/scheme_step + 1);

    figure(1)
    plot(t, rho_hist);
    title('Evolution of average CaMKII state');
    xlabel('Time');
    ylabel('Average CaMKII state');

    % ToDo: add bumps of Ca as colored pins over x-axis

    figure(2)
    plot(t, c_hist);
    title('Evolution of calcium influx');
    xlabel('Time');
    ylabel('Calcium concentration');

    dep_thr = refline([0 theta_dep]);
    dep_thr.Color = 'r';

    pot_thr = refline([0 theta_pot]);
    pot_thr.Color = 'g';
    
    act_thr = refline([0 theta_act]);
    act_thr.Color = 'm';
    
    if strcmp(model, 'pheno')
        figure(3)
        plot(t, w_hist);
        title('Evolution of synaptic strength')
        xlabel('Time');
        ylabel('Average synaptic strength');
    end
end


if strcmp(mode, 'efficient')
    if strcmp(model, 'pheno')
        [rho_int, w_hist, c_int] = pheno_model(pre_spikes_hist, post_spikes_hist, model_params, int_scheme, scheme_step);
        [rho_hist, w_end, c_hist] = pheno_model_efficient(pre_spikes_hist, post_spikes_hist, model_params, int_scheme, scheme_step);
    end

    % Plotting rho as a function of time
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    t = linspace(0, T, T/scheme_step + 1);

    figure(1)
    plot(t, rho_int, 'r');
    hold on
    plot(rho_hist(:,1), rho_hist(:,2), 'xg');
    title('Evolution of average CaMKII state');
    xlabel('Time');
    ylabel('Average CaMKII state');

    % ToDo: add bumps of Ca as colored pins over x-axis

    figure(2)
    plot(t, c_int, 'r');
    hold on
    plot(c_hist(:,1), c_hist(:,2), 'xg')
    title('Evolution of calcium influx');
    xlabel('Time');
    ylabel('Calcium concentration');

    dep_thr = refline([0 theta_dep]);
    dep_thr.Color = 'r';

    pot_thr = refline([0 theta_pot]);
    pot_thr.Color = 'g';
    
    act_thr = refline([0 theta_act]);
    act_thr.Color = 'm';
end

%% mode='STDP') Provides a STDP curve for provided frequency and number of pairings

% Obtaining STDP curve
%%%%%%%%%%%%%%%%%%%%%%
if strcmp(mode, 'STDP') || strcmp(mode, 'all')
    t_min = -100;
    t_max = 100;
    dt = 0.1;

    % Initialization of variables
    rho0_step = 5;
    n_rho0 = 1+floor(rho_max/rho0_step);
    n_dt = 1+floor((t_max - t_min)/dt);
    STDP = [];
    
    muW = 0.5;
    N = 1000;
    
    for rho_0 = 0:rho0_step:rho_max
        model_params(2) = rho_0;
        stdp_params = [model_params, t_min, t_max, dt, n_iter, frequency];
        tmp = get_STDP(model, 'rel', stdp_params, int_scheme, scheme_step);
        STDP = cat(2, STDP, tmp(:,2));
    end
    
    tmp = transfer(0:rho0_step:rho_max, S_attr, noise_lvl);
    rhoDistr = (1/(sqrt(2*pi*muW*(1-muW)))) .* exp(-(tmp - muRho*ones(1, 1+floor(rho_max/rho0_step))).^2 ./(2*muW*(1-muW))); 
    rhoDistr = (1/sum(rhoDistr)).*rhoDistr;
    avgSTDP = STDP * rhoDistr';
    avgSTDP = cat(2,avgSTDP,(t_min:dt:t_max)');
    avgSTDP(:,[1 2]) = avgSTDP(:, [2 1]);
    
    % Validation against simulation
    % [STDP_an, STDP_sim] = get_both_STDP(model, 'rel', stdp_params, int_scheme, scheme_step);

    figure(4)
    plot(avgSTDP(:,1), avgSTDP(:,2), '.b');
    % plot(STDP_an(:,1), STDP_an(:,2), '+r');
    % hold on
    % plot(STDP_sim(:,1), STDP_sim(:,2), 'xg');
   

    neutral_hline = refline([0 1]);
    neutral_hline.Color = 'b';
    
    title('Plasticity as a function of pre-post spike delay')
    xlabel('Pre-post spike delay (ms)')
    ylabel('Relative change in synaptic strength')
%     YL = get(gca,'ylim');
%     YL(1) = 2 - YL(2);
%     set(gca, 'ylim', YL)
    
end

%% mode='freq3') STDP = f(freq,dt) 3D plot

if strcmp(mode, 'freq3')
    dtmin = -75;
    dtmax = 75;
    step_dt = 1;
    dt_params=[dtmin, dtmax, step_dt];
    
    freq_max = 40;
    freq_step = 0.5;
    freq_params = [freq_max, freq_step];
    
    heatmap_params = [model_params, n_iter];
    frq_map = get_freq_heatmap(model, 'rel', heatmap_params, int_scheme, dt_params, freq_params);
    
    n_freq = 1+floor((freq_max-1)/freq_step);
    n_dt = 1+floor((dtmax-dtmin)/step_dt);
    
    frq_heat = zeros(n_freq, n_dt);
    frq_heat(sub2ind([n_freq,n_dt], repelem(1:n_freq,1,n_dt), repmat(1:n_dt,1,n_freq))) = frq_map(:,3);
    
    figure(5)
    
%     scatter3(frq_map(:,1),frq_map(:,2),frq_map(:,3), '.');
    
    imagesc(frq_heat');
    colormap('hot');
    colorbar;
    
    xlabels = 1 + freq_step.*(xticks-1);
    ylabels = dtmin + step_dt.*(yticks-1);
    
    xtickformat('%.1f')
    set(gca, 'XTickLabel', xlabels);
    set(gca, 'YTickLabel', ylabels);

    title('Relative change in syn plast as a function of frequency and dt');
    xlabel('Frequency');
    ylabel('dt');
end

%% mode='freq') STDP = f(freq) 3D plot for two opposite timings
if strcmp(mode, 'freq')
    max_freq = 75;
    step_freq = 0.5;
    freq_params = [model_params, n_iter];

    [freq_prepost, freq_postpre] = get_freqSTDP(model, 'rel', freq_params, int_scheme, d_t, max_freq, step_freq);
    
    figure(6)
    plot(freq_prepost(:,1), freq_prepost(:,2), '+g')
    
    hold on
    plot(freq_postpre(:,1), freq_postpre(:,2), '+r')
    
    xlabel 'Frequency';
    ylabel 'STDP';
end

%% mode='pairs3') STDP = f(n_iter,dt) 3D plot

if strcmp(mode, 'pairs3')
    dtmin = -100;
    dtmax = 100;
    step_dt = 2;
    dt_params=[dtmin, dtmax, step_dt];
    
    pairs_max = 100;
    pairs_step = 1;
    pairs_params = [pairs_max, pairs_step];
    
    heatmap_params = [model_params, frequency];
    pairs_map = get_npairs_heatmap(model, 'rel', heatmap_params, int_scheme, dt_params, pairs_params);
    
    n_pairs = 1+floor((pairs_max-1)/pairs_step);
    n_dt = 1+floor((dtmax-dtmin)/step_dt);
    
    pairs_heat = zeros(n_pairs, n_dt);
    pairs_heat(sub2ind([n_pairs,n_dt], repelem(1:n_pairs,1,n_dt), repmat(1:n_dt,1,n_pairs))) = pairs_map(:,3);
    
    figure(7)
    
    % scatter3(pairs_map(:,1),pairs_map(:,2),pairs_map(:,3), '.');
    
    imagesc(pairs_heat');
    colorbar;
    
    xlabels = 1 + pairs_step.*(xticks-1);
    ylabels = dtmin + step_dt.*(yticks-1);
    
    xtickformat('%.1f')
    set(gca, 'XTickLabel', xlabels);
    set(gca, 'YTickLabel', ylabels);
    
    title('Relative change in syn plast as a function of number of pairings and dt');
    xlabel('Number of pairings');
    ylabel('dt');
end

%% mode='pairs') STDP = f(n_iter) 3D plot for two opposite timings
if strcmp(mode, 'pairs')
    max_pairs = 200;
    step_pairs = 1;
    
    pairs_params = [model_params, frequency];
    [numpairs_prepost, numpairs_postpre] = get_pairsSTDP(model, 'rel', pairs_params, int_scheme, d_t, max_pairs, step_pairs);
    
    figure(8)
    plot(numpairs_prepost(:,1), numpairs_prepost(:,2), '+g')
    xlabel 'Number of pairings';
    ylabel 'STDP';
    hold on
    plot(numpairs_postpre(:,1), numpairs_postpre(:,2), '+r')

end

%% Fitting to data from Venance lab

% Comparing the model to STDP=f(freq,dt) data from L. Venance (INSERM)
freq_data = csvread('STDP_Frequency.csv',1,0);
n_data = size(freq_data,1);

range.S_attr.min = 20; % This will be a function of the other parameters as determined by XPP
range.S_attr.max = 60; % This will be a function of the other parameters as determined by XPP
range.S_attr.step = 5;

range.C_post.min = 0.6;
range.C_post.max = 1.8;
range.C_post.step = 0.4;

range.tau_Ca.min = 10;
range.tau_Ca.max = 60;
range.tau_Ca.step = 10;

range.delay_pre.min = -25;
range.delay_pre.max = 25;
range.delay_pre.step = 5;

range.theta_dep.min = 1;
range.theta_dep.max = 1;
range.theta_dep.step = 1;

range.gamma_dep.min = 200;
range.gamma_dep.max = 200;
range.gamma_dep.step = 200;

range.theta_pot.min = 1.3;
range.theta_pot.max = 1.3;
range.theta_pot.step = 1.3;

range.gamma_pot.min = 155;
range.gamma_pot.max = 155;
range.gamma_pot.step = 155;

range.tau_rho.min = 150000;
range.tau_rho.max = 150000;
range.tau_rho.step = 150000;


noise_lvl = 30; %1/sqrt(N_A*V);

if strcmp(mode, 'dataFit')
    dtmin = -100;
    dtmax = 100;
    step_dt = 1;
    dt_params=[dtmin, dtmax, step_dt];
    
    freq_max = 10;
    freq_step = 0.5;
    freq_params = [freq_max, freq_step];
    
    heatmap_params = [model_params, n_iter];
%     freq_htmp = get_freq_heatmap(model, 'rel', heatmap_params, int_scheme, dt_params, freq_params);
    
%     scatter3(freq_htmp(:,1), freq_htmp(:,2), freq_htmp(:,3), '.g')
    
%     [freq_grid, dt_grid] = meshgrid(1:freq_step:freq_max, dtmin:step_dt:dtmax);
%     STDP_interpol = griddata(freq_htmp(:,1), freq_htmp(:,2), freq_htmp(:,3), freq_grid, dt_grid);
%     surf(freq_grid, dt_grid, STDP_interpol);
%     alpha 0.3
    
    data1Hz = freq_data(floor(freq_data(:,5))==1,:);
    
    stdp_params = [model_params, dtmin, dtmax, step_dt, 100, 1];
    STDP = get_STDP(model, 'rel', stdp_params, int_scheme, scheme_step);
    plot(STDP(:,1), STDP(:,2), '.b')
    hold on
    plot(data1Hz(:,2),data1Hz(:,3)./100,'xr')
    neutral_hline = refline([0 1]);
    neutral_hline.Color = 'b';   
    ('Plasticity as a function of pre-post spike delay');
    xlabel('Pre-post spike delay (ms)');
    ylabel('Relative change in synaptic strength');
    
end


%% Functions and tools

function gamma_pot = balance_coefs(Cpre, Cpost, theta_dep, theta_pot, gamma_dep)
    default = 1;
    if Cpre < theta_dep && Cpost < theta_dep
        r = default;
    elseif Cpre < theta_dep && theta_pot < Cpost
        r = log(Cpost/theta_dep)/log(Cpost/theta_pot);
    elseif theta_dep < Cpre &&  Cpre < theta_pot && theta_pot < Cpost
        r = (log(Cpre/theta_dep) + log(Cpost/theta_dep))/log(Cpost/theta_pot);
    else
        r = (log(Cpre/theta_dep) + log(Cpost/theta_dep))/(log(Cpre/theta_pot) + log(Cpost/theta_pot));
    end
    gamma_pot = r*gamma_dep;
end