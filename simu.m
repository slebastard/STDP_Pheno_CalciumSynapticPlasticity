%% SIMULATION STARTER %%
% 1) Just defines the model, variables and integration scheme
% 2) Runs a single simulation and returns the full evolution of synaptic
% efficacy
% 3) Returns STDP curve by running simulation for several possible timings
% 4) Analyses the impact of frequency and of number of spike pairs on the
% plasticity during an experiment
% 5) Provides a frequency vs dt heatmap
% 6) Provides a STDP = f(n_pairs) for dt=+10ms and dt=-10ms

% Simulation mode
% single    Step 2 only
% STDP      Step 3 only
% freq      Step 4 only
%
% pairs     Step 6 only
% all       Steps 2 to 4

mode = 'efficient';

freq_data = csvread('STDP_Frequency.csv',1,0);


%% 0) Define the environment

T = 10;
rho_0 = 35; % must be between 0 and
rho_max = 200;
S_attr = 40;

C_pre = 0.8;
C_post = 1.5;
tau_Ca = 20;
delay_pre = 5;
Ca_params = [C_pre, C_post, tau_Ca, delay_pre];

theta_dep = 1;
gamma_dep = 200;
dep_params = [theta_dep, gamma_dep];

theta_pot = 1.3;
% gamma_pot = balance_coefs(C_pre, C_post, theta_dep, theta_pot, gamma_dep);
gamma_pot = 200;
pot_params = [theta_pot, gamma_pot];

tau_rho = 150000; % this is larger than I expected. Ask Brunel about this
tau_w = 500000;

N_A = 6.02e17; %mumol^(-1)
V = 2.5e-16; %L
theta_act = theta_dep;

% Modeling noise
noise_lvl = 20; %1/sqrt(N_A*V);
w_0 = transfer(rho_0, S_attr, noise_lvl);

n_iter = 10;
frequency = 1;

model = 'pheno'; %naive or pheno

if strcmp(model, 'naive')
    model_params = [T, rho_0, rho_max, Ca_params, dep_params, pot_params, tau_rho, noise_lvl];
elseif strcmp(model, 'pheno')
    model_params = [T, rho_0, rho_max, Ca_params, dep_params, pot_params, tau_rho, noise_lvl, w_0, tau_w, theta_act];
end

int_scheme = 'euler_expl';
scheme_step = 0.5;

%% 1) Define the stimulation history
d_t = 20;
pre_spikes_hist = linspace(0, 1000*(n_iter-1)./frequency, n_iter);
post_spikes_hist = pre_spikes_hist + d_t;
T = max(1000*(n_iter-1)./frequency + abs(d_t) + 10*tau_Ca);
model_params(1) = T;

%% 2) Full evolution of syn plast on a single simulation
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

%% 3) STDP curve

% Obtaining STDP curve
%%%%%%%%%%%%%%%%%%%%%%
if strcmp(mode, 'STDP') || strcmp(mode, 'all')
    t_min = -75;
    t_max = 75;
    dt = 0.5;

    stdp_params = [model_params, t_min, t_max, dt, n_iter, frequency];
    %STDP = get_STDP(model, 'rel', stdp_params, int_scheme, scheme_step);
    % Validation against simulation
    [STDP_an, STDP_sim] = get_both_STDP(model, 'rel', stdp_params, int_scheme, scheme_step);

    figure(3)
    % plot(STDP(:,1), STDP(:,2), '+r');
    plot(STDP_an(:,1), STDP_an(:,2), '+r');
    hold on
    plot(STDP_sim(:,1), STDP_sim(:,2), 'xg');
    
    title('Plasticity as a function of pre-post spike delay')
    xlabel('Pre-post spike delay (ms)')
    ylabel('Relative change in synaptic strength')

    neutral_hline = refline([0 1]);
    neutral_hline.Color = 'b';
    
%     YL = get(gca,'ylim');
%     YL(1) = 2 - YL(2);
%     set(gca, 'ylim', YL)
    
end

%% Analysis of frequency and number of spike pairs
% if strcmp(mode, 'freq') || strcmp(mode, 'all')
%     dt = 10;
% 
%     freq_def = 100;
%     n_iter_max = 200;
% 
%     n_iter_def = 100;
%     freq_min = 1;
%     freq_max = 50;
% 
%     freq_an_params = [model_params, dt, freq_def, n_iter_max, n_iter_def, freq_min, freq_max];
% 
%     [
%         dw_freq_prepost, ...
%         dw_freq_postpre, ...
%         dw_npairs_prepost, ...
%         dw_npairs_postpre ...    
%     ] ...
%     = get_freq_an(model, freq_an_params, int_scheme, scheme_step);
% 
% 
%     figure(4)
% 
%     plot(dw_freq_prepost(:,1), dw_freq_prepost(:,2), '+')
%     hold on
%     plot(dw_freq_postpre(:,1), dw_freq_postpre(:,2), '+')
% 
%     xlabel('Frequency (Hz)')
%     ylabel('Relative change in synaptic weight')
%     legend('dt = +10ms','dt = -10ms')
%     title('Plasticity as a function of frequency (Hz), for 60 pairings at +-10ms')
% 
%     hold off
% 
% 
%     figure(5)
% 
%     plot(dw_npairs_prepost(:,1), dw_npairs_prepost(:,2), '+')
%     hold on
%     plot(dw_npairs_postpre(:,1), dw_npairs_postpre(:,2), '+')
% 
%     xlabel('Nb of spike pairs at 1Hz');
%     ylabel('Relative change in synaptic weight');
%     legend('dt = +10ms','dt = -10ms')
%     title('Plasticity as a function of the number of spike pairs, for pairings at +-10ms, at 1Hz')
% 
%     hold off
% end
% 
%% Frequency - dt - scatter3

if strcmp(mode, 'freq3')
    dtmin = -50;
    dtmax = 50;
    step_dt = 10;
    dt_params=[dtmin, dtmax, step_dt];
    
    freq_max = 20;
    freq_step = 2;
    freq_params = [freq_max, freq_step];
    
    n_iter = 100;
    
    heatmap_params = [model_params, n_iter];
    frq_htmp = get_freq_heatmap(model, 'rel', heatmap_params, int_scheme, dt_params, freq_params);
    
    figure(6)
    scatter3(frq_htmp(:,1),frq_htmp(:,2),frq_htmp(:,3));
    title = 'Relative change in syn plast as a function of frequency and dt';
    xlabel = 'Frequency';
    ylabel = 'dt';
end

%% NIGHTRUN - Freq-dt scatter 3

if strcmp(mode, 'Nightrun_freq3')
    dtmin = -50;
    dtmax = 50;
    step_dt = 10;
    dt_params=[dtmin, dtmax, step_dt];
    
    freq_max = 20;
    freq_step = 2;
    freq_params = [freq_max, freq_step];
    
    n_iter = 100;
    
    heatmap_params = [model_params, n_iter];
    frq_htmp = get_freq_heatmap(model, 'rel', heatmap_params, int_scheme, dt_params, freq_params);
    
    figure(6)
    scatter3(frq_htmp(:,1),frq_htmp(:,2),frq_htmp(:,3));
    title = 'Relative change in syn plast as a function of frequency and dt';
    xlabel = 'Frequency';
    ylabel = 'dt';
end



%% Nb Pairs - STDP
if strcmp(mode, 'pairs')
    max_pairs = 200;
    step_pairs = 1;
    
    pairs_params = [model_params, frequency];
    [numpairs_prepost, numpairs_postpre] = get_pairsSTDP(model, 'rel', pairs_params, int_scheme, d_t, max_pairs, step_pairs);
    
    figure(4)
    plot(numpairs_prepost(:,1), numpairs_prepost(:,2), '+g')
    
    hold on
    plot(numpairs_postpre(:,1), numpairs_postpre(:,2), '+r')
    
    xlabel('Number of pairings')
    ylabel('STDP')
end

%% Freq - STDP
if strcmp(mode, 'freq')
    max_freq = 75;
    step_freq = 0.5;
    freq_params = [model_params, n_iter];

    [freq_prepost, freq_postpre] = get_freqSTDP(model, 'rel', freq_params, int_scheme, d_t, max_freq, step_freq);
    
    figure(5)
    plot(freq_prepost(:,1), freq_prepost(:,2), '+g')
    
    hold on
    plot(freq_postpre(:,1), freq_postpre(:,2), '+r')
    
    xlabel('Frequency')
    ylabel('STDP')
end

% Tools

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

   