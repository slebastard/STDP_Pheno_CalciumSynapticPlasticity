%% SIMULATION STARTER %%
% 1) Just defines the model, variables and integration scheme
% 2) Runs a single simulation and returns the full evolution of synaptic
% efficacy
% 3) Returns STDP curve by running simulation for several possible timings
% 4) Analyses the impact of frequency and of number of spike pairs on the
% plasticity during an experiment

% Simulation mode
% single    Step 2 only
% STDP      Step 3 only
% freq      Step 4 only
% all       All steps

mode = 'STDP';

%% 0) Define the environment

T = 200;
rho_0 = 0.3; % must be between 0 and 1

C_pre = 1;
C_post = 2;
tau_Ca = 20;
delay_pre = 13.7;
Ca_params = [C_pre, C_post, tau_Ca, delay_pre];

theta_dep = 1;
gamma_dep = 200;
dep_params = [theta_dep, gamma_dep];

theta_pot = 1.3;
gamma_pot = 321;
pot_params = [theta_pot, gamma_pot];

tau = 150000; % this is larger than I expected. Ask Brunel about this
noise_lvl = 0.0;

n_iter = 60;
frequency = 1;

naive_params = [T, rho_0, Ca_params, dep_params, pot_params, tau, noise_lvl];
model = 'naive';

%% 1) Define the stimulation history
d_t = 10;
pre_spikes_hist = linspace(0, 1000*(n_iter-1)/frequency, n_iter);
post_spikes_hist = pre_spikes_hist + d_t;
T = 1000*(n_iter-1)/frequency + 5*tau;
naive_params(1) = T;

%% 2) Full evolution of syn plast on a single simulation
if strcmp(mode, 'single') || strcmp(mode, 'all')
    int_scheme = 'euler_expl';
    scheme_step = 0.1;

    if strcmp(model, 'naive')
        [rho_hist, c_hist] = naive_model(pre_spikes_hist, post_spikes_hist, naive_params, int_scheme, scheme_step);
    end

    % Plotting rho as a function of time
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    t = linspace(0, T, T/scheme_step + 1);

    figure(1)
    plot(t, rho_hist);
    title('Evolution of synaptic strength, naive model');
    xlabel('Time');
    ylabel('Scaled synaptic strength');

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
end

%% 3) STDP curve

% Obtaining STDP curve
%%%%%%%%%%%%%%%%%%%%%%
if strcmp(mode, 'STDP') || strcmp(mode, 'all')
    t_min = -75;
    t_max = 75;
    dt = 3;

    stdp_params = [naive_params, t_min, t_max, dt, n_iter, frequency];
    STDP = get_STDP(model, stdp_params, int_scheme, scheme_step);

    figure(3)
    plot(STDP(:,1), STDP(:,2), '+');

    title('Plasticity as a function of pre-post spike delay')
    xlabel('Pre-post spike delay (ms)')
    ylabel('Relative change in synaptic strength')
    
    neutral_hline = refline([0 1]);
    neutral_hline.Color = 'b';
    
end

%% 4) Analysis of frequency and number of spike pairs
if strcmp(mode, 'freq') || strcmp(mode, 'all')
    dt = 10;

    freq_def = 100;
    n_iter_max = 200;

    n_iter_def = 100;
    freq_min = 1;
    freq_max = 50;

    freq_an_params = [naive_params, dt, freq_def, n_iter_max, n_iter_def, freq_min, freq_max];

    [
        dw_freq_prepost, ...
        dw_freq_postpre, ...
        dw_npairs_prepost, ...
        dw_npairs_postpre ...    
    ] ...
    = get_freq_an(model, freq_an_params, int_scheme, scheme_step);


    figure(4)

    plot(dw_freq_prepost(:,1), dw_freq_prepost(:,2), '+')
    hold on
    plot(dw_freq_postpre(:,1), dw_freq_postpre(:,2), '+')

    xlabel('Frequency (Hz)')
    ylabel('Relative change in synaptic weight')
    legend('dt = +10ms','dt = -10ms')
    title('Plasticity as a function of frequency (Hz), for 60 pairings at +-10ms')

    hold off


    figure(5)

    plot(dw_npairs_prepost(:,1), dw_npairs_prepost(:,2), '+')
    hold on
    plot(dw_npairs_postpre(:,1), dw_npairs_postpre(:,2), '+')

    xlabel('Nb of spike pairs at 1Hz');
    ylabel('Relative change in synaptic weight');
    legend('dt = +10ms','dt = -10ms')
    title('Plasticity as a function of the number of spike pairs, for pairings at +-10ms, at 1Hz')

    hold off
end