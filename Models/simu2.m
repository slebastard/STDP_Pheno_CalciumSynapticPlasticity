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
% dt, number of pairings and frequency in the environment section
% below.
%
% You should be all set!

simu.mode = 'STDP_SWEEP';
simu.model = 'pheno';

% Parameters controlling excitation history
simu.d_t = 30;
simu.n_iter = 100;
simu.frequency = 1;

%% Environment definition

params.T = 10;
params.rho_max = 200;
params.S_attr = 40;

params.C_pre = 0.4;
params.C_post = 0.84;
params.tau_Ca = 80;
params.delay_pre = -15;
params.Ca_params = [C_pre, C_post, tau_Ca, delay_pre];

params.theta_dep = 1;
params.gamma_dep = 200;
params.dep_params = [theta_dep, gamma_dep];

params.theta_pot = 1.08;
params.gamma_pot = 120;
params.pot_params = [theta_pot, gamma_pot];

N_A = 6.02e17; %mumol^(-1)
V = 2.5e-16; %L
params.theta_act = theta_dep;

params.tau_rho = 100000;
params.tau_w = 500000;

params.noise_lvl = 25; %1/sqrt(N_A*V);

% Initialization of variables
params.rho_0 = 25; % must be between 0 and rho_max
params.w_0 = transfer(rho_0, S_attr, noise_lvl);

simu.int_scheme = 'euler_expl';
simu.scheme_step = 0.5;

datapath = 'Data/Venance2016/';

% Defining default excitation timeline
pre_spikes_hist = linspace(0, 1000*(simu.n_iter-1)./simu.requency, simu.n_iter);
post_spikes_hist = pre_spikes_hist + simu.d_t;
simu.T = max(1000*(simu.n_iter-1)./simu.frequency + abs(simu.d_t) + 10*params.tau_Ca);

%% mode='single') Full evolution of syn plast on a single simulation
if strcmp(mode, 'single') || strcmp(mode, 'all')
    if strcmp(model, 'naive')
        [rho_hist, c_hist] = naive_model(pre_spikes_hist, post_spikes_hist, params, simu);
    elseif strcmp(model, 'pheno')
        [rho_hist, w_hist, c_hist] = pheno_model(pre_spikes_hist, post_spikes_hist, params, simu);
    end

    % Plotting rho as a function of time
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    t = linspace(0, simu.T, simu.T/simu.scheme_step + 1);

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

    dep_thr = refline([0 params.theta_dep]);
    dep_thr.Color = 'r';

    pot_thr = refline([0 params.theta_pot]);
    pot_thr.Color = 'g';
    
    act_thr = refline([0 params.theta_act]);
    act_thr.Color = 'm';
    
    if strcmp(model, 'pheno')
        figure(3)
        plot(t, w_hist);
        title('Evolution of synaptic strength')
        xlabel('Time');
        ylabel('Average synaptic strength');
    end
end

%% mode='STDP') Provides a STDP curve for provided frequency and number of pairings

% Obtaining STDP curve
%%%%%%%%%%%%%%%%%%%%%%
if strcmp(mode, 'STDP') || strcmp(mode, 'all')
    STDP.dt.min = -100;
    STDP.dt.max = 100;
    STDP.dt.step = 2;
    STDP.n_iter = simu.n_iter;
    STDP.frequency = simu.frequency;

    STDP.function = get_STDP(model, 'rel', params, STDP);
    % Validation against simulation
    % [STDP_an, STDP_sim] = get_both_STDP(model, 'rel', stdp_params, int_scheme, scheme_step);

    figure(4)
    STDP.plot = plot(STDP.function(:,1), STDP.function(:,2), '.b');
    % plot(STDP_an(:,1), STDP_an(:,2), '+r');
    % hold on
    % plot(STDP_sim(:,1), STDP_sim(:,2), 'xg');
   

    neutral_hline = refline([0 1]);
    neutral_hline.Color = 'b';
    
    title('Plasticity as a function of pre-post spike delay')
    xlabel('Pre-post spike delay (ms)')
    ylabel('Relative change in synaptic strength')
    
end

%% mode='freq3') STDP = f(freq,dt) 3D plot

if strcmp(mode, 'freq3')
    freq3Map.dt.min = -75;
    freq3Map.dt.max = 75;
    freq3Map.dt.step = 1;
    
    freq3Map.freq.max = 40;
    freq3Map.freq.step = 0.5;
    freq3Map.n_iter = simu.n_iter;
    freq3Map.int_scheme = simu.int_scheme;
    freq3Map.dt_params = simu.dt_params;
    freq3Map.freq_params = simu.freq_params;
    
    freq3Map.map = get_freq_heatmap(model, 'rel', freq3Map);
    
    n_freq = 1 + floor((freq_max-1)/freq_step);
    n_dt = 1 + floor((dtmax-dtmin)/step_dt);
    
    freq3Map.heat = zeros(n_freq, n_dt);
    freq3Map.heat(sub2ind([n_freq,n_dt], repelem(1:n_freq,1,n_dt), repmat(1:n_dt,1,n_freq))) = frq_map(:,3);
    
    figure(5)
    
%     scatter3(frq_map(:,1),frq_map(:,2),frq_map(:,3), '.');
    
    freq3Map.plot = imagesc(freq3Map.heat');
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
    freqMap.freq.max = 75;
    freqMap.freq.step = 0.5;
    freqMap.int_scheme = simu.int_scheme;
    freqMap.d_t = simu.d_t;
    freqMap.max_freq = simu.max_freq;
    freqMap.step_freq = simu.step_freq;

    [freqMap.map.prepost, freqMap.map.postpre] = get_freqSTDP(simu.model, 'rel', freqMap);
    
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
freq_data = csvread(strcat(datapath, 'STDP_Frequency.csv'),1,0);
n_data = size(freq_data,1);

if strcmp(mode, 'dataFit')
    dataFit.dt.min = -100;
    dataFit.dt.max = 100;
    dataFit.dt.step = 2;

    dataFit.freq.max = 10;
    dataFit.freq.step = 0.5;
    dataFit.n_iter = simu.n_iter;

    dataFit.heat = get_freq_heatmap(model, 'rel', dataFit);
    
    scatter3(freq_data(:,5), freq_data(:,2), freq_data(:,3)./100, 50*ones(size(freq_data,1),1), '*r')
    hold on
    
    [freq_grid, dt_grid] = meshgrid(1:freq_step:freq_max, dtmin:step_dt:dtmax);
    dataFit.interpol = griddata(dataFit.heat(:,1), dataFit.heat(:,2), dataFit.heat(:,3), freq_grid, dt_grid);
    surf(freq_grid, dt_grid, dataFit.interpol);
    alpha 0.3
    
%     data1Hz = freq_data(floor(freq_data(:,5))==1,:);
%     
%     stdp_params = [model_params, dtmin, dtmax, step_dt, 100, 1];
%     STDP = get_STDP(model, 'rel', stdp_params, int_scheme, scheme_step);
%     plot(STDP(:,1), STDP(:,2), '.b')
%     hold on
%     plot(data1Hz(:,2),data1Hz(:,3)./100,'xr')
%     neutral_hline = refline([0 1]);
%     neutral_hline.Color = 'b';   
%     ('Plasticity as a function of pre-post spike delay');
%     xlabel('Pre-post spike delay (ms)');
%     ylabel('Relative change in synaptic strength');
    
end


%% Find STDP curve attributes as a function of parameters
if strcmp(mode, 'STDP_SWEEP')
    sweep.caBump.min = 0.5;
    sweep.caBump.max = 3;
    sweep.caBump.step = 0.1;
    sweep.caBump.n = 1 + floor((sweep.caBump.max - sweep.caBump.min)/sweep.caBump.step);

    sweep.dt.min = -100;
    sweep.dt.max = 100;
    sweep.dt.step = 2;
    sweep.dt.n = 1 + floor((sweep.dt.max - sweep.dt.min)/sweep.dt.step);

    sweep.STDP.curves = zeros(sweep.caBump.n, sweep.caBump.n, sweep.dt.n);
    sweep.STDP.int = zeros(sweep.caBump.n, sweep.caBump.n);

    sweep.save.curves = 0;
    sweep.save.int = 1;
    sweep.save.path = 'Models/Outputs/Pheno/STDP/paramSweep';

    stdp_params = [model_params, sweep.dt.min, sweep.dt.max, sweep.dt.step, 1, 1];

    for preID=1:sweep.caBump.n
        pre = sweep.caBump.min + preID*sweep.caBump.step;
        for postID=1:preID
            post = sweep.caBump.min + postID*sweep.caBump.step;
            
            model_params(4) = pre;
            model_params(5) = post;
            stdp_params = [model_params, sweep.dt.min, sweep.dt.max, sweep.dt.step, 1, 1];
            
            STDP = get_STDP(model, 'rel', stdp_params);
            sweep.STDP.curves(preID, postID, :) = STDP(:,2);
            sweep.STDP.int(preID, postID) = sweep.dt.step * sum(STDP(:,2));
        end
    end

    if sweep.save.curves
        csvwrite(strcat(sweep.save.path,'_curves.dat'), sweep.STDP.curves);
    end

    if sweep.save.int
        csvwrite(strcat(sweep.save.path,'_int.dat'), sweep.STDP.int);
    end
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