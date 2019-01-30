%% SIMULATION STARTER INSTRUCTIONS %%
% 0) Units
% All times in milliseconds (ms)
%
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
% - poisson   Provide response to correlated poisson processes
%
% - dataFit   Fit to experimental data
% - dataFitCuts   Same, but each frequency is output in a seperate
% subfigure
% -------------------------------------------------------------------------
%
% 2) Choose the model to use for simulation:
% -------------------------------------------------------------------------
% Models
% - naive   Graupner & Brunel 2012 with plasticity break (stopped below
%           calcium threshold, and noise
% - pheno   Graupner & Brunel 2012 with plasticity brake, including a
%           mapping of noise from  the analysis of G&B 2007 realistic model
% - caProd  Same model as pheno, but the calcium is the product of a
%           presynaptic and a postsynaptic component. This avoids
%           plasticity due to unilateral activity
% -------------------------------------------------------------------------
%
% 3) Depending on which modes you use, you will want to change the default
% dt, number of pairings and frequency in the environment section
% below.
%
% You should be all set!
clear all

env = getEnv();
addpath(genpath(env.functionsRoot), env.dataRoot);
data.path = strcat(env.dataRoot,'Venance2016/');            

simu.mode = 'STDP';
simu.model = 'caProd';

% Parameters controlling excitation history
simu.d_t = 2;
simu.n_iter = 10;
simu.frequency = 1;
simu.int_scheme = 'euler_expl';
simu.int_step = 0.5;

%% Environment definition

params = getSynapse();
params.T = 10000;
params.theta_act = params.theta_dep;
params.TD = 0;

params.tau_x = 1e3 .* params.tau_x;
params.tau_Ca = 1e3 .* params.tau_Ca;
params.tau_rho = 1e3 .* params.tau_rho;
params.delay_pre = 1e3 .* params.delay_pre;
params.tau_w = 50000;

prot = params;

N_A = 6.02e17; %mumol^(-1)
V = 2.5e-16; %L

% Defining default excitation timeline
pre_spikes_hist = linspace(0, 1000*(simu.n_iter-1)./simu.frequency, simu.n_iter);
post_spikes_hist = pre_spikes_hist + simu.d_t;
simu.T = max(1000*(simu.n_iter-1)./simu.frequency + abs(simu.d_t) + 10*params.tau_Ca);
post_spikes_hist = post_spikes_hist(post_spikes_hist > 0.48*simu.T & post_spikes_hist < 0.56*simu.T );

set(0,'DefaultFigureWindowStyle','docked')

%% mode='single') Full evolution of syn plast on a single simulation
if strcmp(simu.mode, 'single')
    if strcmp(simu.model, 'naive')
        [rho_hist, c_hist] = naive_model(pre_spikes_hist, post_spikes_hist, params, simu);
    elseif strcmp(simu.model, 'pheno')
        [rho_hist, w_hist, c_hist] = pheno_model(pre_spikes_hist, post_spikes_hist, params, simu);
    elseif strcmp(simu.model, 'caProd')
        [rho_hist, w_hist, c_hist] = caProd_model(pre_spikes_hist, post_spikes_hist, params, simu);
    end

    % Plotting rho as a function of time
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    t = linspace(0, simu.T, simu.T/simu.int_step + 1);

    figure('Name','SYN_SinglePhospho','NumberTitle','off')
    plot(t, rho_hist(1,1:length(t)));
    title('Evolution of average CaMKII state');
    xlabel('Time');
    ylabel('Average CaMKII state');

    % ToDo: add bumps of Ca as colored pins over x-axis

    figure('Name','SYN_SingleCa','NumberTitle','off')
    plot(t, c_hist(1:length(t),1));
    title('Evolution of calcium influx');
    xlabel('Time');
    ylabel('Calcium concentration');

    dep_thr = refline([0 params.theta_dep]);
    dep_thr.Color = 'r';

    pot_thr = refline([0 params.theta_pot]);
    pot_thr.Color = 'g';
    
    act_thr = refline([0 params.theta_act]);
    act_thr.Color = 'm';
    
    if strcmp(simu.model, 'pheno') || strcmp(simu.model, 'caProd')
        figure('Name','SYN_SingleWght','NumberTitle','off')
        plot(t, w_hist(1:length(t),1));
        title('Evolution of synaptic strength')
        xlabel('Time');
        ylabel('Average synaptic strength');
    end
end

%% mode='STDP') Provides a STDP curve for provided frequency and number of pairings

% Obtaining STDP curve
%%%%%%%%%%%%%%%%%%%%%%
if strcmp(simu.mode, 'STDP')
    
    STDP = simu;
    STDP.dt.min = -75;
    STDP.dt.max = 75;
    STDP.dt.step = 1;
    STDP.mode = 'rel';

    if strcmp(simu.model, 'caProd')
        STDP.function = get_STDP_CaProd(STDP, params);
    else
        STDP.function = get_STDP(STDP, params);
    end
    % Validation against simulation
    % [STDP.function, STDP.func_sim] = get_both_STDP_CaProd(STDP, params);

    STDP.integral = sum(STDP.function(:,2))*STDP.dt.step;
    STDP.expectation = STDP.integral/(STDP.dt.max - STDP.dt.min);
    
    
    figure('Name','SYN_STDP','NumberTitle','off')
    STDP.plot = plot(STDP.function(:,1), STDP.function(:,2), '.b');
%     hold on
%     plot(STDP.func_sim(:,1), STDP.func_sim(:,2), 'xg');
    
    dim = [.15 .6 .3 .3];
    sidestr = strcat('STDP expectation: ', num2str(100*STDP.expectation, 3+1), '%');
    annotation('textbox',dim,'String',sidestr,'FitBoxToText','on');

    neutral_hline = refline([0 1]);
    neutral_hline.Color = 'b';
    
    title('Plasticity as a function of pre-post spike delay')
    xlabel('Pre-post spike delay (ms)')
    ylabel('Relative change in synaptic strength')
    
end



%% mode='freq3') STDP = f(freq,dt) 3D plot

if strcmp(simu.mode, 'freq3')
    freq3Map = simu;
    freq3Map.mode = 'rel';
    
    freq3Map.dt.min = -100;
    freq3Map.dt.max = 100;
    freq3Map.dt.step = 4;   
    freq3Map.freq.max = 50;
    freq3Map.freq.step = 5;

    freq3Map.map = get_freq_heatmap(freq3Map, params);
    
    n_freq = 1 + floor((freq3Map.freq.max-1)/freq3Map.freq.step);
    n_dt = 1 + floor((freq3Map.dt.max-freq3Map.dt.min)/freq3Map.dt.step);

    [freq_grid, dt_grid] = meshgrid(1:freq3Map.freq.step:freq3Map.freq.max, freq3Map.dt.min:freq3Map.dt.step:freq3Map.dt.max);
    freq3Map.interpol = griddata(freq3Map.map(:,1), freq3Map.map(:,2), freq3Map.map(:,3), freq_grid, dt_grid);
    % ribboncoloredZ(gca,dt_grid,dataFit.interpol);
    
%     freq3Map.heat = zeros(n_freq, n_dt);
%     freq3Map.heat(sub2ind([n_freq,n_dt], repelem(1:n_freq,1,n_dt), repmat(1:n_dt,1,n_freq))) = freq3Map.map(:,3);
    
    figure('Name','SYN_FreqMap3','NumberTitle','off')
    surf(freq_grid, dt_grid, freq3Map.interpol);
    colormap(bluewhitered), colorbar;
    alpha 0.3
    
%     scatter3(freq3Map.map(:,1),freq3Map.map(:,2),freq3Map.map(:,3), '.');
    
%     freq3Map.plot = imagesc(freq3Map.heat');
%     colormap('hot');
%     colorbar;
%     
%     xlabels = 1 + freq_step.*(xticks-1);
%     ylabels = dtmin + step_dt.*(yticks-1);
%     
%     xtickformat('%.1f')
%     set(gca, 'XTickLabel', xlabels);
%     set(gca, 'YTickLabel', ylabels);
% 
%     title('Relative change in syn plast as a function of frequency and dt');
%     xlabel('Frequency');
%     ylabel('dt');
end

%% mode='freq') STDP = f(freq) 3D plot for two opposite timings
if strcmp(simu.mode, 'freq')
    freqMap.freq.max = 75;
    freqMap.freq.step = 0.5;
    freqMap.int_scheme = simu.int_scheme;
    freqMap.d_t = simu.d_t;
    freqMap.max_freq = simu.max_freq;
    freqMap.step_freq = simu.step_freq;

    [freqMap.map.prepost, freqMap.map.postpre] = get_freqSTDP(simu.model, 'rel', freqMap);
    
    figure('Name','SYN_Freq','NumberTitle','off')
    plot(freq_prepost(:,1), freq_prepost(:,2), '+g')
    
    hold on
    plot(freq_postpre(:,1), freq_postpre(:,2), '+r')
    
    xlabel 'Frequency';
    ylabel 'STDP';
end

%% mode='pairs3') STDP = f(n_iter,dt) 3D plot

if strcmp(simu.mode, 'pairs3')
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
    
    figure('Name','SYN_Pairs3','NumberTitle','off')
    
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
if strcmp(simu.mode, 'pairs')
    max_pairs = 200;
    step_pairs = 1;
    
    pairs_params = [model_params, frequency];
    [numpairs_prepost, numpairs_postpre] = get_pairsSTDP(model, 'rel', pairs_params, int_scheme, d_t, max_pairs, step_pairs);
    
    figure('Name','SYN_Pairs','NumberTitle','off')
    plot(numpairs_prepost(:,1), numpairs_prepost(:,2), '+g')
    xlabel 'Number of pairings';
    ylabel 'STDP';
    hold on
    plot(numpairs_postpre(:,1), numpairs_postpre(:,2), '+r')

end



%% mode='dataFit') Fitting to data from Venance lab

% Comparing the model to STDP=f(freq,dt) data from L. Venance (INSERM)
data.freqSTDP.data = csvread(strcat(data.path, 'STDP_Frequency.csv'),1,0);
data.freqSTDP.length = size(data.freqSTDP.data,1);

if strcmp(simu.mode, 'dataFit')
    dataFit = simu;
    dataFit.mode = 'rel';

    dataFit.dt.min = -100;
    dataFit.dt.max = 100;
    dataFit.dt.step = 4;

    dataFit.freq.max = 200;
    dataFit.freq.step = 5;

    dataFit.heat = get_freq_heatmap(dataFit, params);

    figure('Name','SYN_DataFit','NumberTitle','off')
    data.freqSTDP.freqs=unique(data.freqSTDP.data(:,5));
    n_data_freqs=length(data.freqSTDP.freqs); 

    % scatter3(data.freqSTDP.data(:,5), data.freqSTDP.data(:,2), data.freqSTDP.data(:,3)./100, 50*ones(size(data.freqSTDP.data,1),1), '*r')
    % hold on
    
    [freq_grid, dt_grid] = meshgrid(1:dataFit.freq.step:dataFit.freq.max, dataFit.dt.min:dataFit.dt.step:dataFit.dt.max);
    dataFit.interpol = griddata(dataFit.heat(:,1), dataFit.heat(:,2), dataFit.heat(:,3), freq_grid, dt_grid);
    % ribboncoloredZ(gca,dt_grid,dataFit.interpol);
    surf(freq_grid, dt_grid, dataFit.interpol);
    colormap(bluewhitered), colorbar;
    alpha 0.3
    
%     for f=1:n_data_freqs
%         hold on
%         ids=find(data.freqSTDP.data(:,5)==data.freqSTDP.freqs(f) & data.freqSTDP.data(:,7)~=0);
%         filtered_freq=data.freqSTDP.data(ids,:);
%         [a,b]=sort(filtered_freq(:,2));
%         h = ribbon(filtered_freq(b,2), filtered_freq(b,3)./100, 0.15);
%         set(h, 'XData', filtered_freq(b,5)-1 + get(h, 'XData'));
%     end   
    
%     dataFit.1Hz = freq_data(floor(freq_data(:,5))==1,:);
%     
%     stdp_params = [model_params, dataFit.dt.min, dataFit.dt.max, dataFit.dt.step, 100, 1];
%     dataFit.STDP = get_STDP(model, 'rel', stdp_params, int_scheme, scheme_step);
%     plot(dataFit.STDP(:,1), dataFit.STDP(:,2), '.b')
%     hold on
%     plot(dataFit.1Hz(:,2),dataFit.1Hz(:,3)./100,'xr')
%     neutral_hline = refline([0 1]);
%     neutral_hline.Color = 'b';   
%     ('Plasticity as a function of pre-post spike delay');
%     xlabel('Pre-post spike delay (ms)');
%     ylabel('Relative change in synaptic strength');
    dataFit.paramPos.x = 10.0;
    dataFit.paramPos.y = 65.0;
    dataFit.paramPos.z = 2.6;
    dataFit.paramPos.colSepY = 25.0;
    dataFit.paramPos.colSepZ = 0.15;
    dataFit.paramPos.ftSize = 8;
    stampParams(params, dataFit.paramPos);
    
end

if strcmp(simu.mode, 'dataFitCuts')
    dataFit = simu;
    dataFit.mode = 'rel';

    dataFit.dt.min = -100;
    dataFit.dt.max = 100;
    dataFit.dt.step = 4;

    dataFit.freq.max = 10;
    dataFit.freq.step = 1;

    
    data.freqSTDP.freqs=unique(data.freqSTDP.data(:,5));
    n_data_freqs=length(data.freqSTDP.freqs); 

    dataFit.paramPos.x = 10.0;
    dataFit.paramPos.y = 65.0;
    dataFit.paramPos.z = 2.6;
    dataFit.paramPos.colSepY = 25.0;
    dataFit.paramPos.colSepZ = 0.15;
    dataFit.paramPos.ftSize = 8;
    stampParams(params, dataFit.paramPos);
    
    for f=1:n_data_freqs
        ids=find(data.freqSTDP.data(:,5)==data.freqSTDP.freqs(f) & data.freqSTDP.data(:,7)~=0);
        filtered_freq=data.freqSTDP.data(ids,:);
        [a,b]=sort(filtered_freq(:,2));
        figure('Name',strcat('SYN_DataFitCuts',f),'NumberTitle','off')
        plot(filtered_freq(b,2), filtered_freq(b,3)./100, 's', 'MarkerSize', 8)
        hold on
        dataFit.frequency = data.freqSTDP.freqs(f);
        STDP.function = get_STDP_CaProd(dataFit, params);
        plot(STDP.function(:,1), STDP.function(:,2), 'r', 'LineWidth', 2)
        xlabel('Time delay (ms)')
        ylabel('Relative change in synaptic strength')
        subtitle(strcat('Datafit at ', num2str(data.freqSTDP.freqs(f),1), 'Hz'))
    end   
    
end

%% mode='poisonSingle') Simulate response to Poisson processes
if strcmp(simu.mode, 'poissonSingle')
    % Generating Correlation matrix
    nu_pre = 1.0;       % Hz
    nu_post = 1.0;      % Hz
    c12 = 50;
    C = [c12*1.05*nu_pre/nu_post c12; c12 c12*1.05*nu_post/nu_pre];
    simu.T = 2000;      % ms
    tc = 50;
    
%     % The following simulates correlated Poisson for exponential
%     % correlation functions only (see Brette 2008, section 3)
%     [t, I] = corrPoisson( 2, [nu_pre; nu_post], C, simu.T, tc);
%     
%     preIds = find(I(1,:));
%     pre_spikes_hist = t(preIds);
%     
%     postIds = find(I(2,:));
%     post_spikes_hist = t(postIds);

    % The following simulates two independent Poisson processes
    % Rates are cast back to s^(-1)
    t = indPoisson( 2, [1000/nu_pre; 1000/nu_post], simu.T);
    pre_spikes_hist = t(1,:);
    post_spikes_hist = t(2,:);

    if strcmp(simu.model, 'naive')
        [rho_hist, c_hist] = naive_model(pre_spikes_hist, post_spikes_hist, params, simu);
    elseif strcmp(simu.model, 'pheno')
        [rho_hist, w_hist, c_hist] = pheno_model(pre_spikes_hist, post_spikes_hist, params, simu);
    elseif strcmp(simu.model, 'caProd')
        [rho_hist, w_hist, c_hist] = caProd_model(pre_spikes_hist, post_spikes_hist, params, simu);
    end

    % Plotting rho as a function of time
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    t = linspace(0, simu.T, simu.T/simu.int_step + 1);

%     figure('Name','SYN_CrossCorrPho','NumberTitle','off')
%     plot(t, rho_hist);
%     title('Evolution of average CaMKII state');
%     xlabel('Time');
%     ylabel('Average CaMKII state');
% 
%     % ToDo: add bumps of Ca as colored pins over x-axis
% 
%     figure('Name','SYN_CrossCorrCa','NumberTitle','off')
%     plot(t, c_hist);
%     title('Evolution of calcium influx');
%     xlabel('Time');
%     ylabel('Calcium concentration');
% 
%     dep_thr = refline([0 params.theta_dep]);
%     dep_thr.Color = 'r';
% 
%     pot_thr = refline([0 params.theta_pot]);
%     pot_thr.Color = 'g';
%     
%     act_thr = refline([0 params.theta_act]);
%     act_thr.Color = 'm';
%     
%     if strcmp(simu.model, 'pheno') || strcmp(simu.model, 'caProd')
%         figure('Name','SYN_CrossCorrW','NumberTitle','off')
%         plot(t, w_hist);
%         title('Evolution of synaptic strength')
%         xlabel('Time');
%         ylabel('Average synaptic strength');
%     end
    
    figure('Name','SYN_CrossCorr','NumberTitle','off')
    subplot(2,2,1)
    [rPre, lagsPre] = xcorr(pre_spikes_hist, pre_spikes_hist, 0.2*simu.T);
    plot(lagsPre, rPre)
    title('Pre auto-correlation')
    
    subplot(2,2,2)
    [rX, lagsX] = xcorr(pre_spikes_hist, post_spikes_hist, 0.2*simu.T);
    plot(lagsX, rX)
    title('Pre-post cross-correlation')
    
    subplot(2,2,3)
    [rX, lagsX] = xcorr(pre_spikes_hist, post_spikes_hist, 0.2*simu.T);
    plot(lagsX, rX)
    title('Pre-post cross-correlation')
    
    subplot(2,2,4)
    [rPost, lagsPost] = xcorr(post_spikes_hist, post_spikes_hist, 0.2*simu.T);
    plot(lagsPost, rPost)
    title('Post auto-correlation')
end

%% mode='poisonMap') STDP map for Poisson processes
if strcmp(simu.mode, 'poissonMap')
    % Generating Correlation matrix
    pSTDP = simu;
    pSTDP.T = 1000;
    pSTDP.nuPre.min = 1;
    pSTDP.nuPre.max = 200;
    pSTDP.nuPre.step = 15;
    pSTDP.nuPost.min = 1;
    pSTDP.nuPost.max = 200;
    pSTDP.nuPost.step = 15;
    pSTDP.nTry = 1;
    
    pSTDP.corr.type = 'none';
    pSTDP.corr.c12 = 50;
    pSTDP.corr.tc = 50;
    
    pSTDP.map = poissonMap(params, pSTDP);
    
    % Plotting the STDP surface obtained
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    [pre_grid_intp, post_grid_intp] = meshgrid(pSTDP.nuPre.min:0.2*pSTDP.nuPre.step:pSTDP.nuPre.max, pSTDP.nuPost.min:0.2*pSTDP.nuPost.step:pSTDP.nuPost.max);
    [pre_grid_spl, post_grid_spl] = meshgrid(pSTDP.nuPre.min:pSTDP.nuPre.step:pSTDP.nuPre.max, pSTDP.nuPost.min:pSTDP.nuPost.step:pSTDP.nuPost.max);
    pSTDP.interpol = griddata(pre_grid_spl, post_grid_spl, pSTDP.map, pre_grid_intp, post_grid_intp);
    % ribboncoloredZ(gca,dt_grid,dataFit.interpol);
    figure('Name','SYN_PoissonMap','NumberTitle','off')
    pSTDP.plot = log(pSTDP.interpol);
    surf(pre_grid_intp, post_grid_intp, pSTDP.plot);
    xlabel('Presyn rate')
    ylabel('Postsyn rate')
    zlabel('Log potentiation')
    colormap(bluewhitered), colorbar;
    alpha 0.3
    
    pSTDP.paramPos.x = pSTDP.nuPre.min + 0.85*(pSTDP.nuPre.max - pSTDP.nuPre.min);
    pSTDP.paramPos.y = pSTDP.nuPost.min + 0.85*(pSTDP.nuPost.max - pSTDP.nuPost.min);
    pSTDP.paramPos.z = max(max(pSTDP.plot)) - 0.6*(max(max(pSTDP.plot)) - min(min(pSTDP.plot)));
    pSTDP.paramPos.colSepY = 0.15*(pSTDP.nuPost.max - pSTDP.nuPost.min);
    pSTDP.paramPos.colSepZ = 0.05*(max(max(pSTDP.plot)) - min(min(pSTDP.plot)));
    pSTDP.paramPos.ftSize = 8;
    stampParams(params, pSTDP.paramPos);
end

%% mode='poisonCut') STDP map for Poisson processes
if strcmp(simu.mode, 'poissonCut')
    % Generating Correlation matrix
    pSTDP = simu;
    pSTDP.dir = 'pre';
    pSTDP.val = [1 5 20 40];
    pSTDP.T = 2000;
    pSTDP.nu.min = 1;
    pSTDP.nu.max = 40;
    pSTDP.nu.step = 0.5;
    pSTDP.nTry = 20;

    pSTDP.corr.type = 'none';
    pSTDP.corr.c12 = 50;
    pSTDP.corr.tc = 50;
    
%     % The following simulates correlated Poisson for exponential
%     % correlation functions only (see Brette 2008, section 3)
%     [t, I] = corrPoisson( 2, [nu_pre; nu_post], C, T);
%     
%     preIds = find(I(1,:));
%     pre_spikes_hist = t(preIds);
%     
%     postIds = find(I(2,:));
%     post_spikes_hist = t(postIds);

    % The following simulates two independent Poisson processes
    
    pSTDP.cut = poissonCut(params, pSTDP);

    % Plotting the STDP surface obtained
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    figure('Name','SYN_PoissonCut','NumberTitle','off')
    errorbar(pSTDP.cut.mean(:,1),pSTDP.cut.mean(:,2),pSTDP.cut.var(:,2));
    if strcmp(pSTDP.dir, 'pre')
        title(strcat('Relative plasticity as a function of postsyn rate, at nuPre=', pSTDP.val, 'Hz'))
        xlabel('Postsyn Poisson rate')
    else
        title(strcat('Relative plasticity as a function of postsyn rate, at nuPost', pSTDP.val, 'Hz'))
        xlabel('Presyn Poisson rate')
    end
    ylabel('Rel STDP')
    alpha 0.3

end


%% Find STDP curve attributes as a function of parameters
if strcmp(simu.mode, 'STDP_SWEEP')
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
    default = gamma_dep;
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