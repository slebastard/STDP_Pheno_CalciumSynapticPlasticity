function out = simuNetwork(syn, simu, net, neu, init, plt, gif)
% Simulation of the Brunel JCNS 2000 network with electrode recording
% spikes (for paper Destexhe Touboul on Criticality).

% Modified by S. Lebastard, 2018
% Should accomodate mutliple models of plasticity, including none
% Should allow for heterogeneous synapse rules, ie manipulation of matrix
% of parameters

% ToDo
% - Allow user to start either with uniform weights or with random weights
% - Create clusterization index and look at its evolution
%   For this purpose, use techniques from spectral clustering

% ALL POTENTIALS IN mV
% ALL TIMES IN s
% Add plot tools to path
env = getEnv();
addpath(genpath(env.functionsRoot));


%% Parameterization

% %%%%%%%%%% LOADING PHASE INFORMATION %%%%%%%%

seed = 19;
rng(seed);

if isempty(simu.phases)
    return
else
    simu.T=0;
    simu.nIterTot = 0;
    nPhase = length(simu.phases);
    simu.nIter = zeros(1,nPhase);
    for phID=1:nPhase
        simu.T = simu.T + simu.phases{phID}.T;
        simu.phases{phID}.firstIter = simu.nIterTot + 1;
        simu.phases{phID}.nIter = simu.phases{phID}.T / simu.phases{phID}.dt;
        simu.phases{phID}.lastIter = simu.nIterTot + simu.phases{phID}.nIter;
        simu.nIterTot = simu.nIterTot + floor(simu.phases{phID}.nIter);
        
        % Build PlastOn and PlastInOn matrices
        simu.phases{phID}.PlastON = zeros(net.N,net.N);
        simu.phases{phID}.PlastON(1:net.NE,1:net.NE) = (simu.phases{phID}.EE == 1);
        simu.phases{phID}.PlastON(net.NE+1:net.N,1:net.NE) = (simu.phases{phID}.IE == 1);
        simu.phases{phID}.PlastON(1:net.NE,net.NE+1:net.N) = (simu.phases{phID}.EI == 1);
        simu.phases{phID}.PlastON(net.NE+1:net.N,net.NE+1:net.N) = (simu.phases{phID}.II == 1);
        
        simu.phases{phID}.PlastInON = zeros(net.N,net.NIn);
        simu.phases{phID}.PlastInON(1:net.NE,:) = (simu.phases{phID}.InE == 1);
        simu.phases{phID}.PlastInON(net.NE+1:net.N,:) = (simu.phases{phID}.InE == 1);
    end
    if simu.T == 0
        return
    end
end

% %%%%%%%%   STATIC PARAMETERS  %%%%%%%% 

simu.dt = simu.phases{1}.dt;

% Protocol
prot.TD = 0;
prot.S_attr = syn.S_attr;
prot.noise_lvl = syn.noise_lvl;

% Plotting
plt.timeSpl.inter = 70;
gif.minN = 15;
gif.maxN = 40;
gif.ratesHist = 0;
gif.updIter = 80;
% Histograms
plt.nbins = 100;
plt.edgesRho = linspace(0,syn.rho_max,plt.nbins+1);
plt.edgesW_exc = linspace(0,1,plt.nbins+1);
plt.edgesW_inh = linspace(-net.g,0,plt.nbins+1);

% Network
net.NI=net.N-net.NE;             % # excitatory neurons
CE=round(net.Connectivity*net.NE);          % Number of exc connections
CI=round(net.Connectivity*net.NI);          % Number of inh connections
net.C_ext=CE;

% Neurons
neu.tau=20e-3;
neu.t_rp=2e-3;

% Input
net.nu_thresh=(neu.V_t-neu.V_r)/(syn.J*CE*neu.tau);   % Frequency needed for a neuron to reach the threshold. 
net.nu_ext=net.rExtRel*net.nu_thresh; %./net.NIn;       % external Poisson input rate


%% Creating network

% %%%%%% INITIALIZING SYNAPTIC QTIES %%%%%%

% %%% Synaptic weights %%%
net.W = zeros(net.N,net.N);   % 1st index for postsynaptic neuron, 2nd for presynaptic
net.NIn = net.N;
net.WIn = zeros(net.N,net.NIn);


for i=1:net.N
    ECells=randperm(net.NE);
    ECells=ECells(1:CE);
    net.W(i,ECells) = strcmp(init.mode, 'rand').*rand(1,CE) + ~strcmp(init.mode, 'rand')*0.5*ones(1,CE);
    
    ICells=randperm(net.NI);
    ICells=net.NE+ICells(1:CI);
    net.W(i,ICells)= -strcmp(init.mode, 'rand').*rand(1,CI) - ~strcmp(init.mode, 'rand')*0.5*ones(1,CI);
    
    %InputCells=randperm(net.NIn);
    %InputCells=InputCells(1:ceil(1*net.NIn));
    %net.WIn(i,InputCells) = init.c*0.5;
    net.WIn(i,i) = 0.5;
end

net.synSign = 2*syn.J.*(net.W>0) - 2*net.g.*syn.J.*(net.W<0);
%net.W(net.W<0) = -net.W(net.W<0); % Inhibitory synapses are turned positive, but we have track of who is inhibitory through net.synSign
excNeurons = (1:net.NE)'; inhNeurons = (net.NE+1:net.N)';

% Mean stats for phase plot
net.meanWexc = mean(2*syn.J.*net.W(net.W>0)).*ones(simu.nIterTot+1,1);
net.meanWinh = mean(2*net.g.*syn.J.*net.W(net.W<0)).*ones(simu.nIterTot+1,1);

rateEst = zeros(net.N, simu.nIterTot);
histRates = zeros(plt.nbins, simu.nIterTot);

plt.patchTime = 0.05;
plt.maxPatchSize = 2*floor(0.5*plt.patchTime/simu.dt);
plt.ISI = zeros(simu.nIterTot,2); 
%plt.fourier = zeros(net.N, simu.nIterTot-plt.maxPatchSize);
plt.phaseMetr = zeros(simu.nIterTot,2);
plt.freqMetr = zeros(simu.nIterTot,1);
plt.regMetr = zeros(simu.nIterTot,1);
plt.synMetr = zeros(simu.nIterTot,1);
plt.mainModeFreq = zeros(simu.nIterTot,1);
plt.regime = zeros(simu.nIterTot,2);

subIDsExc = excNeurons( ...
    randperm(size(excNeurons,1),floor(0.8*min(max(gif.minN, 0.3*net.N), gif.maxN))) ...
    );
subIDsInh = inhNeurons( ...
    randperm(size(inhNeurons,1),ceil(0.2*min(max(gif.minN, 0.3*net.N), gif.maxN))) ...
    );
G = digraph(net.W([subIDsExc;subIDsInh]',[subIDsExc;subIDsInh]')');    % Graph for plotting evolution of network
A = abs(net.W + net.W'); % Weight matrix for spectral clustering
% Symmetrization, all information on directionnality is lost
plt.histRho = zeros(plt.nbins, simu.nIterTot);
plt.histW_exc = zeros(plt.nbins, simu.nIterTot);
plt.histW_inh = zeros(plt.nbins, simu.nIterTot);

% %%% Recurrent variables %%%
net.ca = zeros(net.N,net.N);
net.xpre = ones(net.N,net.N);
net.xpost = ones(net.N,net.N);
net.rho = transferinv(abs(net.W), syn.S_attr, syn.noise_lvl, syn.rho_max);
net.W = net.synSign.*transfer(net.rho, prot);
net.actPot = zeros(net.N,net.N);
net.actDep = zeros(net.N,net.N);

% %%% Input variables %%%
net.caIn = zeros(net.N,net.NIn);
net.xpreIn = ones(net.N,net.NIn);
net.xpostIn = ones(net.N,net.NIn);
net.rhoIn = transferinv(net.WIn, syn.S_attr, syn.noise_lvl, syn.rho_max);
%net.WIn = 2.*syn.J.*transfer(net.rhoIn, prot)./net.NIn;
net.WIn = 2.*syn.J.*transfer(net.rhoIn, prot);
net.actPotIn = zeros(net.N,net.NIn);
net.actDepIn = zeros(net.N,net.NIn);

% Sample synapses
splSyn.n = 15;
[synOut, synIn] = find(net.W);
perm = randperm(length(synOut));
rpOut = synOut(perm);
rpIn = synIn(perm);

splSyn.PostNeurons = rpOut(1:splSyn.n);
splSyn.PreNeurons = rpIn(1:splSyn.n);
splSyn.IDs = sub2ind(size(net.ca), splSyn.PostNeurons, splSyn.PreNeurons);
splSyn.ca = zeros(splSyn.n,simu.nIterTot);
splSyn.xpre = zeros(splSyn.n,simu.nIterTot);
splSyn.xpost = zeros(splSyn.n,simu.nIterTot);
splSyn.rho = zeros(splSyn.n,simu.nIterTot);
splSyn.w = zeros(splSyn.n,simu.nIterTot);

% Sample neurons
splNeu.NE = 3;
splNeu.NI = 2;
splNeu.ExcNeurons = excNeurons(randperm(size(excNeurons,1),splNeu.NE));
splNeu.InhNeurons = inhNeurons(randperm(size(inhNeurons,1),splNeu.NI));

plt.regimeThr = zeros(1,2);
plt.regimeThr(1,1) = 0.27;
plt.regimeThr(1,2) = 0.06;

%% Initialization
% %%%%%%%%   SIMULATION PARAMETERS %%%%%%%%  

neu.N_rp=round(neu.t_rp/simu.dt);        % Refractory (in simu.dt)
neu.N_del=round(net.D/simu.dt);          % delay (in simu.dt)

net.V=neu.V_r*ones(net.N,1);            % Voltage vector
net.RI=zeros(net.N,neu.N_del);
net.LS=-(neu.N_del+1)*rand(net.N,1);
net.allspikes=zeros(1,simu.nIterTot);
%ExInput=syn.J*poissrnd(nu_ext*C_ext*simu.dt,net.N,simu.nIterTot);

syn.tau_pre = syn.tau_x;
syn.tau_post = syn.tau_x;

plt.Rasterplot=zeros(net.N,simu.nIterTot);
plt.RasterplotIn=zeros(net.NIn,simu.nIterTot);

% %%%%%%% INITIALIZING GIFS %%%%%%%
LWisimu.dths = (5.*(G.Edges.Weight>0) + 0.8.*(G.Edges.Weight<0)).*abs(G.Edges.Weight);
EColors = [0 1 0].*(G.Edges.Weight>0) + [1 0 0].*(G.Edges.Weight<0);
if gif.graph
    figure(1)
    plot(G,'LineWidth',LWisimu.dths,'EdgeColor',EColors)
    axis tight
    set(gca,'nextplot','replacechildren','visible','off')
    f1 = getframe;
    [im1,map1] = rgb2ind(f1.cdata,256,'nodither');
end

if gif.lapl
    A = net.W + net.W';
    D = diag(sum(A,1));
    L_rw = eye(net.N) - D^(-1)*A;
    [eVals, ~] = eig(L_rw);
    
    figure(2)
    spc = sort(diag(eVals));
    plot(1:20, spc(1:20), '+r')
    xlabel('Index')
    ylabel('Eigenvalue')
    title('Eigendecomposition at time 0s')
    set(gca,'nextplot','replacechildren')
    f2 = getframe;
    [im2,map2] = rgb2ind(f2.cdata,256,'nodither');
end

%% Simulation
tic();

% This part will execute the simulate() function, that takes as arguments:
% - a network net
% - a phase P defined by a duration and PlastON & PlastInON matrices
% The function modifies net.net.W and net.Win, feeds a rasterplot, etc

for phID = 1:nPhase
    [net, plt] = runPhase(net, neu, syn, prot, simu.phases{phID}, plt, splSyn, gif, simu.nIterTot);
end

toc()
progressbar(1)

%% Post-simulation analysis
G = digraph(net.W');
LWisimu.dths = (5.*(G.Edges.Weight>0) + 0.8.*(G.Edges.Weight<0)).*abs(G.Edges.Weight);
EColors = [0 1 0].*(G.Edges.Weight>0) + [1 0 0].*(G.Edges.Weight<0);

A = abs(net.W + net.W');
D = diag(sum(A,1));
L_rw = eye(net.N) - D^(-1)*A;
[eVals, ~] = eig(L_rw);

% filter = repmat(gausswin(100,2.5e-2),net.N,1);
filter = gausswin(100,3);
filter = (1/sum(filter)).*filter;
%for n=1:net.N
%    rateEst(n,:) = (1/simu.dt).*conv(plt.Rasterplot(n,:), filter, 'same');
%end

fMax = floor(500/plt.nbins)*plt.nbins;
edgesRates = 0:fMax/plt.nbins:fMax;

if gif.ratesHist
    rates.min = 145;
    rates.max = 170;
    rates.nbins = 10;
    rates.edges = linspace(rates.min,rates.max,rates.nbins+1);
    rates.synapses = zeros((net.N).^2,2);
    
    figure(4)
    rates.synapses(:,1) = repelem(rateEst(:,1),net.N,1);
    rates.synapses(:,2) = repmat(rateEst(:,1),net.N,1);
    noNulSyns = rates.synapses(abs(net.W)>1e-7,:);
    hist3(noNulSyns,'Edges',{rates.edges, rates.edges},'CDataMode','auto','FaceColor','interp')
    view(0, 90);
    f4 = getframe;
    [im4,map4] = rgb2ind(f4.cdata,256,'nodither');
end

for i=1:simu.nIterTot
    % Rates estimation
    histRates(:,i) = histcounts(rateEst(:,i),edgesRates);  
    if mod(i,gif.updIter)==0 && gif.ratesHist
        rates.synapses(:,1) = repelem(rateEst(:,i),net.N,1);
        rates.synapses(:,2) = repmat(rateEst(:,i),net.N,1);
        noNulSyns = rates.synapses(abs(net.W)>1e-7,:);
        figure(4)
        hist3(noNulSyns,'Edges',{rates.edges, rates.edges},'CDataMode','auto','FaceColor','interp')
        view(0, 90);
        f4 = getframe;
        im4(:,:,1,floor(i/gif.updIter)) = rgb2ind(f4.cdata,map4,'nodither');
    end
    
    % Fourier transform
    if i>floor(0.5*plt.maxPatchSize) && i<=simu.nIterTot - floor(0.5*plt.maxPatchSize)
        % Find the last two spiking events per neuron. Use this to find a
        % second estimate of ISI and determine a patch window
        
        meanISI = 0;
        varISI = 0;
        missM = 0; 
        missS = 0;
        for nerID=1:net.N
            raster = plt.Rasterplot(nerID,1:i);
            lastSpikes = find(raster,11,'last');
            ISI = lastSpikes - circshift(lastSpikes,1);
            ISI = ISI(2:end);

            
            if length(ISI)>1
                varISI = varISI + (std(ISI)^2)./net.N; 
            else
                missS = missS+1;
            end
            
            if length(ISI)>=1
                meanISI = meanISI + mean(ISI)./net.N;
            else
                missM = missM + 1;               
            end
            
        end
        
        plt.ISI(i,1) = meanISI*(net.N-missM)/net.N; 
        plt.ISI(i,2) = sqrt(varISI)*(net.N-missS)/net.N;
        
        patch = plt.Rasterplot(:,i-floor(0.5*plt.maxPatchSize):i+floor(0.5*plt.maxPatchSize));
        
        y = fft(patch,plt.maxPatchSize,2);
        ys = fftshift(y(:,2:end), 2);
        modY = abs(ys);
        phsY = angle(ys)./pi;
        %plt.fourier(:,i-floor(0.5*plt.maxPatchSize)) = fft(patch,plt.maxPatchSize,2);
        
        sigPhase = 0.5*std(phsY); % Measures synchrony
        meanPhase = mean(phsY);
        plt.phaseMetr(i,1) = mean(meanPhase);
        plt.phaseMetr(i,2) = mean(sigPhase);
        sigFreq = std(modY,0,2);
        plt.freqMetr(i,1) = mean(sigFreq);
        %[~,plt.mainModeFreq(i-floor(0.5*plt.maxPatchSize),1)] = max(modY,[],2); % Measures regularity
    end
end

for i=floor(0.5*plt.maxPatchSize)+1:simu.nIterTot - floor(0.5*plt.maxPatchSize)
    
	plt.synMetr(i,1) = mean(plt.phaseMetr(max(i-floor(0.2*plt.maxPatchSize),floor(0.5*plt.maxPatchSize)):min(i+floor(0.2*plt.maxPatchSize),simu.nIterTot-floor(0.5*plt.maxPatchSize)),2));
    plt.regMetr(i,1) = mean(plt.ISI(max(i-floor(0.2*plt.maxPatchSize),floor(0.5*plt.maxPatchSize)):min(i+floor(0.2*plt.maxPatchSize),simu.nIterTot-floor(0.5*plt.maxPatchSize)),2)./plt.ISI(max(i-floor(0.2*plt.maxPatchSize),floor(0.5*plt.maxPatchSize)):min(i+floor(0.2*plt.maxPatchSize),simu.nIterTot-floor(0.5*plt.maxPatchSize)),1));
        % Determining regime
        
    plt.regime(i,1) = (plt.synMetr(i,1) < plt.regimeThr(1,1)); % SYNCHRONICITY
    plt.regime(i,2) = (plt.regMetr(i,1) < plt.regimeThr(1,2)); % REGULARITY
    
end

% Cutting regions per regime
SR = and(plt.regime(:,1)==1, plt.regime(:,2)==1);
AR = and(plt.regime(:,1)==0, plt.regime(:,2)==1);
SI = and(plt.regime(:,1)==1, plt.regime(:,2)==0);
AI = and(plt.regime(:,1)==0, plt.regime(:,2)==0);

%% Creating graph for visualization
if gif.graph
    imwrite(im1, map1, strcat(env.outputRoot, 'Figures/Network/sampleGraph_randinit.gif'), 'DelayTime',0, 'LoopCount',inf)
end

if gif.lapl
    imwrite(im2,map2, strcat(env.outputRoot, 'Figures/Network/specClust_asym.gif'), 'DelayTime',0, 'LoopCount',inf)
end

if gif.ratesHist
    imwrite(im4,map4, strcat(env.outputRoot, 'Figures/Network/ratesHist.gif'), 'DelayTime',0, 'LoopCount',inf)
end

%% Plotting network stats

set(0,'DefaultFigureWindowStyle','docked')

RasterIn = plt.RasterplotIn(1:ceil(0.1*net.NIn),:);
RasterNet = plt.Rasterplot(ceil(0.1*net.N):net.N,:);
Rasterplot = [RasterIn; zeros(15,simu.nIterTot); RasterNet];

switch plt.all.raster
    case 1
        fig.raster = figure('Name','NET_Rasterplot_FULL','NumberTitle','off');
        [I1,I2] = find(Rasterplot);
        totActSnaps = sum(Rasterplot);
        ax1 = subplot(2,1,1);
        imagesc(Rasterplot)
        xticklabels(simu.dt.*xticks)
        ax2 = subplot(2,1,2);
        
        SRspikes = totActSnaps;
        SRspikes(1,1:max(1,floor(0.5*plt.maxPatchSize))) = NaN;
        SRspikes(1,max(end-floor(0.5*plt.maxPatchSize),1):end) = NaN;
        SRspikes(SI) = NaN;
        SRspikes(AR) = NaN;
        SRspikes(AI) = NaN;
        
        ARspikes = totActSnaps; 
        ARspikes(1,1:max(1,floor(0.5*plt.maxPatchSize))) = NaN;
        ARspikes(1,max(end-floor(0.5*plt.maxPatchSize),1):end) = NaN;        
        ARspikes(SI) = NaN;
        ARspikes(SR) = NaN;
        ARspikes(AI) = NaN;
        
        SIspikes = totActSnaps; 
        SIspikes(1,1:max(1,floor(0.5*plt.maxPatchSize))) = NaN;
        SIspikes(1,max(end-floor(0.5*plt.maxPatchSize),1):end) = NaN;  
        SIspikes(SR) = NaN;
        SIspikes(AR) = NaN;
        SIspikes(AI) = NaN;
        
        AIspikes = totActSnaps;    
        AIspikes(1,1:max(1,floor(0.5*plt.maxPatchSize))) = NaN;
        AIspikes(1,max(end-floor(0.5*plt.maxPatchSize),1):end) = NaN;         
        AIspikes(SI) = NaN;
        AIspikes(AR) = NaN;
        AIspikes(SR) = NaN;
        
        p = plot((1:simu.nIterTot), SRspikes, (1:simu.nIterTot), ARspikes, (1:simu.nIterTot), SIspikes, (1:simu.nIterTot), AIspikes); 
        p(1).LineWidth = 1;
        p(1).Color = 'r';
        p(2).LineWidth = 1;
        p(2).Color = 'g';
        p(3).LineWidth = 1;
        p(3).Color = 'b';
        p(4).LineWidth = 1;
        p(4).Color = 'k';        
        
        xticklabels(simu.dt.*xticks)
        title('Spiking activity in neural population')
        ax1.Position(1,2) = ax1.Position(1,2) - 0.7*ax2.Position(1,4);
        ax1.Position(1,4) = 1.7*ax1.Position(1,4);
        ax2.Position(1,4) = 0.3*ax2.Position(1,4);
 
        text(0.04,-0.45, strcat('N: ', num2str(net.N)), 'FontSize', 15, 'Units', 'normalized');
        text(0.12,-0.45, strcat('C_E: ', num2str(net.Connectivity*net.NE)), 'FontSize', 15, 'Units', 'normalized');
        text(0.04,-0.75, strcat('g: ', num2str(net.g)), 'FontSize', 15, 'Units', 'normalized');
        text(0.12,-0.75, strcat('nu_{rel}: ', num2str(net.rExtRel)), 'FontSize', 15, 'Units', 'normalized'); 
        text(0.2,-0.45, strcat('N_{In}: ', num2str(net.NIn)), 'FontSize', 15, 'Units', 'normalized');
        text(0.2,-0.75, strcat('D: ', num2str(net.D)), 'FontSize', 15, 'Units', 'normalized');
        text(0.28,-0.45, strcat('C_{pre}: ', num2str(syn.C_pre)), 'FontSize', 15, 'Units', 'normalized');
        text(0.28,-0.75, strcat('C_{post}: ', num2str(syn.C_post)), 'FontSize', 15, 'Units', 'normalized');
        text(0.36,-0.45, strcat('\gamma_p: ', num2str(syn.gamma_pot)), 'FontSize', 15, 'Units', 'normalized');
        text(0.36,-0.75, strcat('\gamma_d: ', num2str(syn.gamma_dep)), 'FontSize', 15, 'Units', 'normalized');
        text(0.44,-0.45, strcat('S_{attr}: ', num2str(syn.S_attr)), 'FontSize', 15, 'Units', 'normalized');         
        
        %case 1
        rasterSnaps = zeros(ceil(0.1*net.NIn) + net.N + 15 + 1 -ceil(0.1*net.N), plt.timeSpl.n*ceil(plt.timeSpl.dur/simu.dt) + (plt.timeSpl.n-1)*plt.timeSpl.inter);
        totActSnaps = zeros(1, plt.timeSpl.n*ceil(plt.timeSpl.dur/simu.dt) + (plt.timeSpl.n-1)*plt.timeSpl.inter);
        xticksList=[]; xtickvalsList=[];
        for tSplID = 1:plt.timeSpl.n
            rasterSnaps(:, 1+(tSplID-1)*(ceil(plt.timeSpl.dur/simu.dt)+plt.timeSpl.inter):tSplID*ceil(plt.timeSpl.dur/simu.dt)+(tSplID-1)*plt.timeSpl.inter) ...
                = Rasterplot(:,1+floor((tSplID-1)/(plt.timeSpl.n-1)*(simu.nIterTot-ceil(plt.timeSpl.dur/simu.dt)-1)):floor((tSplID-1)/(plt.timeSpl.n-1)*(simu.nIterTot-ceil(plt.timeSpl.dur/simu.dt)-1)+ceil(plt.timeSpl.dur/simu.dt)));           

            rasterSnaps(:, 1 + tSplID*ceil(plt.timeSpl.dur/simu.dt)+(tSplID-1)*plt.timeSpl.inter : plt.timeSpl.inter + tSplID*ceil(plt.timeSpl.dur/simu.dt)+(tSplID-1)*plt.timeSpl.inter) ... 
                = 1;
            
            xticksList = cat(2, xticksList, ...
                1+(tSplID-1)*(ceil(plt.timeSpl.dur/simu.dt)+plt.timeSpl.inter) + ceil((1/3)*(plt.timeSpl.dur/simu.dt).*(0:1:3))...
                );
            
            xtickvalsList = cat(2, xtickvalsList, ...
                simu.dt.*(1+floor((tSplID-1)/(plt.timeSpl.n-1)*(simu.nIterTot-ceil(plt.timeSpl.dur/simu.dt)-1)) + ceil((1/3)*(plt.timeSpl.dur/simu.dt).*(0:1:3)))...
                );
        end    
        totActSnaps = sum(rasterSnaps);
        for tSplID = 1:plt.timeSpl.n
            totActSnaps(1, 1 + tSplID*ceil(plt.timeSpl.dur/simu.dt)+(tSplID-1)*plt.timeSpl.inter : plt.timeSpl.inter + tSplID*ceil(plt.timeSpl.dur/simu.dt)+(tSplID-1)*plt.timeSpl.inter) = 0;
        end
        fig.rasterSpl = figure('Name','NET_Rasterplot_SAMPLES','NumberTitle','off');
        ax1 = subplot(2,1,1);
        imagesc(rasterSnaps)
        xticks(xticksList)
        xticklabels(xtickvalsList)
        ax2 = subplot(2,1,2);
        bar(totActSnaps)
        xticks(xticksList)
        xticklabels(xtickvalsList)
        ax1.Position(1,2) = ax1.Position(1,2) - 0.7*ax2.Position(1,4);
        ax1.Position(1,4) = 1.7*ax1.Position(1,4);
        ax2.Position(1,4) = 0.3*ax2.Position(1,4);
        
        
        text(0.04,-0.45, strcat('N: ', num2str(net.N)), 'FontSize', 15, 'Units', 'normalized');
        text(0.12,-0.45, strcat('C_E: ', num2str(net.Connectivity*net.NE)), 'FontSize', 15, 'Units', 'normalized');
        text(0.04,-0.75, strcat('g: ', num2str(net.g)), 'FontSize', 15, 'Units', 'normalized');
        text(0.12,-0.75, strcat('nu_{rel}: ', num2str(net.rExtRel)), 'FontSize', 15, 'Units', 'normalized'); 
        text(0.2,-0.45, strcat('N_{In}: ', num2str(net.NIn)), 'FontSize', 15, 'Units', 'normalized');
        text(0.2,-0.75, strcat('D: ', num2str(net.D)), 'FontSize', 15, 'Units', 'normalized');
        text(0.28,-0.45, strcat('C_{pre}: ', num2str(syn.C_pre)), 'FontSize', 15, 'Units', 'normalized');
        text(0.28,-0.75, strcat('C_{post}: ', num2str(syn.C_post)), 'FontSize', 15, 'Units', 'normalized');
        text(0.36,-0.45, strcat('\gamma_p: ', num2str(syn.gamma_pot)), 'FontSize', 15, 'Units', 'normalized');
        text(0.36,-0.75, strcat('\gamma_d: ', num2str(syn.gamma_dep)), 'FontSize', 15, 'Units', 'normalized');
        text(0.44,-0.45, strcat('S_{attr}: ', num2str(syn.S_attr)), 'FontSize', 15, 'Units', 'normalized');           
        
        getRegimeMetrics(simu, plt, fig);
        
end

if plt.spl.ca  % DEPRECATED
    fig.ca = figure('Name','NET_Calcium','NumberTitle','off');
    imagesc(splSyn.ca)
    title('Synaptic calcium activity in sample synapses')
    colorbar
end

if plt.spl.rho  % DEPRECATED
    fig.rho = figure('Name','NET_Phospho','NumberTitle','off');
    imagesc(splSyn.rho)
    title('Phosphorylation state at sample synapses')
    colorbar
end

if plt.spl.w  % DEPRECATED
    splSyn.wHist = synSign(splSyn.IDs).*transfer(splSyn.rho, prot);
    fig.w = figure('Name','NET_Weights','NumberTitle','off');
    imagesc(splSyn.wHist)
    title('Synaptic weight at sample synapses')
    colorbar
end

splSyn.stats = cat(2, linspace(1,splSyn.n,splSyn.n)', splSyn.PreNeurons, splSyn.PostNeurons, splSyn.w(:,1), splSyn.w(:,end));

switch plt.spl.pres
    case 1
        splSyn.pres = zeros(3*splSyn.n, simu.nIterTot);
        for i=1:splSyn.n
            splSyn.pres(3*(i-1)+1,:) = 100*plt.Rasterplot(splSyn.PreNeurons(i,1),:);
            splSyn.pres(3*(i-1)+2,:) = 100*plt.Rasterplot(splSyn.PostNeurons(i,1),:);
            splSyn.pres(3*(i-1)+3,:)= (1/syn.J).*synSign(splSyn.IDs(i)).*splSyn.rho(i,:);
        end
        fig.pres = figure('Name','NET_Summary','NumberTitle','off');
        imagesc(splSyn.pres)
        colorbar
    case 2
        splSyn.pres = zeros(5*floor(splSyn.n/2), simu.nIterTot);
        for i=1:splSyn.n
            splSyn.pres(5*(i-1)+1,:) = plt.Rasterplot(splSyn.PreNeurons(i,1),:);
            splSyn.pres(5*(i-1)+2,:) = plt.Rasterplot(splSyn.PostNeurons(i,1),:);
            splSyn.pres(5*(i-1)+3,:)= splSyn.ca(i,:);
            splSyn.pres(5*(i-1)+4,:)= (1/(syn.J*syn.rho_max)).*net.synSign(splSyn.IDs(i)).*splSyn.rho(i,:);
            splSyn.pres(5*(i-1)+5,:)= (1/syn.J).*splSyn.w(i,:);
        end
        fig.pres = figure('Name','NET_Summary','NumberTitle','off');
        imagesc(splSyn.pres)
        xticklabels(simu.dt.*xticks)
        colorbar
    case 3
        splSyn.pres = zeros(5*floor(splSyn.n/2), plt.timeSpl.n*ceil(plt.timeSpl.dur/simu.dt) + (plt.timeSpl.n-1)*plt.timeSpl.inter);
        xticksList=[]; xtickvalsList=[];
        for tSplID = 1:plt.timeSpl.n
            for i=1:splSyn.n
                splSyn.pres(5*(i-1)+1, 1+(tSplID-1)*(ceil(plt.timeSpl.dur/simu.dt)+plt.timeSpl.inter):tSplID*ceil(plt.timeSpl.dur/simu.dt)+(tSplID-1)*plt.timeSpl.inter) ...
                    = plt.Rasterplot(splSyn.PreNeurons(i,1),1+floor((tSplID-1)/(plt.timeSpl.n-1)*(simu.nIterTot-ceil(plt.timeSpl.dur/simu.dt)-1)):floor((tSplID-1)/(plt.timeSpl.n-1)*(simu.nIterTot-ceil(plt.timeSpl.dur/simu.dt)-1)+ceil(plt.timeSpl.dur/simu.dt)));           

                splSyn.pres(5*(i-1)+2, 1+(tSplID-1)*(ceil(plt.timeSpl.dur/simu.dt)+plt.timeSpl.inter):tSplID*ceil(plt.timeSpl.dur/simu.dt)+(tSplID-1)*plt.timeSpl.inter) ... 
                    = plt.Rasterplot(splSyn.PostNeurons(i,1),1+floor((tSplID-1)/(plt.timeSpl.n-1)*(simu.nIterTot-ceil(plt.timeSpl.dur/simu.dt)-1)):floor((tSplID-1)/(plt.timeSpl.n-1)*(simu.nIterTot-ceil(plt.timeSpl.dur/simu.dt)-1)+ceil(plt.timeSpl.dur/simu.dt)));

                splSyn.pres(5*(i-1)+3, 1+(tSplID-1)*(ceil(plt.timeSpl.dur/simu.dt)+plt.timeSpl.inter):tSplID*ceil(plt.timeSpl.dur/simu.dt)+(tSplID-1)*plt.timeSpl.inter) ...
                    = splSyn.ca(i,1+floor((tSplID-1)/(plt.timeSpl.n-1)*(simu.nIterTot-ceil(plt.timeSpl.dur/simu.dt)-1)):floor((tSplID-1)/(plt.timeSpl.n-1)*(simu.nIterTot-ceil(plt.timeSpl.dur/simu.dt)-1)+ceil(plt.timeSpl.dur/simu.dt)));

                splSyn.pres(5*(i-1)+4, 1+(tSplID-1)*(ceil(plt.timeSpl.dur/simu.dt)+plt.timeSpl.inter):tSplID*ceil(plt.timeSpl.dur/simu.dt)+(tSplID-1)*plt.timeSpl.inter) ...
                    = (1/(syn.J*syn.rho_max)).*net.synSign(splSyn.IDs(i)).*splSyn.rho(i,1+floor((tSplID-1)/(plt.timeSpl.n-1)*(simu.nIterTot-ceil(plt.timeSpl.dur/simu.dt)-1)):floor((tSplID-1)/(plt.timeSpl.n-1)*(simu.nIterTot-ceil(plt.timeSpl.dur/simu.dt)-1)+ceil(plt.timeSpl.dur/simu.dt)));

                splSyn.pres(5*(i-1)+5, 1+(tSplID-1)*(ceil(plt.timeSpl.dur/simu.dt)+plt.timeSpl.inter):tSplID*ceil(plt.timeSpl.dur/simu.dt)+(tSplID-1)*plt.timeSpl.inter) ...
                    = (1/syn.J).*splSyn.w(i,1+floor((tSplID-1)/(plt.timeSpl.n-1)*(simu.nIterTot-ceil(plt.timeSpl.dur/simu.dt)-1)):floor((tSplID-1)/(plt.timeSpl.n-1)*(simu.nIterTot-ceil(plt.timeSpl.dur/simu.dt)-1)+ceil(plt.timeSpl.dur/simu.dt)));
                
            end
            
            xticksList = cat(2, xticksList, ...
                1+(tSplID-1)*(ceil(plt.timeSpl.dur/simu.dt)+plt.timeSpl.inter) + ceil((1/3)*(plt.timeSpl.dur/simu.dt).*(0:1:3))...
                );
            
            xtickvalsList = cat(2, xtickvalsList, ...
                simu.dt.*(1+floor((tSplID-1)/(plt.timeSpl.n-1)*(simu.nIterTot-ceil(plt.timeSpl.dur/simu.dt)-1)) + ceil((1/3)*(plt.timeSpl.dur/simu.dt).*(0:1:3)))...
                );
        end
        fig.pres = figure('Name','NET_Summary','NumberTitle','off');
        imagesc(splSyn.pres)
        xticks(xticksList)
        xticklabels(xtickvalsList)
        colorbar
end

if plt.spl.hist
    anyPlastInput = 0;
    anyPlastExc = 0;
    anyPlastInh = 0;
    for phID = 1:nPhase
        anyPlastInput = anyPlastInput + simu.phases{phID}.InE + simu.phases{phID}.InI;
        anyPlastExc = anyPlastExc + simu.phases{phID}.EE + simu.phases{phID}.EI;
        anyPlastInh = anyPlastInh + simu.phases{phID}.IE + simu.phases{phID}.II;
    end
    
    tickList = ceil(plt.nbins.*(0.1:0.1:1));
    if anyPlastExc
        fig.excHist = figure('Name','NET_WghHistExcit','NumberTitle','off');
        ax1 = subplot(2,1,1);
        imagesc(log(plt.histW_exc));
        set(gca,'YDir','normal')
        title('Evolution of weight distribution for excitatory synapses')
        xlabel('Time')
        xticklabels(simu.dt.*xticks)
        ylabel('Synaptic weight')
        yticks(tickList)
        yticklabels(plt.edgesW_exc(1, tickList))
        colorbar
        Wexc_hStep = plt.edgesW_exc(2) - plt.edgesW_exc(1);
        valsWexc =plt. edgesW_exc + Wexc_hStep;
        valsWexc = repmat(valsWexc(:,1:end-1)',1,simu.nIterTot);
        sumWexc = sum(valsWexc.*plt.histW_exc);
        ax2 = subplot(2,1,2);
        bar(plt.edgesW_exc(2:end), plt.histW_exc(:,end))
        ax1.Position(1,2) = ax1.Position(1,2) - 0.5*ax2.Position(1,4);
        ax1.Position(1,4) = 1.5*ax1.Position(1,4);
        ax2.Position(1,4) = 0.5*ax2.Position(1,4);
    end

    if anyPlastInh
        fig.inhHist = figure('Name','NET_WghHistInhib','NumberTitle','off');
        ax1 = subplot(2,1,1);
        imagesc(log(plt.histW_inh));
        set(gca,'YDir','normal')
        title('Evolution of weight distribution for inhibitory synapses')
        xlabel('Time')
        xticklabels(simu.dt.*xticks)
        ylabel('Synaptic weight')
        yticks(tickList)
        yticklabels(plt.edgesW_inh(1, tickList))
        colorbar
        Winh_hStep = plt.edgesW_inh(2) - plt.edgesW_inh(1);
        valsWinh = plt.edgesW_inh + Winh_hStep;
        valsWinh = repmat(valsWinh(:,1:end-1)',1,simu.nIterTot);
        sumWinh = sum(valsWinh.*plt.histW_inh);
        ax2 = subplot(2,1,2);
        bar(plt.edgesW_inh(2:end), plt.histW_inh(:,end))
        ax1.Position(1,2) = ax1.Position(1,2) - 0.5*ax2.Position(1,4);
        ax1.Position(1,4) = 1.5*ax1.Position(1,4);
        ax2.Position(1,4) = 0.5*ax2.Position(1,4);
    end
    
    fig.rhoHist = figure('Name','NET_PhosphoHist','NumberTitle','off');
    ax1 = subplot(2,1,1);
    imagesc(log(plt.histRho));
    set(gca,'YDir','normal')
    title('Evolution of phosphorylation distribution (all synapses)')
    xlabel('Time')
    xticklabels(simu.dt.*xticks)
    ylabel('Phosphorylation level')
    yticks(tickList)
    yticklabels(plt.edgesRho(1, tickList))
    colorbar
    ax2 = subplot(2,1,2);
    bar(plt.edgesRho(2:end), plt.histRho(:,end))
    ax1.Position(1,2) = ax1.Position(1,2) - 0.5*ax2.Position(1,4);
    ax1.Position(1,4) = 1.5*ax1.Position(1,4);
    ax2.Position(1,4) = 0.5*ax2.Position(1,4);
    
    
    % sprintf('Excitatory: %0.3f - Inhibitory: %0.3f - Total: %0.3f', sumWexc(1,end)/syn.J, -sumWinh(1,end)/syn.J, (sumWexc(1,end) + sumWinh(1,end))/syn.J)
    fig.rateHist = figure('Name','NET_RatesHist','NumberTitle','off');
    ax1 = subplot(2,1,1);
    imagesc(log(histRates));
    set(gca,'YDir','normal')
    for n=1:splNeu.NE
        hold on
        plot((plt.nbins/edgesRates(1,end)).*rateEst(splNeu.ExcNeurons(n,1),:), 'g', 'LineWidth', 1)
    end
    for n=1:splNeu.NI
        hold on
        plot((plt.nbins/edgesRates(1,end)).*rateEst(splNeu.InhNeurons(n,1),:), 'r', 'LineWidth', 1)
    end
    title('Evolution of firing rates distribution')
    xlabel('Time')
    xticklabels(simu.dt.*xticks)
    ylabel('Firing rate')
    yticks(tickList)
    yticklabels(edgesRates(1, tickList))
    colorbar
    
    ax2 = subplot(2,1,2);
    bar(edgesRates(2:end), mean(histRates(:,end-50:end),2))
    ax1.Position(1,2) = ax1.Position(1,2) - 0.5*ax2.Position(1,4);
    ax1.Position(1,4) = 1.5*ax1.Position(1,4);
    ax2.Position(1,4) = 0.5*ax2.Position(1,4);
    
end

switch plt.spl.phase
    case 2
        
        fig.phase = figure('Name','NET_BrunelPhase','NumberTitle','off');
        c = rand(1,3);
        plot((1/net.meanWexc(2,1)).*net.meanWexc(2:end,1), (1/net.meanWinh(2,1)).*net.meanWinh(2:end,1), '-', 'Color', c)
        for i=1:20:min(simu.nIterTot,201)
            text((1/net.meanWexc(2,1)).*net.meanWexc(1+i,1), (1/net.meanWinh(2,1)).*net.meanWinh(i+1,1),num2str(1+i),'Color', c)
        end
        for i=301:100:min(simu.nIterTot,1001)
            text((1/net.meanWexc(2,1)).*net.meanWexc(1+i,1), (1/net.meanWinh(2,1)).*net.meanWinh(i+1,1),num2str(1+i),'Color', c)
        end
        for i=2001:1000:min(simu.nIterTot,30001)
            text((1/net.meanWexc(2,1)).*net.meanWexc(1+i,1), (1/net.meanWinh(2,1)).*net.meanWinh(i+1,1),num2str(1+i),'Color', c)
        end
        title('Evolution of average synaptic weights')
        xlabel('Relative mean excitatory weight')
        ylabel('Relative mean inhibitory weight')
        phaseParamPos.x = 10; phaseParamPos.y = 10;
        phaseParamPos.colSep = 5; phaseParamPos.rowSep = 0.5; phaseParamPos.ftSize = 3;
        stampNetworkParams( fig.phase, ...
            ['g','\nu_{rel}','\epsilon','C_{pre}','C_{post}','\theta_p','\theta_d', 'S_attr'], ...
            [net.g, net.rExtRel, net.Connectivity, syn.C_pre, syn.C_post, syn.theta_pot, syn.theta_dep, syn.S_attr], ...
            phaseParamPos);
        
    case 1
        fig.phase = figure('Name','NET_BrunelPhase','NumberTitle','off');
        c = rand(1,3);
        plot((1/net.meanWexc(2,1)).*net.meanWexc(2:end,1), -net.meanWinh(2:end,1)./net.meanWexc(2:end,1), '-', 'Color', c)
        for i=1:20:min(simu.nIterTot,201)
            text((1/net.meanWexc(2,1)).*net.meanWexc(1+i,1), -net.meanWinh(1+i,1)./net.meanWexc(1+i,1),num2str(1+i),'Color', c)
        end
        for i=301:100:min(simu.nIterTot,1001)
            text((1/net.meanWexc(2,1)).*net.meanWexc(1+i,1), -net.meanWinh(1+i,1)./net.meanWexc(1+i,1),num2str(1+i),'Color', c)
        end
        for i=2001:1000:min(simu.nIterTot,20001)
            text((1/net.meanWexc(2,1)).*net.meanWexc(1+i,1), -net.meanWinh(1+i,1)./net.meanWexc(1+i,1),num2str(1+i),'Color', c)
        end
        title('Evolution of average synaptic weights')
        xlabel('Relative mean excitatory weight')
        ylabel('Relative mean inhibitory weight')
        phaseParamPos.x = 10; phaseParamPos.y = 10;
        phaseParamPos.colSep = 50; phaseParamPos.rowSep = 8; phaseParamPos.ftSize = 8;
        % ylim([min(3, min(-net.meanWinh(2:end,1)./net.meanWexc(2:end,1))) max(8, max(-net.meanWinh(2:end,1)./net.meanWexc(2:end,1)))])

%         SRtoAI_thr = refline([0 4.2]);
%         SRtoAI_thr.Color = 'r';
% 
%         AItoSI_thr = refline([0 5.3]);
%         AItoSI_thr.Color = 'g';
        
        stampNetworkParams( fig.phase, ...
            ["g","\nu_{rel}","\epsilon","C_{pre}","C_{post}","\theta_p","\theta_d","S_attr"], ...
            [net.g, net.rExtRel, net.Connectivity, syn.C_pre, syn.C_post, syn.theta_pot, syn.theta_dep, syn.S_attr], ...
            phaseParamPos);
end

%% Data logging

% Output the parameters as well as some measures in a CSV
outputFolder = '../Outputs/Numerical/Network/';
outputFile = strcat(outputFolder,'out.csv');

% Parameter structure
out.ID = getSimuID(outputFile);
out.date = date;
out = concatStruct(out,syn);
out = concatStruct(out,simu);
out = concatStruct(out,net);
out = concatStruct(out,neu);
out = concatStruct(out,init);

out.finalMeanWExc = net.meanWexc(end,1);
out.finalMeanWinh = net.meanWinh(end,1);
out.finalG = -net.meanWinh(end,1)./net.meanWexc(end,1);
%out.finalNu = -net.meanWinh(end,1)./net.meanWexc(end,1);
out.finalFiringRate = mean(rateEst(:,end));

if 1 % Removing some fields that we do not want in the output
    out = rmfield(out, 'phases');
    out = rmfield(out, 'W');
    out = rmfield(out, 'WIn');
    out = rmfield(out, 'synSign');
    out = rmfield(out, 'meanWexc');
    out = rmfield(out, 'meanWinh');
    out = rmfield(out, 'ca');
    out = rmfield(out, 'xpre');
    out = rmfield(out, 'xpost');
    out = rmfield(out, 'rho');
    out = rmfield(out, 'actPot');
    out = rmfield(out, 'actDep');
    out = rmfield(out, 'caIn');
    out = rmfield(out, 'xpreIn');
    out = rmfield(out, 'xpostIn');
    out = rmfield(out, 'rhoIn');
    out = rmfield(out, 'actPotIn');
    out = rmfield(out, 'actDepIn');
    out = rmfield(out, 'V');
    out = rmfield(out, 'RI');
    out = rmfield(out, 'LS');
    out = rmfield(out, 'allspikes');
end

% Writing to CSV
struct2csv(out,outputFile)

%% Output macroscopic variables
out.meanWexc = net.meanWexc(end,1);
out.meanWinh = net.meanWinh(end,1);
sqWexc = (net.W - net.meanWexc(end,1)).^2;
sqWinh = (net.W - net.meanWinh(end,1)).^2;
out.stdWexc = sqrt(mean(sqWexc(net.W>0)));
out.stdWinh = sqrt(mean(sqWinh(net.W<0)));

out.synScore = plt.synMetr(end,1);
out.regScore = plt.regMetr(end,1);

out.synRegime = plt.regime(end,1);
out.regRegime = plt.regime(end,2);