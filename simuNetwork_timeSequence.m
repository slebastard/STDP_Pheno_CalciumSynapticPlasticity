% Simulation of the Brunel JCNS 2000 network with electrode recording
% spikes (for paper Destexhe Touboul on Criticality).

% Modified by S. Lebastard, 2018
% Should accomodate mutliple models of plasticity, including none
% Should allow for heterogeneous synapse rules, ie manipulation of matrix
% of parameters

% ToDo
% - Refactorize code to gain modularity, be more readable
% - Allow user to start either with uniform weights or with random weights
% - Create clusterization index and look at its evolution
%   For this purpose, use techniques from spectral clustering

clear all
close all

% Add plot tools to path
env = getEnv();
addpath(genpath(env.functionsRoot));


%% Parameterization

% %%%%%%%%   PARAMETERS OF THE SIMULATION  %%%%%%%% 

syn = getSynapse();

prot.TD = 0;
prot.S_attr = syn.sAttr;
prot.noise_lvl = syn.sigma;

simu.T = 15;
simu.dt=5e-4;
Iterations=ceil(simu.T/simu.dt);


plt.all.raster = 0;
plt.spl.ca = 0;
plt.spl.rho = 0;
plt.spl.w = 0;
plt.spl.pres = 0;
plt.spl.hist = 2;
plt.spl.phase = 0;

plt.timeSpl.n = 5;
plt.timeSpl.dur = 0.02*simu.T;
plt.timeSpl.inter = 70;


gif.graph = 0;
    gif.minN = 15;
    gif.maxN = 40;
gif.lapl = 0;
gif.W = 1;

init.w = 'norm';
init.c = 0.7;
init.mode = 'rand';
init.strap = 0.2;     % Time (s) for which the system should run before monitoring starts
init.strapIters = ceil(init.strap/simu.dt);


% %%%%%%%%   PARAMETERS OF THE NETWORK  %%%%%%%%   

net.N=100;
net.NE=ceil(0.8*net.N);      % # excitatory neurons
net.NI=net.N-net.NE;             % # excitatory neurons

net.Connectivity=0.1;           % Connectivity coefficients, parameter epsilon in Brunel 2000.
net.D=3e-3;                        % Transmission Delay
CE=round(net.Connectivity*net.NE);          % Number of exc connections
CI=round(net.Connectivity*net.NI);          % Number of inh connections
C_ext=CE;


% %%%%%%%%   PARAMETERS OF THE NEURONS  %%%%%%%%   

% Threshold, reset, time constant and refractory period (p.185, 2nd column,
% Brunel 2000)
neu.V_t=20e-3;
neu.V_r=10e-3;
neu.tau=20e-3;
neu.t_rp=2e-3;

% Bifurcation parameters
syn.J=(1/sqrt(net.N))*0.2e-3; %0.1e-3;                   % Strength of exc. Connections
net.g=1.5;                      % Strength of inh/exc (inh connections -gJ)
net.rExtRel=3.5;         % nu_ext/nu_thresh

nu_thresh=neu.V_t/(syn.J*CE*neu.tau);   % Frequency needed for a neuron to reach the threshold. 
nu_ext=net.rExtRel*nu_thresh;       % external Poisson input rate


% %%%%%%%%   PARAMETERS OF THE ELECTRODES  %%%%%%%%  
L=100;          % Size of the cortex
NCX=8;          % Number of electrodes rows
NCY=8;          % Number of electrodes columns
NC=NCX*NCY;     % Total number of electrodes

% Parameters of histograms
nbins = 100;
edgesRho = linspace(0,syn.rhoMax,nbins+1);
edgesW_exc = linspace(0,1,nbins+1);
edgesW_inh = linspace(-net.g,0,nbins+1);

%% Creating network

% %%%%%% INITIALIZING SYNAPTIC QTIES %%%%%%

% %%% Synaptic weights %%%
W = zeros(net.N,net.N);   % 1st index for postsynaptic neuron, 2nd for presynaptic

for i=1:net.N
    ECells=randperm(net.NE);
    ECells=ECells(1:CE);
    W(i,ECells) = strcmp(init.mode, 'rand').*2.*init.c.*syn.J.*rand(1,CE) + ~strcmp(init.mode, 'rand').*init.c.*syn.J;
    
    ICells=randperm(net.NI);
    ICells=net.NE+ICells(1:CI);
    W(i,ICells)=-strcmp(init.mode, 'rand').*2.*net.g.*init.c.*syn.J.*rand(1,CI) -~strcmp(init.mode, 'rand').*init.c .*net.g*syn.J;
end

synSign = syn.J.*(W>0) - net.g*syn.J.*(W<0);
[~,excNeurons] = ind2sub(size(W), find(W>0));
[~,inhNeurons] = ind2sub(size(W), find(W<0));

% Mean stats for phase plot
meanWexc = mean(W((1/syn.J).*W>0)).*ones(Iterations+1,1);
meanWinh = mean(W((1/syn.J).*W<0)).*ones(Iterations+1,1);

rateEst = zeros(net.N, Iterations);
histRates = zeros(nbins, Iterations);

subIDsExc = excNeurons( ...
    randperm(size(excNeurons,1),floor(0.8*min(max(gif.minN, 0.3*net.N), gif.maxN))) ...
    );
subIDsInh = inhNeurons( ...
    randperm(size(inhNeurons,1),ceil(0.2*min(max(gif.minN, 0.3*net.N), gif.maxN))) ...
    );
G = digraph(W([subIDsExc;subIDsInh]',[subIDsExc;subIDsInh]')');    % Graph for plotting evolution of network
A = abs(W + W'); % Weight matrix for spectral clustering
% Symmetrization, all information on directionnality is lost
histRho = zeros(nbins, Iterations);
histW_exc = zeros(nbins, Iterations);
histW_inh = zeros(nbins, Iterations);

% %%% Other qties %%%
ca = zeros(net.N,net.N);
xpre = ones(net.N,net.N);
xpost = ones(net.N,net.N);
% rho = zeros(net.N,net.N);
rho = transferinv(W./synSign, syn.sAttr, syn.sigma, syn.rhoMax);
actPot = zeros(net.N,net.N);
actDep = zeros(net.N,net.N);


% %%% Electrodes wiring %%%
% Electrodes are regularly located, we compute the attenuation coefficient
% due to the distance. 

[X,Y] = meshgrid((1:8)',(1:8)');
POSC =10* [X(:) Y(:)];              % Positions of electrodes
POSN=L*rand(net.N,2);                   % Positions of the neurons

DCN=zeros(NCX*NCY,net.N);
for i=1:(NCX*NCY)
    for j=1:net.N
        DCN(i,j)=sqrt((POSC(i,1)-POSN(j,1))^2+(POSC(i,2)-POSN(j,2))^2);
    end
end
DCN=(DCN.^(-2)).*(DCN<10);
%W=rand(net.N,net.N)<=net.Connectivity;


% Sample synapses
splSyn.n = 15;
[synOut, synIn] = find(W);
perm = randperm(length(synOut));
rpOut = synOut(perm);
rpIn = synIn(perm);

splSyn.PostNeurons = rpOut(1:splSyn.n);
splSyn.PreNeurons = rpIn(1:splSyn.n);
splSyn.IDs = sub2ind(size(ca),splSyn.PostNeurons,splSyn.PreNeurons);
splSyn.ca = zeros(splSyn.n,Iterations);
splSyn.xpre = zeros(splSyn.n,Iterations);
splSyn.xpost = zeros(splSyn.n,Iterations);
splSyn.rho = zeros(splSyn.n,Iterations);
splSyn.w = zeros(splSyn.n,Iterations);

% Sample neurons
splNeu.NE = 3;
splNeu.NI = 2;
splNeu.ExcNeurons = excNeurons(randperm(size(excNeurons,1),splNeu.NE));
splNeu.InhNeurons = inhNeurons(randperm(size(inhNeurons,1),splNeu.NI));

%% Processing time signal, initializing filters

% Reading input signal
% in.fileName = 'test.signal';
% in = readSignal(in.fileName); % k (dim) x M (timeSteps) signal

in.dim = 3;
in.timeSteps = 2200;
in.data = zeros(in.dim, in.timeSteps);
in.data(1,1:100) = gausswin(100,3);
in.data(1,1001:1100) = gausswin(100,3);
in.data(1,2001:2100) = gausswin(100,3);
in.data(2,51:150) = gausswin(100,3);
in.data(2,1051:1150) = gausswin(100,3);
in.data(2,2051:2150) = gausswin(100,3);
in.data(3,101:200) = gausswin(100,3);
in.data(3,1101:1200) = gausswin(100,3);
in.data(3,2101:2200) = gausswin(100,3);

if Iterations > in.timeSteps
    in.data = cat(1, in.data, zeros(in.dim, Iterations-in.timeSteps));
elseif Iterations < in.timeSteps
    Iterations = in.timeSteps;
    simu.T = simu.dt * Iterations;
end
% Append time with steady signal (to reach param-defined network steady state)
% + scale with base nu_ext
in.data = nu_ext.*cat(1, init.in.steadyVal.*ones(in.dim, init.strapIters), in.data);

% Initializing neuronal filters
in.tune = uniformSimplex(net.N, in.dim);
idx = kmeans(in.tune,in.dim);
[~,map] = sort(idx);

% Sort matrix W by cluster identification
W = W(map,map);

%% Initialization
% %%%%%%%%   SIMULATION PARAMETERS %%%%%%%%  

N_rp=round(neu.t_rp/simu.dt);        % Refractory (in simu.dt)
N_del=round(net.D/simu.dt);          % delay (in simu.dt)

V=neu.V_r*ones(net.N,1);            % Voltage vector
RI=zeros(net.N,N_del);
LS=-(N_del+1)*rand(net.N,1);
allspikes=zeros(1,Iterations);
%ExInput=syn.J*poissrnd(nu_ext*C_ext*simu.dt,net.N,Iterations);
Current=zeros(NCX*NCY,Iterations);

tau_pre = syn.tauDamp;
tau_post = syn.tauDamp;

Rasterplot=zeros(net.N,Iterations);

% %%%%%%% INITIALIZING GIFS %%%%%%%
LWisimu.dths = (5.*(G.Edges.Weight>0) + 0.8.*(G.Edges.Weight<0)).*abs(G.Edges.Weight);
EColors = [0 1 0].*(G.Edges.Weight>0) + [1 0 0].*(G.Edges.Weight<0);
if gif.graph
    figure(1)
    plot(G,'LineWisimu.dth',LWisimu.dths,'EdgeColor',EColors)
    axis tight
    set(gca,'nextplot','replacechildren','visible','off')
    f1 = getframe;
    [im1,map1] = rgb2ind(f1.cdata,256,'nodither');
end

if gif.lapl
    A = W + W';
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

if gif.W
    figure(3)
    imagesc(W)
    f3 = getframe;
    [im3,map3] = rgb2ind(f3.cdata,256,'nodither');
end

% Building external stimuli history
in.nu = in.tune*in.data;

%% Simulation
tic();

if init.strap > 0
    for i=1:init.strapIters
        if (1+mod(i-1,1e2))==1
            ExInput=syn.J*poissrnd(C_ext*simu.dt .* in.nu(:,i),net.N,1e2);
        end
        V=(1-simu.dt/neu.tau)*V+ExInput(:,1+mod(i-1,1e2))+RI(:,1+mod(i-1,N_del));        % Voltage update
        ca = ca.*exp(-simu.dt/syn.tauCa);
        xpre = 1 - exp(-simu.dt/tau_pre).*(1-xpre);
        xpost = 1 - exp(-simu.dt/tau_post).*(1-xpost);
        spike = (V>=neu.V_t);
        spikingNeurons = find(spike);
        ca(spikingNeurons,:) = ca(spikingNeurons,:) + syn.Cpost.*xpost(spikingNeurons,:);
        ca(:,spikingNeurons) = ca(:,spikingNeurons) + syn.Cpre.*xpre(:,spikingNeurons);
        xpost(spikingNeurons,:) = xpost(spikingNeurons,:)*(1 - syn.dampFactor);
        xpre(:,spikingNeurons) = xpre(:,spikingNeurons)*(1 - syn.dampFactor);
        rho = rho + simu.dt./syn.tauRho .* (syn.gPot.*(syn.rhoMax - rho).*actPot - syn.gDep.*rho.*actDep);
        W = synSign.*transfer(rho, prot);
        V(LS>i-N_del)=neu.V_r;
    end
end

for i=1:Iterations
    if (1+mod(i-1,1e2))==1
        ExInput=syn.J*poissrnd(C_ext*simu.dt .* in.nu(:,init.strapIters + i),net.N,1e2);
    end
    
    % nu_ext= max(0,net.rExtRel*nu_thresh*(1+5*cos(i*simu.dt/(0.5*2*pi))));
    
    V=(1-simu.dt/neu.tau)*V+ExInput(:,1+mod(i-1,1e2))+RI(:,1+mod(i-1,N_del));        % Voltage update
    ca = ca.*exp(-simu.dt/syn.tauCa);
    xpre = 1 - exp(-simu.dt/tau_pre).*(1-xpre);
    xpost = 1 - exp(-simu.dt/tau_post).*(1-xpost);
    
    Current(:,i)=DCN*V;         % Current to the electrodes
    spike = (V>=neu.V_t);               % Spiking neurons have a "1"
    spikingNeurons = find(spike);
    
    ca(spikingNeurons,:) = ca(spikingNeurons,:) + syn.Cpost.*xpost(spikingNeurons,:);
    ca(:,spikingNeurons) = ca(:,spikingNeurons) + syn.Cpre.*xpre(:,spikingNeurons);
    
    % xpost(spikingNeurons,:) = 0;
    % xpre(:,spikingNeurons) = 0;
    xpost(spikingNeurons,:) = xpost(spikingNeurons,:)*(1 - syn.dampFactor);
    xpre(:,spikingNeurons) = xpre(:,spikingNeurons)*(1 - syn.dampFactor);
    
    actPot = (ca >= syn.tPot).*(W~=0);
    actDep = (ca >= syn.tDep).*(W~=0);
    
    rho = rho + simu.dt./syn.tauRho .* (syn.gPot.*(syn.rhoMax - rho).*actPot - syn.gDep.*rho.*actDep);
    W = synSign.*transfer(rho, prot);
    
    meanWexc(i+1,1) = mean(W((1/syn.J).*W>0));
    meanWinh(i+1,1) = mean(W((1/syn.J).*W<0));
    
    histRho(:,i) = (1/net.N^2).*histcounts(rho(rho>0),edgesRho);
    histW_exc(:,i) = (1/net.N^2).*histcounts((1/syn.J).*W((1/syn.J).*W>=5e-3),edgesW_exc);
    histW_inh(:,i) = (1/net.N^2).*histcounts((1/syn.J).*W((1/syn.J).*W<=-5e-3),edgesW_inh);
    
    V(LS>i-N_del)=neu.V_r;          % Refractory period. 
    
    % Current generated
    %RI(:,1+mod(i-1,N_del))=W*spike; 
    RI(:,1+mod(i,N_del))=W*spike;
    
    LS(spike)=i;                % Time of last spike
    V(spike)=neu.V_r;              % Reset membrane potential
    
    % Store spike times
    Rasterplot(:,i)=spike;
    allspikes(1,i)=sum(spike); % Each row is (neuron number,spike time)
    
    % Updating synapse sample history
    splSyn.ca(:,i) =  ca(splSyn.IDs);
    splSyn.xpre(:,i) =  xpre(splSyn.IDs);
    splSyn.xpost(:,i) =  xpost(splSyn.IDs);
    splSyn.rho(:,i) =  rho(splSyn.IDs);
    splSyn.w(:,i) = W(splSyn.IDs);
    
    if mod(i,20)==0
        if gif.W
            figure(3)
            imagesc(W)
            f3 = getframe;
            im3(:,:,1,floor(i/20)) = rgb2ind(f3.cdata,map3,'nodither');
        end
        
        % Printing network to GIF
        if gif.graph
            G = digraph(W([subIDsExc subIDsInh],[subIDsExc subIDsInh])');
            LWisimu.dths = (5.*(G.Edges.Weight>0) + 0.8.*(G.Edges.Weight<0)).*abs(G.Edges.Weight);
            EColors = [0 1 0].*(G.Edges.Weight>0) + [1 0 0].*(G.Edges.Weight<0);
            figure(1)
            plot(G,'LineWisimu.dth',LWisimu.dths,'EdgeColor',EColors)
            f1 = getframe;
            im1(:,:,1,floor(i/20)) = rgb2ind(f1.cdata,map1,'nodither');
        end
        
        % Eigendecomposition to find clusters
        if gif.lapl
            A = abs(W + W');
            D = diag(sum(A,1));
            L_rw = eye(net.N) - D^(-1)*A;
            [eVals, ~] = eig(L_rw);
            figure(2)
            spc = sort(diag(eVals));
            plot(1:20, spc(1:20), '+r')
            xlabel('Index')
            ylabel('Eigenvalue')
            title(strcat('Eigendecomposition at time ', num2str(simu.dt*i,3),'s'))
            f2 = getframe;
            im2(:,:,1,floor(i/20)) = rgb2ind(f2.cdata,map2,'nodither');
        end
    end
    
    progressbar(i/Iterations);
    
end
toc()

%% Post-simulation analysis
G = digraph(W');
LWisimu.dths = (5.*(G.Edges.Weight>0) + 0.8.*(G.Edges.Weight<0)).*abs(G.Edges.Weight);
EColors = [0 1 0].*(G.Edges.Weight>0) + [1 0 0].*(G.Edges.Weight<0);

A = abs(W + W');
D = diag(sum(A,1));
L_rw = eye(net.N) - D^(-1)*A;
[eVals, ~] = eig(L_rw);

% filter = repmat(gausswin(100,2.5e-2),net.N,1);
filter = gausswin(100,3);
filter = (1/sum(filter)).*filter;
for n=1:net.N
    rateEst(n,:) = (1/simu.dt).*conv(Rasterplot(n,:), filter, 'same');
end
edgesRates = 0:2:200;
for i=1:Iterations
    histRates(:,i) = histcounts(rateEst(:,i),size(edgesRates,2)-1);
end


%% Creating graph for visualization
if gif.graph
    imwrite(im1, map1, strcat(env.outputsRoot, 'Figures/Network/sampleGraph_randinit.gif'), 'DelayTime',0, 'LoopCount',inf)
end

if gif.lapl
    imwrite(im2,map2, strcat(env.outputsRoot, 'Figures/Network/specClust_asym.gif'), 'DelayTime',0, 'LoopCount',inf)
end

if gif.W
    imwrite(im1, map1, strcat(env.outputsRoot, 'Figures/Network/sampleGraph_randinit.gif'), 'DelayTime',0, 'LoopCount',inf)
end

%% Plotting network stats

switch plt.all.raster
    case 1
        figure(1)
        [I1,I2] = find(Rasterplot);
        totActSnaps = sum(Rasterplot);
        ax1 = subplot(2,1,1);
        imagesc(Rasterplot)
        xticklabels(simu.dt.*xticks)
        ax2 = subplot(2,1,2);
        plot(totActSnaps)
        xticklabels(simu.dt.*xticks)
        title('Spiking activity in neural population')
        ax1.Position(1,2) = ax1.Position(1,2) - 0.7*ax2.Position(1,4);
        ax1.Position(1,4) = 1.7*ax1.Position(1,4);
        ax2.Position(1,4) = 0.3*ax2.Position(1,4);
    case 2
        rasterSnaps = zeros(net.N, plt.timeSpl.n*ceil(plt.timeSpl.dur/simu.dt) + (plt.timeSpl.n-1)*plt.timeSpl.inter);
        totActSnaps = zeros(1, plt.timeSpl.n*ceil(plt.timeSpl.dur/simu.dt) + (plt.timeSpl.n-1)*plt.timeSpl.inter);
        xticksList=[]; xtickvalsList=[];
        for tSplID = 1:plt.timeSpl.n
            rasterSnaps(:, 1+(tSplID-1)*(ceil(plt.timeSpl.dur/simu.dt)+plt.timeSpl.inter):tSplID*ceil(plt.timeSpl.dur/simu.dt)+(tSplID-1)*plt.timeSpl.inter) ...
                = Rasterplot(:,1+floor((tSplID-1)/(plt.timeSpl.n-1)*(Iterations-ceil(plt.timeSpl.dur/simu.dt)-1)):floor((tSplID-1)/(plt.timeSpl.n-1)*(Iterations-ceil(plt.timeSpl.dur/simu.dt)-1)+ceil(plt.timeSpl.dur/simu.dt)));           

            rasterSnaps(:, 1 + tSplID*ceil(plt.timeSpl.dur/simu.dt)+(tSplID-1)*plt.timeSpl.inter : plt.timeSpl.inter + tSplID*ceil(plt.timeSpl.dur/simu.dt)+(tSplID-1)*plt.timeSpl.inter) ... 
                = 1;
            
            xticksList = cat(2, xticksList, ...
                1+(tSplID-1)*(ceil(plt.timeSpl.dur/simu.dt)+plt.timeSpl.inter) + ceil((1/3)*(plt.timeSpl.dur/simu.dt).*(0:1:3))...
                );
            
            xtickvalsList = cat(2, xtickvalsList, ...
                simu.dt.*(1+floor((tSplID-1)/(plt.timeSpl.n-1)*(Iterations-ceil(plt.timeSpl.dur/simu.dt)-1)) + ceil((1/3)*(plt.timeSpl.dur/simu.dt).*(0:1:3)))...
                );
        end    
        totActSnaps = sum(rasterSnaps);
        for tSplID = 1:plt.timeSpl.n
            totActSnaps(1, 1 + tSplID*ceil(plt.timeSpl.dur/simu.dt)+(tSplID-1)*plt.timeSpl.inter : plt.timeSpl.inter + tSplID*ceil(plt.timeSpl.dur/simu.dt)+(tSplID-1)*plt.timeSpl.inter) = 0;
        end
        figure(1)
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
end

if plt.spl.ca
    figure(2)
    imagesc(splSyn.ca)
    title('Synaptic calcium actvivity in sample synapses')
    colorbar
end

if plt.spl.rho
    figure(3)
    imagesc(splSyn.rho)
    title('Phosphorylation state at sample synapses')
    colorbar
end

if plt.spl.w
    splSyn.wHist = synSign(splSyn.IDs).*transfer(splSyn.rho, prot);
    figure(4)
    imagesc(splSyn.wHist)
    title('Synaptic weight at sample synapses')
    colorbar
end

splSyn.stats = cat(2, linspace(1,splSyn.n,splSyn.n)', splSyn.PreNeurons, splSyn.PostNeurons, splSyn.w(:,1), splSyn.w(:,end));

switch plt.spl.pres
    case 1
        splSyn.pres = zeros(3*splSyn.n, Iterations);
        for i=1:splSyn.n
            splSyn.pres(3*(i-1)+1,:) = 100*Rasterplot(splSyn.PreNeurons(i,1),:);
            splSyn.pres(3*(i-1)+2,:) = 100*Rasterplot(splSyn.PostNeurons(i,1),:);
            splSyn.pres(3*(i-1)+3,:)= (1/syn.J).*synSign(splSyn.IDs(i)).*splSyn.rho(i,:);
        end
        figure(5)
        imagesc(splSyn.pres)
        colorbar
    case 2
        splSyn.pres = zeros(5*floor(splSyn.n/2), Iterations);
        for i=1:splSyn.n
            splSyn.pres(5*(i-1)+1,:) = Rasterplot(splSyn.PreNeurons(i,1),:);
            splSyn.pres(5*(i-1)+2,:) = Rasterplot(splSyn.PostNeurons(i,1),:);
            splSyn.pres(5*(i-1)+3,:)= splSyn.ca(i,:);
            splSyn.pres(5*(i-1)+4,:)= (1/(syn.J*syn.rhoMax)).*synSign(splSyn.IDs(i)).*splSyn.rho(i,:);
            splSyn.pres(5*(i-1)+5,:)= (1/syn.J).*splSyn.w(i,:);
        end
        figure(5)
        imagesc(splSyn.pres)
        xticklabels(simu.dt.*xticks)
        colorbar
    case 3
        splSyn.pres = zeros(5*floor(splSyn.n/2), plt.timeSpl.n*ceil(plt.timeSpl.dur/simu.dt) + (plt.timeSpl.n-1)*plt.timeSpl.inter);
        xticksList=[]; xtickvalsList=[];
        for tSplID = 1:plt.timeSpl.n
            for i=1:splSyn.n
                splSyn.pres(5*(i-1)+1, 1+(tSplID-1)*(ceil(plt.timeSpl.dur/simu.dt)+plt.timeSpl.inter):tSplID*ceil(plt.timeSpl.dur/simu.dt)+(tSplID-1)*plt.timeSpl.inter) ...
                    = Rasterplot(splSyn.PreNeurons(i,1),1+floor((tSplID-1)/(plt.timeSpl.n-1)*(Iterations-ceil(plt.timeSpl.dur/simu.dt)-1)):floor((tSplID-1)/(plt.timeSpl.n-1)*(Iterations-ceil(plt.timeSpl.dur/simu.dt)-1)+ceil(plt.timeSpl.dur/simu.dt)));           

                splSyn.pres(5*(i-1)+2, 1+(tSplID-1)*(ceil(plt.timeSpl.dur/simu.dt)+plt.timeSpl.inter):tSplID*ceil(plt.timeSpl.dur/simu.dt)+(tSplID-1)*plt.timeSpl.inter) ... 
                    = Rasterplot(splSyn.PostNeurons(i,1),1+floor((tSplID-1)/(plt.timeSpl.n-1)*(Iterations-ceil(plt.timeSpl.dur/simu.dt)-1)):floor((tSplID-1)/(plt.timeSpl.n-1)*(Iterations-ceil(plt.timeSpl.dur/simu.dt)-1)+ceil(plt.timeSpl.dur/simu.dt)));

                splSyn.pres(5*(i-1)+3, 1+(tSplID-1)*(ceil(plt.timeSpl.dur/simu.dt)+plt.timeSpl.inter):tSplID*ceil(plt.timeSpl.dur/simu.dt)+(tSplID-1)*plt.timeSpl.inter) ...
                    = splSyn.ca(i,1+floor((tSplID-1)/(plt.timeSpl.n-1)*(Iterations-ceil(plt.timeSpl.dur/simu.dt)-1)):floor((tSplID-1)/(plt.timeSpl.n-1)*(Iterations-ceil(plt.timeSpl.dur/simu.dt)-1)+ceil(plt.timeSpl.dur/simu.dt)));

                splSyn.pres(5*(i-1)+4, 1+(tSplID-1)*(ceil(plt.timeSpl.dur/simu.dt)+plt.timeSpl.inter):tSplID*ceil(plt.timeSpl.dur/simu.dt)+(tSplID-1)*plt.timeSpl.inter) ...
                    = (1/(syn.J*syn.rhoMax)).*synSign(splSyn.IDs(i)).*splSyn.rho(i,1+floor((tSplID-1)/(plt.timeSpl.n-1)*(Iterations-ceil(plt.timeSpl.dur/simu.dt)-1)):floor((tSplID-1)/(plt.timeSpl.n-1)*(Iterations-ceil(plt.timeSpl.dur/simu.dt)-1)+ceil(plt.timeSpl.dur/simu.dt)));

                splSyn.pres(5*(i-1)+5, 1+(tSplID-1)*(ceil(plt.timeSpl.dur/simu.dt)+plt.timeSpl.inter):tSplID*ceil(plt.timeSpl.dur/simu.dt)+(tSplID-1)*plt.timeSpl.inter) ...
                    = (1/syn.J).*splSyn.w(i,1+floor((tSplID-1)/(plt.timeSpl.n-1)*(Iterations-ceil(plt.timeSpl.dur/simu.dt)-1)):floor((tSplID-1)/(plt.timeSpl.n-1)*(Iterations-ceil(plt.timeSpl.dur/simu.dt)-1)+ceil(plt.timeSpl.dur/simu.dt)));
                
            end
            
            xticksList = cat(2, xticksList, ...
                1+(tSplID-1)*(ceil(plt.timeSpl.dur/simu.dt)+plt.timeSpl.inter) + ceil((1/3)*(plt.timeSpl.dur/simu.dt).*(0:1:3))...
                );
            
            xtickvalsList = cat(2, xtickvalsList, ...
                simu.dt.*(1+floor((tSplID-1)/(plt.timeSpl.n-1)*(Iterations-ceil(plt.timeSpl.dur/simu.dt)-1)) + ceil((1/3)*(plt.timeSpl.dur/simu.dt).*(0:1:3)))...
                );
        end
        figure(5)
        imagesc(splSyn.pres)
        xticks(xticksList)
        xticklabels(xtickvalsList)
        colorbar
end

if plt.spl.hist
    tickList = ceil(nbins.*(0.1:0.1:1));
    figure(6)
    ax1 = subplot(2,1,1);
    imagesc(log(histW_exc));
    set(gca,'YDir','normal')
    title('Evolution of weight distribution for excitatory synapses')
    xlabel('Time')
    xticklabels(simu.dt.*xticks)
    ylabel('Synaptic weight')
    yticks(tickList)
    yticklabels(edgesW_exc(1, tickList))
    colorbar
    Wexc_hStep = edgesW_exc(2) - edgesW_exc(1);
    valsWexc = edgesW_exc + Wexc_hStep;
    valsWexc = repmat(valsWexc(:,1:end-1)',1,Iterations);
    sumWexc = sum(valsWexc.*histW_exc);
    ax2 = subplot(2,1,2);
    bar(edgesW_exc(2:end), histW_exc(:,end))
    ax1.Position(1,2) = ax1.Position(1,2) - 0.5*ax2.Position(1,4);
    ax1.Position(1,4) = 1.5*ax1.Position(1,4);
    ax2.Position(1,4) = 0.5*ax2.Position(1,4);
    
    
    figure(7)
    ax1 = subplot(2,1,1);
    imagesc(log(histW_inh));
    set(gca,'YDir','normal')
    title('Evolution of weight distribution for inhibitory synapses')
    xlabel('Time')
    xticklabels(simu.dt.*xticks)
    ylabel('Synaptic weight')
    yticks(tickList)
    yticklabels(edgesW_inh(1, tickList))
    colorbar
    Winh_hStep = edgesW_inh(2) - edgesW_inh(1);
    valsWinh = edgesW_inh + Winh_hStep;
    valsWinh = repmat(valsWinh(:,1:end-1)',1,Iterations);
    sumWinh = sum(valsWinh.*histW_inh);
    ax2 = subplot(2,1,2);
    bar(edgesW_inh(2:end), histW_inh(:,end))
    ax1.Position(1,2) = ax1.Position(1,2) - 0.5*ax2.Position(1,4);
    ax1.Position(1,4) = 1.5*ax1.Position(1,4);
    ax2.Position(1,4) = 0.5*ax2.Position(1,4);
    
    
    figure(8)
    ax1 = subplot(2,1,1);
    imagesc(log(histRho));
    set(gca,'YDir','normal')
    title('Evolution of phosphorylation distribution (all synapses)')
    xlabel('Time')
    xticklabels(simu.dt.*xticks)
    ylabel('Phosphorylation level')
    yticks(tickList)
    yticklabels(edgesRho(1, tickList))
    colorbar
    ax2 = subplot(2,1,2);
    bar(edgesRho(2:end), histRho(:,end))
    ax1.Position(1,2) = ax1.Position(1,2) - 0.5*ax2.Position(1,4);
    ax1.Position(1,4) = 1.5*ax1.Position(1,4);
    ax2.Position(1,4) = 0.5*ax2.Position(1,4);
    
    
    % sprintf('Excitatory: %0.3f - Inhibitory: %0.3f - Total: %0.3f', sumWexc(1,end)/syn.J, -sumWinh(1,end)/syn.J, (sumWexc(1,end) + sumWinh(1,end))/syn.J)
    figure(9)
    ax1 = subplot(2,1,1);
    imagesc(log(histRates));
    set(gca,'YDir','normal')
    for n=1:splNeu.NE
        hold on
        plot((nbins/edgesRates(1,end)).*rateEst(splNeu.ExcNeurons(n,1),:), 'g', 'LineWidth', 1)
    end
    for n=1:splNeu.NI
        hold on
        plot((nbins/edgesRates(1,end)).*rateEst(splNeu.InhNeurons(n,1),:), 'r', 'LineWidth', 1)
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
    case 1
        
        figure(10)
        c = rand(1,3);
        plot((1/meanWexc(2,1)).*meanWexc(2:end,1), (1/meanWinh(2,1)).*meanWinh(2:end,1), '-', 'Color', c)
        for i=1:20:201
            text((1/meanWexc(2,1)).*meanWexc(1+i,1), (1/meanWinh(2,1)).*meanWinh(i+1,1),num2str(1+i),'Color', c)
        end
        for i=301:100:1001
            text((1/meanWexc(2,1)).*meanWexc(1+i,1), (1/meanWinh(2,1)).*meanWinh(i+1,1),num2str(1+i),'Color', c)
        end
        for i=2001:1000:30001
            text((1/meanWexc(2,1)).*meanWexc(1+i,1), (1/meanWinh(2,1)).*meanWinh(i+1,1),num2str(1+i),'Color', c)
        end
        title('Evolution of average synaptic weights')
        xlabel('Relative mean excitatory weight')
        ylabel('Relative mean inhibitory weight')
    case 2
        figure(10)
        c = rand(1,3);
        plot((1/meanWexc(2,1)).*meanWexc(2:end,1), -meanWinh(2:end,1)./meanWexc(2:end,1), '-', 'Color', c)
        for i=1:20:201
            text((1/meanWexc(2,1)).*meanWexc(1+i,1), -meanWinh(1+i,1)./meanWexc(1+i,1),num2str(1+i),'Color', c)
        end
        for i=301:100:1001
            text((1/meanWexc(2,1)).*meanWexc(1+i,1), -meanWinh(1+i,1)./meanWexc(1+i,1),num2str(1+i),'Color', c)
        end
        for i=2001:1000:20001
            text((1/meanWexc(2,1)).*meanWexc(1+i,1), -meanWinh(1+i,1)./meanWexc(1+i,1),num2str(1+i),'Color', c)
        end
        title('Evolution of average synaptic weights')
        xlabel('Relative mean excitatory weight')
        ylabel('Relative mean inhibitory weight')
        % ylim([min(3, min(-meanWinh(2:end,1)./meanWexc(2:end,1))) max(8, max(-meanWinh(2:end,1)./meanWexc(2:end,1)))])

        SRtoAI_thr = refline([0 4.2]);
        SRtoAI_thr.Color = 'r';

        AItoSI_thr = refline([0 5.3]);
        AItoSI_thr.Color = 'g';
end