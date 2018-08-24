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
addpath(genpath('Functions'))


%% Parameterization

% %%%%%%%%   PARAMETERS OF THE SIMULATION  %%%%%%%% 

syn = get_synapse();
prot.TD = 0;
prot.S_attr = syn.sAttr;
prot.noise_lvl = syn.sigma;

Duration=5;
dt=5e-4;
N=50;
Iterations=ceil(Duration/dt);

plt.all.raster = 0;
plt.spl.ca = 0;
plt.spl.rho = 0;
plt.spl.w = 0;
plt.spl.pres = 2;
plt.spl.hist = 1;

gif.graph = 0;
gif.lapl = 0;

init.w = 'norm';
init.c = 0.5;

% %%%%%%%%   PARAMETERS OF THE NETWORK  %%%%%%%%   

NE=ceil(0.8*N);      % # excitatory neurons
NI=N-NE;             % # excitatory neurons

Connectivity=0.4;           % Connectivity coefficients, parameter epsilon in Brunel 2000.
D=3e-3;                        % Transmission Delay
CE=round(Connectivity*NE);          % Number of exc connections
CI=round(Connectivity*NI);          % Number of inh connections
C_ext=CE;


% %%%%%%%%   PARAMETERS OF THE NEURONS  %%%%%%%%   

% Threshold, reset, time constant and refractory period (p.185, 2nd column,
% Brunel 2000)
V_t=20e-3;
V_r=10e-3;
tau=20e-3;
t_rp=2e-3;

% Bifurcation parameters
J=0.2e-3; %0.1e-3;                   % Strength of exc. Connections
g=8;                      % Strength of inh/exc (inh connections -gJ)
ratioextthresh=0.9;         % nu_ext/nu_thresh

nu_thresh=V_t/(J*CE*tau);   % Frequency needed for a neuron to reach the threshold. 
nu_ext=ratioextthresh*nu_thresh;       % external Poisson input rate


% %%%%%%%%   PARAMETERS OF THE ELECTRODES  %%%%%%%%  
L=100;          % Size of the cortex
NCX=8;          % Number of electrodes rows
NCY=8;          % Number of electrodes columns
NC=NCX*NCY;     % Total number of electrodes

% Parameters of histograms
nbins = 100;
edgesRho = linspace(0,syn.rhoMax,nbins+1);
edgesW_exc = linspace(0,J,nbins+1);
edgesW_inh = linspace(-g*J,0,nbins+1);

%% Creating network

% %%%%%% INITIALIZING SYNAPTIC QTIES %%%%%%

% %%% Synaptic weights %%%
W = zeros(N,N);   % 1st index for postsynaptic neuron, 2nd for presynaptic

for i=1:N
    ECells=randperm(NE);
    ECells=ECells(1:CE);
    W(i,ECells)=init.c .* J;
    % W(i,ECells)=J.*rand(1,CE);
    
    ICells=randperm(NI);
    ICells=NE+ICells(1:CI);
    W(i,ICells)=-init.c .*g*J;
    % W(i,ICells)=-g*J.*rand(1,CI);
end

synSign = J.*(W>0) - g*J.*(W<0);


G = digraph(W');    % Graph for plotting evolution of network
A = abs(W + W'); % Weight matrix for spectral clustering
% Symmetrization, all information on directionnality is lost
histRho = zeros(nbins, Iterations);
histW_exc = zeros(nbins, Iterations);
histW_inh = zeros(nbins, Iterations);

% %%% Other qties %%%
ca = zeros(N,N);
xpre = ones(N,N);
xpost = ones(N,N);
% rho = zeros(N,N);
rho = transferinv(W./synSign, syn.sAttr, syn.sigma, syn.rhoMax);
actPot = zeros(N,N);
actDep = zeros(N,N);


% %%% Electrodes wiring %%%
% Electrodes are regularly located, we compute the attenuation coefficient
% due to the distance. 

[X,Y] = meshgrid((1:8)',(1:8)');
POSC =10* [X(:) Y(:)];              % Positions of electrodes
POSN=L*rand(N,2);                   % Positions of the neurons

DCN=zeros(NCX*NCY,N);
for i=1:(NCX*NCY)
    for j=1:N
        DCN(i,j)=sqrt((POSC(i,1)-POSN(j,1))^2+(POSC(i,2)-POSN(j,2))^2);
    end
end
DCN=(DCN.^(-2)).*(DCN<10);
%W=rand(N,N)<=Connectivity;


% Sample synapses
splSynapses.n = 15;
[synOut, synIn] = find(W);
perm = randperm(length(synOut));
rpOut = synOut(perm);
rpIn = synIn(perm);

splSynapses.PostNeurons = rpOut(1:splSynapses.n);
splSynapses.PreNeurons = rpIn(1:splSynapses.n);
splSynapses.IDs = sub2ind(size(ca),splSynapses.PostNeurons,splSynapses.PreNeurons);
splSynapses.ca = zeros(splSynapses.n,Iterations);
splSynapses.xpre = zeros(splSynapses.n,Iterations);
splSynapses.xpost = zeros(splSynapses.n,Iterations);
splSynapses.rho = zeros(splSynapses.n,Iterations);
splSynapses.w = zeros(splSynapses.n,Iterations);

%% Initialization
% %%%%%%%%   SIMULATION PARAMETERS %%%%%%%%  

N_rp=round(t_rp/dt);        % Refractory (in dt)
N_del=round(D/dt);          % delay (in dt)

V=V_r*ones(N,1);            % Voltage vector
RI=zeros(N,N_del);
LS=-(N_del+1)*rand(N,1);
allspikes=zeros(1,Iterations);
%ExInput=J*poissrnd(nu_ext*C_ext*dt,N,Iterations);
Current=zeros(NCX*NCY,Iterations);

tau_pre = syn.tauDamp;
tau_post = syn.tauDamp;

Rasterplot=zeros(N,Iterations);

% %%%%%%% INITIALIZING GIFS %%%%%%%
LWidths = (5.*(G.Edges.Weight>0) + 0.8.*(G.Edges.Weight<0)).*abs(G.Edges.Weight);
EColors = [0 1 0].*(G.Edges.Weight>0) + [1 0 0].*(G.Edges.Weight<0);
if gif.graph
    figure(1)
    plot(G,'LineWidth',LWidths,'EdgeColor',EColors)
    axis tight
    set(gca,'nextplot','replacechildren','visible','off')
    f1 = getframe;
    [im1,map1] = rgb2ind(f1.cdata,256,'nodither');
end

if gif.lapl
    A = W + W';
    D = diag(sum(A,1));
    L_rw = eye(N) - D^(-1)*A;
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
for i=1:Iterations
    if (1+mod(i-1,1e4))==1
        ExInput=J*poissrnd(nu_ext*C_ext*dt,N,1e4);
    end
    
    V=(1-dt/tau)*V+ExInput(:,1+mod(i-1,1e4))+RI(:,1+mod(i-1,N_del));        % Voltage update
    ca = ca.*exp(-dt/syn.tauCa);
    xpre = 1 - exp(-dt/tau_pre).*(1-xpre);
    xpost = 1 - exp(-dt/tau_post).*(1-xpost);
    
    Current(:,i)=DCN*V;         % Current to the electrodes
    spike = (V>=V_t);               % Spiking neurons have a "1"
    spikingNeurons = find(spike);
    
    ca(spikingNeurons,:) = ca(spikingNeurons,:) + syn.Cpost.*xpost(spikingNeurons,:);
    ca(:,spikingNeurons) = ca(:,spikingNeurons) + syn.Cpre.*xpre(:,spikingNeurons);
    
    % xpost(spikingNeurons,:) = 0;
    % xpre(:,spikingNeurons) = 0;
    xpost(spikingNeurons,:) = xpost(spikingNeurons,:)*(1 - syn.dampFactor);
    xpre(:,spikingNeurons) = xpre(:,spikingNeurons)*(1 - syn.dampFactor);
    
    actPot = (ca >= syn.tPot).*(W~=0);
    actDep = (ca >= syn.tDep).*(W~=0);
    
    % W = syn.plast(W, spike, dt);
    rho = rho + dt./syn.tauRho .* (syn.gPot.*(syn.rhoMax - rho).*actPot - syn.gDep.*rho.*actDep);
    W = synSign.*transfer(rho, prot);
    
    
    histRho(:,i) = (1/N^2).*histcounts(rho,edgesRho);
    histW_exc(:,i) = (1/N^2).*histcounts(W(W>=0),edgesW_exc);
    histW_inh(:,i) = (1/N^2).*histcounts(W(W<=0),edgesW_inh);
    
    V(LS>i-N_del)=V_r;          % Refractory period. 
    
    % Current generated
    %RI(:,1+mod(i-1,N_del))=W*spike; 
    RI(:,1+mod(i,N_del))=W*spike;
    
    LS(spike)=i;                % Time of last spike
    V(spike)=V_r;              % Reset membrane potential
    
    % Store spike times
    Rasterplot(:,i)=spike;
    allspikes(1,i)=sum(spike); % Each row is (neuron number,spike time)
    
    % Updating synapse sample history
    splSynapses.ca(:,i) =  ca(splSynapses.IDs);
    splSynapses.xpre(:,i) =  xpre(splSynapses.IDs);
    splSynapses.xpost(:,i) =  xpost(splSynapses.IDs);
    splSynapses.rho(:,i) =  rho(splSynapses.IDs);
    splSynapses.w(:,i) = W(splSynapses.IDs);
    
    if mod(i,20)==0
        % Printing network to GIF
        if gif.graph
            G = digraph(W');
            LWidths = (5.*(G.Edges.Weight>0) + 0.8.*(G.Edges.Weight<0)).*abs(G.Edges.Weight);
            EColors = [0 1 0].*(G.Edges.Weight>0) + [1 0 0].*(G.Edges.Weight<0);
            figure(1)
            plot(G,'LineWidth',LWidths,'EdgeColor',EColors)
            f1 = getframe;
            im1(:,:,1,floor(i/20)) = rgb2ind(f1.cdata,map1,'nodither');
        end
        
        % Eigendecomposition to find clusters
        if gif.lapl
            A = abs(W + W');
            D = diag(sum(A,1));
            L_rw = eye(N) - D^(-1)*A;
            [eVals, ~] = eig(L_rw);
            figure(2)
            spc = sort(diag(eVals));
            plot(1:20, spc(1:20), '+r')
            xlabel('Index')
            ylabel('Eigenvalue')
            title(strcat('Eigendecomposition at time ', num2str(dt*i,3),'s'))
            f2 = getframe;
            im2(:,:,1,floor(i/20)) = rgb2ind(f2.cdata,map2,'nodither');
        end
    end
    
    progressbar(i/Iterations);
    
end
toc()
beep;

G = digraph(W');
LWidths = (5.*(G.Edges.Weight>0) + 0.8.*(G.Edges.Weight<0)).*abs(G.Edges.Weight);
EColors = [0 1 0].*(G.Edges.Weight>0) + [1 0 0].*(G.Edges.Weight<0);

A = abs(W + W');
D = diag(sum(A,1));
L_rw = eye(N) - D^(-1)*A;
[eVals, ~] = eig(L_rw);

%% Correlation analysis
splSynapses.corr = zeros();

%% Creating graph for visualization
if gif.graph
    imwrite(im1,map1,'Outputs/Figures/Network/sampleGraph_randinit.gif','DelayTime',0,'LoopCount',inf)
end

if gif.lapl
    imwrite(im2,map2,'Outputs/Figures/Network/specClust_asym.gif','DelayTime',0,'LoopCount',inf)
end

%% Plotting network stats

if plt.all.raster
    figure()
    % subplot(2,1,1);
    [I1,I2] = find(Rasterplot);
    plot(I2,I1,'.','MarkerSize',1)
    title('Spiking activity in neural population')
end

if plt.spl.ca
    figure()
    imagesc(splSynapses.ca)
    title('Synaptic calcium actvivity in sample synapses')
    colorbar
end

if plt.spl.rho
    figure()
    imagesc(splSynapses.rho)
    title('Phosphorylation state at sample synapses')
    colorbar
end

if plt.spl.w
    splSynapses.wHist = synSign(splSynapses.IDs).*transfer(splSynapses.rho, prot);
    figure()
    imagesc(splSynapses.wHist)
    title('Synaptic weight at sample synapses')
    colorbar
end

splSynapses.stats = cat(2, linspace(1,splSynapses.n,splSynapses.n)', splSynapses.PreNeurons, splSynapses.PostNeurons, splSynapses.w(:,1), splSynapses.w(:,end))

if plt.spl.pres == 1
    splSynapses.pres = zeros(3*splSynapses.n, Iterations);
    for i=1:splSynapses.n
        splSynapses.pres(3*(i-1)+1,:) = 100*Rasterplot(splSynapses.PreNeurons(i,1),:);
        splSynapses.pres(3*(i-1)+2,:) = 100*Rasterplot(splSynapses.PostNeurons(i,1),:);
        splSynapses.pres(3*(i-1)+3,:)= (1/J).*synSign(splSynapses.IDs(i)).*splSynapses.rho(i,:);
    end
    figure()
    imagesc(splSynapses.pres)
    colorbar
elseif plt.spl.pres == 2
    splSynapses.pres = zeros(5*floor(splSynapses.n/2), Iterations);
    for i=1:splSynapses.n
        splSynapses.pres(5*(i-1)+1,:) = Rasterplot(splSynapses.PreNeurons(i,1),:);
        splSynapses.pres(5*(i-1)+2,:) = Rasterplot(splSynapses.PostNeurons(i,1),:);
        splSynapses.pres(5*(i-1)+3,:)= splSynapses.ca(i,:);
        splSynapses.pres(5*(i-1)+4,:)= (1/(J*syn.rhoMax)).*synSign(splSynapses.IDs(i)).*splSynapses.rho(i,:);
        splSynapses.pres(5*(i-1)+5,:)= (1/J).*splSynapses.w(i,:);
    end
    figure()
    imagesc(splSynapses.pres)
    colorbar
end

if plt.spl.hist    
    figure()
    imagesc(log(histW_exc));
    set(gca,'YDir','normal')
    title('Evolution of weight distribution for excitatory synapses')
    xlabel('Time')
    ylabel('Synaptic weight')
    yticks(tickList)
    yticklabels(edgesW_exc(1, tickList))
    colormap(bluewhitered)
    colorbar
    Wexc_hStep = edgesW_exc(2) - edgesW_exc(1);
    valsWexc = edgesW_exc + Wexc_hStep;
    valsWexc = repmat(valsWexc(:,1:end-1)',1,Iterations);
    sumWexc = sum(valsWexc.*histW_exc);
    
    figure()
    imagesc(log(histW_inh));
    set(gca,'YDir','normal')
    title('Evolution of weight distribution for inhibitory synapses')
    xlabel('Time')
    ylabel('Synaptic weight')
    yticks(tickList)
    yticklabels(edgesW_inh(1, tickList))
    colormap(bluewhitered)
    colorbar
    Winh_hStep = edgesW_inh(2) - edgesW_inh(1);
    valsWinh = edgesW_inh + Winh_hStep;
    valsWinh = repmat(valsWinh(:,1:end-1)',1,Iterations);
    sumWinh = sum(valsWinh.*histW_inh);
    
    figure()
    imagesc(log(histRho));
    set(gca,'YDir','normal')
    title('Evolution of phosphorylation distribution (all synapses)')
    xlabel('Time')
    ylabel('Phosphorylation level')
    tickList = [10 20 30 40 50 60 70 80 90 100];
    yticks(tickList)
    yticklabels(edgesRho(1, tickList))
    colormap(bluewhitered)
    colorbar
    
    figure()
    imagesc(repmat(sumWexc + sumWinh, 8, 1))
    colormap(bluewhitered)
    colorbar
    sprintf('Excitatory: %0.3f - Inhibitory: %0.3f - Total: %0.3f', sumWexc(1,end)/J, -sumWinh(1,end)/J, (sumWexc(1,end) + sumWinh(1,end))/J)
end