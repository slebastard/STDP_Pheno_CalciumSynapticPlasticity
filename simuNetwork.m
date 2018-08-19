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

% Parameters to adjust for the simulation:

syn = get_synapse();

Duration=1;
dt=5e-4;
N=20;

Iterations=ceil(Duration/dt);


% %%%%%%%%   PARAMETERS OF THE NETWORK  %%%%%%%%   

NE=ceil(0.8*N);      % # excitatory neurons
NI=N-NE;             % # excitatory neurons



% Threshold, reset, time constant and refractory period (p.185, 2nd column,
% Brunel 2000)
V_t=20e-3;
V_r=10e-3;
tau=20e-3;
t_rp=2e-3;



% %%%%%% BIFURCATION PARAMETERS %%%%%%%
J=0.1e-3;                   % Strength of exc. Connections
g=8;                      % Strength of inh/exc (inh connections -gJ)
ratioextthresh=0.9;         % nu_ext/nu_thresh

Connectivity=0.4;           % Connectivity coefficients, parameter epsilon in Brunel 2000.
D=3e-3;                        % Transmission Delay
CE=round(Connectivity*NE);          % Number of exc connections
CI=round(Connectivity*NI);          % Number of inh connections
C_ext=CE;


nu_thresh=V_t/(J*CE*tau);   % Frequency needed for a neuron to reach the threshold. 
nu_ext=ratioextthresh*nu_thresh;       % external Poisson input rate






% %%%%%%%%   PARAMETERS OF THE ELECTRODES  %%%%%%%%  
L=100;          % Size of the cortex
NCX=8;          % Number of electrodes rows
NCY=8;          % Number of electrodes columns
NC=NCX*NCY;     % Total number of electrodes



W=zeros(N,N);   % 1st index for postsynaptic neuron, 2nd for presynaptic

for i=1:N
    ECells=randperm(NE);
    ECells=ECells(1:CE);
    W(i,ECells)=J;
    % W(i,ECells)=J*(1+rand);
    
    ICells=randperm(NI);
    ICells=NE+ICells(1:CI);
    W(i,ICells)=-g*J;
    % W(i,ICells)=-g*J*rand;
end

synSign = (W>0) - g.*(W<0);

% Graph for plotting evolution of network
G = digraph(W');

% Weight matrix for spectral clustering
A = abs(W + W'); % Symmetrization, all information on directionnality is lost

ca = zeros(N,N);
xpre = ones(N,N);
xpost = ones(N,N);
rho = zeros(N,N);
%rho = transferinv(1e3*J.*abs(synSign), syn.sAttr, syn.sigma, syn.rhoMax).*ones(N,N);
actPot = zeros(N,N);
actDep = zeros(N,N);

nbins = 100;
histRho = zeros(nbins, Iterations);

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

splNeurons.n = 50;
splNeurons.IDs = rand(N,splNeurons.n,1);
splNeurons.V = zeros();

splSynapses.n = 15;
[synOut, synIn] = find(W);
perm = randperm(length(synOut));
rpOut = synOut(perm);
rpIn = synIn(perm);

splSynapses.PostNeurons = rpOut(1:splSynapses.n);
splSynapses.PreNeurons = rpIn(1:splSynapses.n);
splSynapses.IDs = sub2ind(size(ca),splSynapses.PostNeurons,splSynapses.PreNeurons);
splSynapses.ca = ca(splSynapses.IDs);
splSynapses.xpre = xpre(splSynapses.IDs);
splSynapses.xpost = xpost(splSynapses.IDs);
splSynapses.rho = rho(splSynapses.IDs);

% %%%%%%%%   SIMULATION PARAMETERS %%%%%%%%  

N_rp=round(t_rp/dt);        % Refractory (in dt)
N_del=round(D/dt);          % delay (in dt)

V=V_r*ones(N,1);            % Voltage vector
RI=zeros(N,N_del);
LS=-(N_del+1)*rand(N,1);
allspikes=zeros(1,Iterations);
%ExInput=J*poissrnd(nu_ext*C_ext*dt,N,Iterations);
Current=zeros(NCX*NCY,Iterations);

tau_pre = 3*syn.tauCa;
tau_post = 3*syn.tauCa;

Rasterplot=zeros(N,Iterations);
tic();

LWidths = (5.*(G.Edges.Weight>0) + 0.8.*(G.Edges.Weight<0)).*abs(G.Edges.Weight);
EColors = [0 1 0].*(G.Edges.Weight>0) + [1 0 0].*(G.Edges.Weight<0);
figure(1)
plot(G,'LineWidth',LWidths,'EdgeColor',EColors)
axis tight

set(gca,'nextplot','replacechildren','visible','off')
f1 = getframe;
[im1,map1] = rgb2ind(f1.cdata,256,'nodither');

figure(2)
A = W + W';
D = diag(sum(A,1));
L_rw = eye(N) - D^(-1)*A;
[eVals, ~] = eig(L_rw);
spc = sort(diag(eVals));
plot(1:20, spc(1:20), '+r')
xlabel('Index')
ylabel('Eigenvalue')
title('Eigendecomposition at time 0s')

set(gca,'nextplot','replacechildren')
f2 = getframe;
[im2,map2] = rgb2ind(f2.cdata,256,'nodither');

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
    xpost(spikingNeurons,:) = xpost(spikingNeurons,:) - 0.3*xpost(spikingNeurons,:);
    xpre(:,spikingNeurons) = xpre(:,spikingNeurons) - 0.3*xpre(:,spikingNeurons);
    
    
    actPot = (ca >= syn.tPot).*(W~=0);
    actDep = (ca >= syn.tDep).*(W~=0);
    
    % W = syn.plast(W, spike, dt);
    rho = rho + dt./syn.tauRho .* (syn.gPot.*(syn.rhoMax - rho).*actPot - syn.gDep.*rho.*actDep);
    W = synSign.*transfer(rho, syn.sAttr, syn.sigma);
    
    
    histRho(:,i) = (1/N^2).*histcounts(rho,nbins);
    
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
    
    if mod(i,20)==0
        % Printing network to GIF
        G = digraph(W');
        LWidths = (5.*(G.Edges.Weight>0) + 0.8.*(G.Edges.Weight<0)).*abs(G.Edges.Weight);
        EColors = [0 1 0].*(G.Edges.Weight>0) + [1 0 0].*(G.Edges.Weight<0);
        figure(1)
        plot(G,'LineWidth',LWidths,'EdgeColor',EColors)
        f1 = getframe;
        im1(:,:,1,floor(i/20)) = rgb2ind(f1.cdata,map1,'nodither');
        
        % Eigendecomposition to find clusters
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
    progressbar(i/Iterations);
    
end
toc()
beep;

%% Correlation analysis
splSynapses.corr = zeros();

%% Creating graph for visualization
% imwrite(im1,map1,'sampleGraph_randinit.gif','DelayTime',0,'LoopCount',inf)
imwrite(im2,map2,'specClust_asym.gif','DelayTime',0,'LoopCount',inf)

%% Plotting network stats

figure(1)
% subplot(2,1,1);
[I1,I2] = find(Rasterplot);
plot(I2,I1,'.','MarkerSize',1)
title('Spiking activity in neural population')

figure(2)
imagesc(splSynapses.ca)
title('Synaptic calcium actvivity in sample synapses')
colorbar

figure(3)
imagesc(splSynapses.rho)
title('Phosphorylation state at sample synapses')
colorbar

splSynapses.wHist = synSign(splSynapses.IDs).*transfer(splSynapses.rho, syn.sAttr, syn.sigma);
figure(4)
imagesc(splSynapses.wHist)
title('Synaptic weight at sample synapses')
colorbar
%

splSynapses.stats = cat(2, linspace(1,splSynapses.n,splSynapses.n)', splSynapses.PreNeurons, splSynapses.PostNeurons, splSynapses.wHist(:,1), splSynapses.wHist(:,end))

splSynapses.pres = zeros(3*splSynapses.n, Iterations);
for i=1:splSynapses.n
    splSynapses.pres(3*(i-1)+1,:) = 100*Rasterplot(splSynapses.PreNeurons(i,1),:);
    splSynapses.pres(3*(i-1)+2,:) = 100*Rasterplot(splSynapses.PostNeurons(i,1),:);
    splSynapses.pres(3*(i-1)+3,:)= synSign(splSynapses.IDs(i)).*splSynapses.rho(i,:);
end

figure(5)
imagesc(splSynapses.pres)
colorbar

figure(7)
imagesc(histRho);
colorbar

% imagesc(Rasterplot)
% colormap (1-gray)
% figure(2)
% plot(Cumulate(allspikes,3));
% Analyse_Avalanches(allspikes,8)
% Collapse_var_Auto(allspikes,15)