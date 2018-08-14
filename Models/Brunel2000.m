% Simulation of the Brunel JCNS 2000 network with electrode recording
% spikes (for paper Destexhe Touboul on Criticality).

clear all
close all

% Add plot tools to path
addpath('PlotTools')

% Parameters to adjust for the simulation:

Duration=20;
dt=0.1e-3;
N=1000;

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



W=zeros(N,N);

for i=1:N
    ECells=randperm(NE);
    ECells=ECells(1:CE);
    W(i,ECells)=J;
    
    ICells=randperm(NI);
    ICells=NE+ICells(1:CI);
    W(i,ICells)=-g*J;
end
    

% Electrodes are regularly located, we compute the attenuation coefficient
% due to the distance. 

[X,Y] = meshgrid((1:8)',(1:8)');
POSC =10* [X(:) Y(:)];              % Positions of electrodes
POSN=L*rand(N,2);                   % Position of the neurons

DCN=zeros(NCX*NCY,N);
for i=1:(NCX*NCY)
    for j=1:N
        DCN(i,j)=sqrt((POSC(i,1)-POSN(j,1))^2+(POSC(i,2)-POSN(j,2))^2);
    end
end
DCN=(DCN.^(-2)).*(DCN<10);
%W=rand(N,N)<=Connectivity;




% %%%%%%%%   SIMULATION PARAMETERS %%%%%%%%  

N_rp=round(t_rp/dt);        % Refractory (in dt)
N_del=round(D/dt);          % delay (in dt)

V=V_r*ones(N,1);            % Voltage vector
RI=zeros(N,N_del);
LS=-(N_del+1)*rand(N,1);
allspikes=zeros(1,Iterations);
%ExInput=J*poissrnd(nu_ext*C_ext*dt,N,Iterations);
Current=zeros(NCX*NCY,Iterations);

Rasterplot=zeros(N,Iterations);
tic();
for i=1:Iterations
    if (1+mod(i-1,1e4))==1
        ExInput=J*poissrnd(nu_ext*C_ext*dt,N,1e4);
    end
    
    V=(1-dt/tau)*V+ExInput(:,1+mod(i-1,1e4))+RI(:,1+mod(i-1,N_del));        % Voltage update
    Current(:,i)=DCN*V;         % Current to the electrodes
    spike=V>=V_t;               % Spiking neurons have a "1"
    
    V(LS>i-N_del)=V_r;          % Refractory period. 
    
    % Current generated
    %RI(:,1+mod(i-1,N_del))=W*spike; 
    RI(:,1+mod(i,N_del))=W*spike;
    
    LS(spike)=i;                % Time of last spike
    V(spike)=V_r;              % Reset membrane potential
    
    % Store spike times
    Rasterplot(:,i)=spike;
    allspikes(1,i)=sum(spike); % Each row is (neuron number,spike time)
    
    
    progressbar(i/Iterations);
    
end
toc()
beep;

%%
figure(1)
% subplot(2,1,1);
[I1,I2]=find(Rasterplot);
plot(I2,I1,'.','MarkerSize',1)
%~

% imagesc(Rasterplot)
% colormap (1-gray)
% figure(2)
% plot(Cumulate(allspikes,3));
% Analyse_Avalanches(allspikes,8)
% Collapse_var_Auto(allspikes,15)
