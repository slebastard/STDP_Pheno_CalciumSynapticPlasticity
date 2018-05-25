% Loading data

data = csvread('data_MSN_Simon.csv',1,1);

% Col       Field               Unit
% 1    1    dt                  ms
% 2    2    STDP                %rel%
% 4    3    Init EPSP ampl      mV
% 6    4    Final EPSP ampl     mV
% 7    5    jitter              ms
% 9    6    Plasticity test     %cat%

data = data(:,[1 2 4 6 7 9]);

% First try: 3D plot of dEPSP/EPSP = f(dt, EPSP_0)

figure(1)
scatter3(data(:,1), data(:,3), data(:,4))
xlabel('dt')
ylabel('Init EPSP ampl')
zlabel('Final EPSP ampl')
% zlabel('Rel var in EPSP ampl')