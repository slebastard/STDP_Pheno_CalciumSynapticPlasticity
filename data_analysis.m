function data_analysis()
%% Loading data

data = csvread('data_MSN_Simon.csv',1,1);

% Col       Field               Unit
% 1    1    dt                  ms
% 2    2    STDP                %rel%
% 4    3    Init EPSP ampl      mV
% 6    4    Final EPSP ampl     mV
% 7    5    jitter              ms
% 9    6    Plasticity test     %cat%

data = data(:,[1 2 4 6 7 9]);

% Available modes:
%    EPSPf to visualize final EPSP amplitude as a function of init EPSP
%  amplitude and dt
%    relSTDP to vizualise the amplitude ratio EPSPf/EPSPi as a function of
%  init EPSP amplitude and dt
mode = 'relSTDP';

% First try: 3D plot of dEPSP/EPSP = f(dt, EPSP_0)

figure(1)
if strcmp(mode, 'EPSPf')
    scatter3(data(:,1), data(:,3), data(:,4))
elseif strcmp(mode, 'relSTDP')
    scatter3(data(:,1), data(:,3), data(:,2))
end
hold on

%% Implementing the analytical formulation predicting STDP from naive model

C_pre = 1;
C_post = 2;
theta_dep = 1;
gamma_dep = 200;
theta_pot = 1.3;
gamma_pot = 321;

freq = 1;
tau = 150000;
tau_Ca = 20;


n_dt = 50;

dt_min = min(data(:,1));
dt_max = max(data(:,1));

STDP_naive = [];

for dt = data(:,1)
    [a, b] = STDP(dt);
    STDP_naive = cat(1, STDP_naive, [dt, a, b]);
end

lc_pot = tau_Ca * log((theta_pot - C_pre)/C_post);
hc_pot = tau_Ca * log(theta_pot/C_post);
lc_dep = tau_Ca * log((theta_dep - C_pre)/C_post);
hc_dep = tau_Ca * log(theta_dep/C_post);

XL = get(gca, 'XLim');
x_min = XL(1);

YL = get(gca, 'YLim');
y_min = YL(1);
y_max = YL(2);

ZL = get(gca, 'ZLim');
z_min = ZL(1);
z_max = ZL(2);

%Low dep - Low pot
patch([lc_dep, lc_dep, lc_pot, lc_pot], [y_min, y_max, y_max, y_min], [0, 0, 0, 0], 'FaceColor', [0.5 0 0.5]);

%Low pot - High dep
patch([lc_pot, lc_pot, hc_dep, hc_dep], [y_min, y_max, y_max, y_min], [0, 0, 0, 0], 'FaceColor', [0.5 0.5 0]);

%High dep - High pot
patch([hc_dep, hc_dep, hc_pot, hc_pot], [y_min, y_max, y_max, y_min], [0, 0, 0, 0], 'FaceColor', [0 0.5 0.5]);
alpha(0.2)

if strcmp(mode, 'EPSPf')
    scatter3(data(:,1), data(:,3), STDP_naive(:,2).*data(:,3) + STDP_naive(:,3), 'r+')
    xlabel('dt')
    ylabel('Init EPSP ampl')
    zlabel('Final EPSP ampl')
elseif strcmp(mode, 'relSTDP')
    scatter3(data(:,1), data(:,3), STDP_naive(:,2) + STDP_naive(:,3) ./ data(:,3), 'r+')
    xlabel('dt')
    ylabel('Init EPSP ampl')
    zlabel('Rel var in EPSP ampl')
end

hold off

figure(2)
plot(data(:,1), a, 'x')
title('Slope as a fct of dt')
xlabel('dt')
ylabel('Slope')

figure(3)
plot(data(:,1), Ca_topTheta_rate(theta_pot, data(:,1)), 'x')
title('Rate pot as a fct of dt')
xlabel('dt')
ylabel('r_{pot}')

figure(4)
plot(data(:,1), b./(1-a), 'x')
title('Limit efficacy as a function of \delta_t')
xlabel('\delta_t')
ylabel('Limit efficacy')

%% Functions definition
function [r, dt_crit_low, dt_crit_high] = Ca_topTheta_rate(theta, dt)
    dt_crit_low = log((theta - C_pre)/C_post);
    dt_crit_high = log(theta/C_post);

    r = tau_Ca * (freq/1000) * (...
        log((C_post * exp(dt/tau_Ca) + C_pre)./(theta*exp(dt/tau_Ca))) .* (dt/tau_Ca > dt_crit_high) ...
        + (log(C_post/theta) + log((C_post*exp(dt/tau_Ca) + C_pre)/theta)) .* (dt/tau_Ca > dt_crit_low) .* (dt/tau_Ca <= dt_crit_high) ...
        + log(C_post/theta) .* (dt/tau_Ca <= dt_crit_low) ...
        );   
end

function [a, b] = STDP(dt)
    r_pot= Ca_topTheta_rate(theta_pot, dt);
    r_dep = Ca_topTheta_rate(theta_dep, dt) - r_pot;
    
    a = exp(-(r_pot*(gamma_pot+gamma_dep)+r_dep*gamma_dep)/(tau*freq/1000));
    b = gamma_pot/(gamma_pot + gamma_dep)*(exp(-(gamma_dep*r_dep)/(tau*freq/1000))-exp(-(gamma_pot+gamma_dep)*r_pot/(tau*freq/1000)));
end

end