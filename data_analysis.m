function data_analysis()
%% Loading data

data = csvread('data_MSN.csv',1,0);

% CSV col   Array col   Field               Unit
% 1         1           sourceID            Laurent=1, Elodie=2, Hao=3, Yihui=4
% 2         2           dt                  ms
% 3                     n_pairs
% 4                     frequency
% 5         3           STDP                %rel%
% 6         4           Init EPSP ampl      mV
% 7         5           Final EPSP ampl     mV

data = data(:,[1 2 5 6 7]);

% Data mapping
source = data(:,1);
dt = data(:,2);
r_rho = data(:,3);
w_i = data(:,4);
w_f = data(:,5);

% Experimental values of efficacy need to be normalized first

w_max = [];
w_min = [];
for i=1:4
    w_max = cat(1, w_max, max(cat(1, w_i .* (source==i), w_f .* (source==i))));
    w_min = cat(1, w_min, min(cat(1, w_i .* (source==i) + 100000 .* (source~=i), w_f .* (source==i) + 100000 .* (source~=i))));
end
rho_i = (w_i - w_min(source))./(w_max(source) - w_min(source));
rho_f = (w_f - w_min(source))./(w_max(source) - w_min(source));

% Available modes:
%    EPSPf to visualize final EPSP amplitude as a function of init EPSP
%  amplitude and dt
%    relSTDP to vizualise the amplitude ratio EPSPf/EPSPi as a function of
%  init EPSP amplitude and dt
mode = 'EPSPf';

% First try: 3D plot of dEPSP/EPSP = f(dt, EPSP_0)

figure(1)
if strcmp(mode, 'EPSPf')
    N = size(dt,1);
    colorMap = [zeros(N, 1), zeros(N, 1), ones(N,1)];
    pos_dt = (dt >= 0);
    colorMap(pos_dt, :) = repmat([1,0,0],nnz(pos_dt),1);
    %scatter3(dt, rho_i, (rho_f-rho_i)./(rho_i.^1.5), 8, colorMap);
    scatter3(dt, rho_i, rho_f, 8, colorMap);
    view(90,0)
%     hold on
%     [x,y] = meshgrid(-30:30, 0:1);
%     surf(x,y,0.85*y-0.18)
%     zlim([0,1])
%     shading interp;
%     alpha(0.5)
elseif strcmp(mode, 'relSTDP')
    N = size(dt,1);
    colorMap = [zeros(N, 1), zeros(N, 1), ones(N,1)];
    pos_dt = (dt >= 0);
    colorMap(pos_dt, :) = repmat([1,0,0],nnz(pos_dt),1);
    scatter3(dt, rho_i, r_rho, 8, colorMap); 
end
hold on

%% Implementing the analytical formulation predicting STDP from naive model

C_pre = 1.1;
C_post = 0.7;
theta_dep = 1;
gamma_dep = 2;
theta_pot = 1.3;
gamma_pot = 3.2;

freq = 1;
n_iter = 100;
tau = 150000;
tau_Ca = 20;
delay_pre = 5;

n_dt = 50;

dt_min = min(dt);
dt_max = max(dt);

[a, b] = STDP(dt);
STDP_naive = cat(2, dt, a, b);

lc_pot = tau_Ca * log((theta_pot - C_pre)/C_post);
hc_pot = tau_Ca * log(C_pre/(theta_pot-C_post));
lc_dep = tau_Ca * log(C_pre/theta_dep);
hc_dep = tau_Ca * log(C_pre/(theta_dep - C_post));

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

rho_lim = STDP_naive(:,3) ./ (1-STDP_naive(:,2));
rho_lim(isnan(rho_lim)) = rho_i(isnan(rho_lim));
rho_f_model = (rho_i - rho_lim).*STDP_naive(:,2).^(n_iter) + rho_lim;
r_rho_model = rho_f_model ./ rho_i;

if strcmp(mode, 'EPSPf')
    colorMap = [zeros(N, 1), ones(N, 1), zeros(N,1)];
    %scatter3(dt, rho_i, (rho_f_model - rho_i)./rho_i, 8, colorMap)
    scatter3(dt, rho_i, rho_f_model, 8, colorMap)
    xlabel('dt')
    ylabel('Init EPSC ampl')
    zlabel('Final EPSC ampl')
elseif strcmp(mode, 'relSTDP')
    scatter3(dt, rho_i, r_rho_model, 8, 'g')
    xlabel('dt')
    ylabel('Init EPSC ampl')
    zlabel('Rel var in EPSC ampl')
end

hold off

figure(2)
plot(dt, a, 'x')
title('Slope as a fct of dt')
xlabel('dt')
ylabel('Slope')
% 
% figure(3)
% plot(dt, Ca_topTheta_rate(theta_pot, dt), 'x')
% title('Rate pot as a fct of dt')
% xlabel('dt')
% ylabel('r_{pot}')
% 
% figure(4)
% plot(dt, b./(1-a), 'x')
% title('Limit efficacy as a function of \delta_t')
% xlabel('\delta_t')
% ylabel('Limit efficacy')

%% Functions definition
function r = Ca_topTheta_rate(theta, dt)

    if C_pre<=theta && C_post<=theta
        dt_crit_low = log((theta-C_pre)/C_post);
        dt_crit_high = log(C_pre/(theta-C_post));
        r = tau_Ca * (freq/1000) * (...
            log((C_pre+C_post*exp(dt/tau_Ca))/theta) .* (dt/tau_Ca  > dt_crit_low) .* (dt/tau_Ca  <= 0) ...
            + log((C_post+C_pre*exp(-dt/tau_Ca))/theta) .* (dt/tau_Ca  > 0) .* (dt/tau_Ca  <= dt_crit_high) ...
        );

    elseif theta<=C_pre && theta<=C_post
        dt_crit_low = log(theta/C_post);
        dt_crit_high = log(C_pre/theta);

        r = tau_Ca * (freq/1000) * (...
            ( log(C_pre/theta) + log((C_pre*exp(-dt/tau_Ca) + C_post)/theta) ) .* (dt/tau_Ca > dt_crit_high) ...
            +  (log((C_pre+C_post*exp(dt/tau_Ca))/theta)-(dt/tau_Ca)) .* (dt/tau_Ca  > dt_crit_low) .* (dt/tau_Ca  <= dt_crit_high) ...
            + ( log(C_post/theta) + log((C_post*exp(dt/tau_Ca) + C_pre)/theta) ) .* (dt/tau_Ca  <= dt_crit_low) ...
            );

    elseif C_pre<theta
        dt_crit_low = log((theta - C_pre)/C_post);
        dt_crit_high = log(theta/C_post);

        r = tau_Ca * (freq/1000) * (...
            log((C_post * exp(dt/tau_Ca) + C_pre)./(theta*exp(dt/tau_Ca))) .* (dt/tau_Ca > dt_crit_high) ...
            + (log(C_post/theta) + log((C_post*exp(dt/tau_Ca) + C_pre)/theta)) .* (dt/tau_Ca  > dt_crit_low) .* (dt/tau_Ca  <= dt_crit_high) ...
            + log(C_post/theta) .* (dt/tau_Ca  <= dt_crit_low) ...
            );

    else
        dt_crit_low = log(C_pre/theta);
        dt_crit_high = log(C_pre/(theta-C_post));

        r = tau_Ca * (freq/1000) * (...
            log((C_pre * exp(-dt/tau_Ca) + C_post)./(theta*exp(-dt/tau_Ca))) .* (dt/tau_Ca  <= dt_crit_low) ...
            + (log(C_pre/theta) + log((C_pre*exp(-dt/tau_Ca) + C_post)/theta)) .* (dt/tau_Ca  > dt_crit_low) .* (dt/tau_Ca  <= dt_crit_high) ...
            + log(C_pre/theta) .* (dt/tau_Ca > dt_crit_high) ...
            );

    end
end

function [a, b] = STDP(dt)
    r_pot= Ca_topTheta_rate(theta_pot, dt-delay_pre);
    r_dep = Ca_topTheta_rate(theta_dep, dt-delay_pre) - r_pot;
    
    a = exp(-(r_dep*gamma_dep + r_pot*(gamma_dep+gamma_pot))/((freq/1000)*tau));
    b = (gamma_pot/(gamma_pot + gamma_dep)) * exp(-(r_dep*gamma_dep)/(tau*(freq/1000))) .* (1 - exp(-(r_pot*(gamma_pot+gamma_dep))/(tau*(freq/1000))));
end

function rho_f = predict(rho_i, dt)
    [a,b] = STDP(dt);
    rho_lim = b./(1-a);
    rho_f = (rho_i-rho_lim)*a.^n_iter + rho_lim;
end

end