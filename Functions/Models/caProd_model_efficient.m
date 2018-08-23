function [rho_hist, w_end, c_hist] = caProd_model_efficient( pre_spikes_hist, post_spikes_hist, params, simu )
%NAIVE_MODEL Simulates the behavior of a synapse whose behavior follows the
%naive Calcium_based dynamics
%   Detailed explanation goes here

%% Unpacking params
%%%%%%%%%%%%%%%%%%%


switch nargin
    case 0
        pre_spikes_hist = [];
        post_spikes_hist = [];
        params = default_params();
        simu = default.simu();
    case 1
        post_spikes_hist = [];
        params = default_params();
        simu = default.simu();
    case 2
        params = default_params();
        simu = default.simu();
    case 3
        simu = default.simu();
    case 4
    otherwise
        error('4 inputs max are accepted')
end

%%%%%%%%%%%%%%%%%%%%
% Unpacking params %
%%%%%%%%%%%%%%%%%%%%

T = simu.T;
step = simu.int_step;
n_steps = T / step;

rho_0 = params.rho_0;
rho_max = params.rho_max;
C_pre = params.C_pre;
C_post = params.C_post;
tau_Ca = params.tau_Ca;
tau_x = params.tau_x;
delay_pre = params.delay_pre;
theta_dep = params.theta_dep;
gamma_dep = params.gamma_dep;
theta_pot = params.theta_pot;
gamma_pot = params.gamma_pot;
tau_rho = params.tau_rho;
sigma = params.noise_lvl;
S_attr = params.S_attr;
tau_w = params.tau_w;
theta_act = params.theta_act;
w_0 = params.w_0;

tau_CaPre = tau_Ca;
tau_CaPost = tau_Ca;
dampFactor = params.dampFactor;

prot = params;
prot.frequency = 1000/T;
prot.n_iter = 1;

%% Simulating process
%%%%%%%%%%%%%%%%%%%%%

% Initiate simulation %
%%%%%%%%%%%%%%%%%%%%%%%
rho = rho_0;
w = w_0;
t = 0;
x_pre = 1;
x_post = 1;
c_pre = 0;
c_post = 0;
c = 0;

rho_hist = rho;
c_hist = [];
c_histPre = [];
c_histPost = [];
times = [];

evts = [];

len_pre = size(pre_spikes_hist,2);
len_post = size(post_spikes_hist,2);

if len_pre > 0
    for pre_evtID=1:len_pre
        evts = cat(1, evts, [pre_spikes_hist(pre_evtID) + delay_pre, -1]);
    end
end

if len_post > 0
    for post_evtID=1:len_post
        evts = cat(1, evts, [post_spikes_hist(post_evtID),1]);
    end
end

evts = sortrows(evts, 1);

% Check whether simulation is trivial %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isempty(evts)
    rho_hist = cat(1, [0, rho_0], [T, rho_0]);
    w_end = w_0;
    return
end
    
% Figure out the calcium history %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for bump_id = 1:size(evts,1)
    tn = evts(bump_id, 1);
    times = [times; tn];
    x_pre = 1 - exp(-(tn-t)/tau_x)*(1 - x_pre); % comment this out in prod model
    x_post = 1 - exp(-(tn-t)/tau_x)*(1 - x_post); % comment this out in prod model
    c_post = C_post.*(x_post.^x_pre).*(evts(bump_id,2)==1) + c_post.*exp(-(tn-t)/tau_CaPost);
    c_pre = C_pre.*(x_pre.^x_post).*(evts(bump_id,2)==-1) + c_pre.*exp(-(tn-t)/tau_CaPre);
    x_pre = x_pre.*(1 - (evts(bump_id,2)==-1).*dampFactor);
    x_post = x_post.*(1 - (evts(bump_id,2)==1)*dampFactor);
    % PRODUCT MODEL
    % c = c_pre.*c_post;
    % ADDITIVE MODELS
    c = c_pre + c_post;
    
    c_hist = [c_hist;c];
    t = tn;
end
c_hist = cat(2, times, c_hist);

% Extract period objects for which the calcium is above thresholds %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if any(c_hist(:,2)>theta_pot)
    t_openPot = evts(c_hist(:,2)>theta_pot, 1);
    c_openPot = c_hist(c_hist(:,2)>theta_pot, 2);
    t_maxDurPot = tau_Ca.*log(c_openPot/theta_pot);
    t_closePot = t_openPot + min(circshift(t_openPot, -1) - t_openPot, tau_Ca.*log(c_openPot./theta_pot));
    t_closePot(end) =  t_openPot(end) + t_maxDurPot(end);
    t_closePot = cat(2, t_openPot, t_closePot, ones(length(t_closePot),1));
end

if any(c_hist(:,2)>theta_dep)
    t_openDep = evts(c_hist(:,2)>theta_dep,1);
    c_openDep = c_hist(c_hist(:,2)>theta_dep, 2);
    t_maxDurDep = tau_Ca.*log(c_openDep/theta_dep);
    t_closeDep = t_openDep + min(circshift(t_openDep,-1) - t_openDep, tau_Ca.*log(c_openDep./theta_dep));
    t_closeDep(end) =  t_openDep(end) + t_maxDurDep(end);
    t_closeDep = cat(2, t_openDep, t_closeDep, -1.*ones(length(t_closeDep),1));
end

if ~(exist('t_closePot', 'var') || exist('t_closeDep', 'var'))
    rho_hist = cat(1, [0, rho_0], [T, rho_0]);
    w_end = w_0;
    return
elseif ~exist('t_closePot', 'var')
    t_closePot = [0.0 0.0 1];
elseif ~exist('t_closeDep', 'var')
    t_closeDep = [0.0 0.0 -1];
end

t_dur = cat(1, t_closePot, t_closeDep);
t_dur = sortrows(t_dur, [1 2]);

t_durOverlap = (circshift(t_dur(:,1),1) - t_dur(:,1) == 0);
t_durOverlap(1) = 0;
substitute = circshift(t_dur(:,2),1);
t_dur(t_durOverlap,1) = substitute(t_durOverlap);

% Finding rho %
%%%%%%%%%%%%%%%%%%%%%%

for id=1:size(t_dur,1)
    rho_f = ((rho - rho_max*gamma_pot/(gamma_pot+gamma_dep))*exp(-(t_dur(id,2)-t_dur(id,1))*(gamma_pot+gamma_dep)/tau_rho) + rho_max*gamma_pot/(gamma_pot+gamma_dep)) .* (t_dur(id,3)==1) ...
        + rho.* exp(-gamma_dep*(t_dur(id,2)-t_dur(id,1))/tau_rho) .*(t_dur(id,3)==-1);
    rho = rho_f;
    rho_hist = cat(1, rho_hist, rho);
end

intvls_pot = t_dur(t_dur(:,3)==1,:);
dur_pot = sum(intvls_pot(:,2) - intvls_pot(:,1));
prot.alpha_pot = dur_pot/T;

intvls_dep = t_dur(t_dur(:,3)==-1,:);
dur_dep = dur_pot + sum(intvls_dep(:,2) - intvls_dep(:,1));    
prot.alpha_dep = dur_dep/T;

times = cat(1, t_dur(:,1), t_dur(end,2));
rho_hist = cat(2, times, rho_hist);
w_end = transfer(rho_hist(end), prot);

end    
