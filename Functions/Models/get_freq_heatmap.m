function heatmap = get_freq_heatmap(dataFit, params)
% STDP EXPERIMENT
% - Runs a battery of model simulation with Calcium bumps
% reflecting different temporal differences. Uses those simulation to build
% a dw = f(delta_t) curve
% - All time params are in ms, all frequencies are in Hz
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Default parameter values + unpacking params
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch nargin
    case 0
        error('Please specify parameters in dataFit object')
    case 1
        params = default_params();
    case 2
    otherwise
        error('2 inputs max are accepted. Please provide freq and dt parameters as arrays')
end

%%%%%%%%%%%%%%%%%%%%
% Unpacking params %
%%%%%%%%%%%%%%%%%%%%

STDP = dataFit;

%% Running simulations, returning STDP curve
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n_points_freq = 1+floor((dataFit.freq.max-1)/dataFit.freq.step);
freqs = linspace(1, dataFit.freq.max, n_points_freq);

heatmap = [];

for freq_id = 1:n_points_freq
    frq = freqs(1,freq_id);
    STDP.frequency = frq;
    if strcmp(dataFit.model, 'caProd')
        std = get_STDP_CaProd(STDP, params);
    else
        std = get_STDP(STDP, params);
    end
    std = cat(2, frq*ones(size(std,1),1), std);
    heatmap = cat(1, heatmap, std);
end

end