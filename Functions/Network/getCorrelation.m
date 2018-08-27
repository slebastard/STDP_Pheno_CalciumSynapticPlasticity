function [outputArg1,outputArg2] = getCorrelation(pre_spikes_hist, post_spikes_hist, params)
%GETCORRELATION Summary of this function goes here
%   Detailed explanation goes here
% Components of params required
% dt: bin size
% T: experiment duration

T = params.T;
dt = params.dt;
% Find empirical rates
pre.rate = length(pre_spikes_hist)/pre_spikes_hist(end,1);
post.rate = length(post_spikes_hist)/post_spikes_hist(end,1);

% Computing correlation
for k = floor(-T/(3*dt)):1:ceil(T/(3*dt))
    
end

end

