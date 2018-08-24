function t = indPoisson( N, mu, T )
%CORRPOISSON Generates uncorrelated Poisson processes
% Currently only supports homogeneous Poisson processes

% %% Inputs %%
%   N  - Number of Poisson processes to simulate
%   mu - Contains the spike rates of all processes
%     N*1 matrix
%   T  - Total duration of simulation

% %% Outputs %%
%   t  - Spike times
%   I  - N*length(t) matrix containing spike flags at each times, for each
%   of the N processes to simulate
    
    t = zeros(N, 1);  % Contains the min of time reached by generation
    
    while min(t(:,end)) < T
        intsp = randexp(mu);
        t = cat(2, t, t(:,end) + intsp);
    end
    
end

