function [t, I] = corrPoisson( N, mu, C, T )
%CORRPOISSON Generates correlated Poisson processes
% Currently only supports homogeneous Poisson processes

% %% Inputs %%
%   N  - Number of Poisson processes to simulate
%   mu - Contains the spike rates of all processes
%     N*1 matrix
%   C  - Correlation matrix
%     N*N matrix
%   T  - Total duration of simulation

% %% Outputs %%
%   t  - Spike times
%   I  - N*length(t) matrix containing spike flags at each times, for each
%   of the N processes to simulate

    tc = 1; % Time constant for Ornstein-Uhlenbeck

    L = chol(C, 'lower');
    M = sum(mu) + 10*sum(C)^(0.5);
    A = zeros(N,N);
    b = zeros(N,1);
    for i=N:1
        A(i:end,:) = A(i:end,:) + repmat(L(i,:),N+1-i,1);
        b(i:end,1) = b(i:end,1) + mu(i,1);
    end
    r = sum(mu);
    
    t = 0;
    Y = zeros(N, 1);
    I = zeros(N, 1);
    
    while t(1,end) < T
        dur = exprnd(1/M);
        t = cat(2, t, t(1,end)+dur);
        Y = Y.*exp(-dur/tc) + sqrt(1 - exp(-2*dur/tc)).*rand(N,1);
        c = b + A*Y;
        s = r + A*Y;
        if M*rand() <= s
            v = s*rand();
            % Binary search
            
            left = 1;
            right = N;
            index = -1;
            while index==-1
                if right==1
                    index=1;
                    break
                elseif left==N
                    index = N;
                    break
                end
                
                mid = ceil((left + right) / 2);
                if c(mid,1) < v
                    left = mid;
                elseif c(mid-1,1) > v
                    right = mid-1;
                else
                    index = mid;
                end
            end
            spike =  zeros(N,1);
            spike(index,1) = 1;
            I = cat(2, I, spike);
        end
    end
    
end

