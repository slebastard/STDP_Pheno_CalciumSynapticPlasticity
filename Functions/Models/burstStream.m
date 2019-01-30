function out = burstStream(signal, tuning, T, dt, res)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

N = size(tuning,2);
K = floor(T/dt);

out = zeros(K,N);

for i=1:K
    out(i,:) = tuning(1+floor(res*signal(i)/pi),:);
end

end

