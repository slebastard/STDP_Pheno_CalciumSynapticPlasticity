function wf = zeta(rho, sigma)
%ZETA Summary of this function goes here
%   Detailed explanation goes here

S_attr = 1; % this is a function of theta_alpha

if nargin == 1
    sigma = 1; % replace with default sigma from G&B 2007
end

wf = 0.5*erfc((S_attr-rho)/sqrt(2*sigma^2));

end

