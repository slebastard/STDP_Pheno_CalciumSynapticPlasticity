function wf = transfer(rho, S_attr, sigma)
%ZETA Summary of this function goes here
%   Detailed explanation goes here

if nargin == 1
    S_attr = 1;
    sigma = 0.05; % replace with default sigma from G&B 2007
end

wf = 0.5*erfc((S_attr-real(rho))./sqrt(2*sigma^2));

end

