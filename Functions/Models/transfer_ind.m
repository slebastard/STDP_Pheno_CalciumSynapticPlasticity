function wf = transfer_ind(rho, protocParams)
%ZETA Summary of this function goes here
%   Detailed explanation goes here

S_attr = protocParams.S_attr;
sigma = protocParams.noise_lvl ./ 12;

wf = 0.5*erfc((S_attr-real(rho))./sqrt(2*sigma^2));

end

