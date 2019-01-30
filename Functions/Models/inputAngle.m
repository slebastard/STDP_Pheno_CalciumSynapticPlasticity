function out = inputAngle(sigma, T ,dt)
%INPUTANGLE Summary of this function goes here
%   Detailed explanation goes here

i = 1;
out = 0;
for t=dt:dt:T
    i = i+1;
    r1 = randn(1);
    out(i) = mod(out(i-1) + sigma*r1, pi);
end

end

