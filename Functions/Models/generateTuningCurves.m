function out = generateTuningCurves(sigma,N,res,maxRate,minRate)
%GENERATETUNINGCURVES Summary of this function goes here
%   Detailed explanation goes here
out = zeros(res,N);
for i=1:N
    mean = pi*rand;
    angles = pi/res:pi/res:pi;
    out(:,i) = minRate + circshift(maxRate.*exp(-((angles-pi/2).^2)./(2*sigma^2)),floor(res*(mean-pi/2)/pi));
end

end

