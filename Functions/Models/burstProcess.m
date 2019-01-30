function out = burstProcess(tau, c, T)
%BURSTPROCESS Generates N independent inhomogeneous Poisson processes,
% according to a OU rate
%
% out.proc contains the spike times
% out.rate contains the rate function
%
% All times in ms
% All amplitudes in mV
%
% tau is the relaxation time constant
% c is the diffussion constant

dt = 0.01; % integrate step

r0 = 0; % initial value for rate process
itg0 = 0;

% 1) Implement an Ornstein-Uhlenbeck process
i = 1;
out.rate = r0;
for t=dt:dt:T
    i = i+1;
    r1 = randn(1);
    out.rate(i) = out.rate(i-1)*exp(-dt/tau) + sqrt((c*tau*0.5)*(1-(exp(-dt/tau))^2))*r1;
end

%2) Redress this process to have a smooth, positive function of time
out.rate(out.rate<0)=0;
maxRate = max(max(out.rate));

% Find integral or redressed signal
% i = 1;
% itg(1) = itg0;
% for t=dt:dt:T
%     i = i+1;
%     itg(i) = itg(i-1) + rate(i)*dt;
% end

t = 0;
i=1;

while t < T
    t = t - log(rand)./maxRate;
    if t > T
        break
    end
    if rand <= out.rate(floor(t/dt))/maxRate
        out.proc(i) = t;
        i = i+1;
    end
end

end

