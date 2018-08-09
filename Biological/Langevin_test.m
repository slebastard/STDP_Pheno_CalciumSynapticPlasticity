Ca_stim = 1;
Ca_bas= 0.1;
tau=1;

t0=0;
dt=0.05;
tstim=1;
tfinal=5;

y0 = Ca_stim;

eq_det1 = @(t,y) 0;
eq_det2 = @(t,y) -(1/tau)*(y-Ca_bas);
eq_stoch = @(t,y) 0.1;
 
opts = sdeset('NonNegative',1,...
              'RandSeed',1,...
              'SDEType','Ito');

t = t0:dt:tfinal;
y = sde_euler(eq_det1, eq_stoch, t0:dt:tstim, y0', opts);
y = [y; sde_euler(eq_det2, eq_stoch, tstim+dt:dt:tfinal, y(end,:), opts)];

plot(t, y);