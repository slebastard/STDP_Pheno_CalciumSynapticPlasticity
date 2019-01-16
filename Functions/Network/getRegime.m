function plt = getRegime(plt, step)
% getRegime Summary of this function goes here
%   Regimes description:
% 1 - SR
% 2 - AR
% 3 - SI
% 4 - AI

plt.regime(step,1) = 0;

% MEASURING REGULARITY
% - clustering on 1d Fourier signals?
% - finding modes on each row, using the number of modes
% - using the width of the fft on each row, and averaging these number over
% all rows to get an index of irregularity
%  * Do we need to know the number of clusters?
% MEASURING SYNCHRONICITY
% - Don't just capture fft, as phase is the most important aspect here.
% Therefore find a way to capture the phase
% Then compute a std deviation of phase for each time step, in a way that
% takes into account 2*pi periodicity for phase
% 
% - 
y = fftshift()

end

