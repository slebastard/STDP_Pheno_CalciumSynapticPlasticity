function met = getRegimeMetrics(simu, plt, fig)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    fig.regime = figure('Name','NET_Regime_Metrics','NumberTitle','off');
    ax1 = subplot(2,1,1);
    plot(plt.synMetr(:,1));
    synch_thr = refline([0 plt.regimeThr(1,1)]);
    synch_thr.Color = 'r';
    simu()
    title('Synchrony metric')
    xlabel('Iteration')
    ylabel('StdDev of phase transform')
    xticklabels(simu.dt.*xticks)
    
    ax2 = subplot(2,1,2);
    plot(plt.regMetr(:,1));
    reg_thr = refline([0 plt.regimeThr(1,2)]);
    reg_thr.Color = 'r';
    title('Regularity metric')
    xlabel('Iteration')
    ylabel('Coefficient of variation of ISI')
    xticklabels(simu.dt.*xticks)
    
end

