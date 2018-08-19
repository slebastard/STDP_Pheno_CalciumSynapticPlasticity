function plotSynapseHistory(dt, Duration, Rasterplot, splSynapses, synID, damping)
    figure(6)
    t = 0:dt:Duration-dt;
    Iterations=ceil(Duration/dt);   
    
    subplot(6,1,1)
    [Ipre1,Ipre2] = find(Rasterplot(splSynapses.PreNeurons(synID),:));
    plot(Ipre2, Ipre1,'.','MarkerSize',4)
    xlim([1 Iterations])
    ylim([0.75 1.25])
    % plot(t, splSynapses.PreNeurons(synID,:))
    ylabel('Presyn spikes')    

    subplot(6,1,2)
    [Ipost1,Ipost2] = find(Rasterplot(splSynapses.PostNeurons(synID),:));
    plot(Ipost2, Ipost1,'.','MarkerSize',4)
    xlim([1 Iterations])
    ylim([0.75 1.25])
    % plot(t, splSynapses.PreNeurons(synID,:))
    ylabel('Postsyn spikes')   
    
    if damping
        subplot(6,1,3)
        plot(t, splSynapses.xpre(synID,:))
        ylabel('Ca presyn sat var')

        subplot(6,1,4)
        plot(t, splSynapses.xpost(synID,:))
        ylabel('Ca postsyn sat var')         
        
        subplot(6,1,5)
        plot(t, splSynapses.ca(synID,:))
        ylabel('Calcium concentration')

        subplot(6,1,6)
        plot(t, splSynapses.rho(synID,:))
        ylabel('Phosphorylation level')         
    else
        subplot(6,1,3)
        plot(t, splSynapses.ca(synID,:))
        ylabel('Calcium concentration')

        subplot(6,1,4)
        plot(t, splSynapses.rho(synID,:))
        ylabel('Phosphorylation level')        
    end
    
    xlabel('Time (s)')
end