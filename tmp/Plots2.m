    figure(9)
    data_freqs=unique(freq_data(:,5));
    n_data_freqs=length(data_freqs);
    for f=1:n_data_freqs
        ids=find(freq_data(:,5)==data_freqs(f) & freq_data(:,7)~=0);
        filtered_freq=freq_data(ids,:);
        [a,b]=sort(filtered_freq(:,2));
        plot3(filtered_freq(b,5), filtered_freq(b,2), filtered_freq(b,3)./100, 'r','linewidth',2)
        hold on
    end
    
    hold on
    
    [freq_grid, dt_grid] = meshgrid(1:freq_step:freq_max, dtmin:step_dt:dtmax);
    STDP_interpol = griddata(freq_htmp(:,1), freq_htmp(:,2), freq_htmp(:,3), freq_grid, dt_grid);
    ribbon(dt_grid,STDP_interpol)
    % surf(freq_grid, dt_grid, STDP_interpol);
    alpha 0.3