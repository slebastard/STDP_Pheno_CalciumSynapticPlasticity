figure;
if strcmp(mode, 'dataFit')
    dtmin = -100;
    dtmax = 100;
    step_dt = 1;
    dt_params=[dtmin, dtmax, step_dt];
    
    freq_max = 10;
    freq_step = 0.5;
    freq_params = [freq_max, freq_step];
    
    heatmap_params = [model_params, n_iter];
    freq_htmp = get_freq_heatmap(model, 'rel', heatmap_params, int_scheme, dt_params, freq_params);
    
    AAA=unique(freq_data(:,5));
    nf=length(AAA);
    for u=1:nf
        coeffs=find(freq_data(:,5)==AAA(u));
        dats=freq_data(coeffs,:);
        [a,b]=sort(dats(:,2));
        plot3(dats(b,5), dats(b,2), dats(b,3)./100, 'r','linewidth',2)
        hold on
    end
    
    [freq_grid, dt_grid] = meshgrid(.1:freq_step:freq_max, dtmin:step_dt:dtmax);
    STDP_interpol = griddata(freq_htmp(:,1), freq_htmp(:,2), freq_htmp(:,3), freq_grid, dt_grid);
    surf(freq_grid, dt_grid, STDP_interpol);
    alpha 0.3
    surf(freq_grid, dt_grid, ones(size(STDP_interpol)));
    alpha 0.3
end
%%
close all

dtmin = -100;
    dtmax = 100;
    step_dt = 10;
figure;
hold on;
surf(freq_grid, dt_grid, log(STDP_interpol));
alpha 0.3
surf(freq_grid, dt_grid, zeros(size(STDP_interpol)));
alpha 0.3
shading interp

for u=1:nf
        figure;
        coeffs=find(freq_data(:,5)==AAA(u));
        dats=freq_data(coeffs,:);
        dats=dats(abs(dats(:,3)-100)>20,:);
        [a,b]=sort(dats(:,2));
        plot(dats(b,2), (dats(b,3)./100), 'r','linewidth',2)
        hold on
        stdp_params = [model_params, dtmin, dtmax, step_dt, 100, AAA(u)];
        STDP = get_STDP(model, 'rel', stdp_params, int_scheme, scheme_step);
        plot(STDP(:,1), STDP(:,2), '.b')
        ylim([0.5,3])
end
        
%         data1Hz = freq_data(floor(freq_data(:,5))==1,:);
%     
%     stdp_params = [model_params, dtmin, dtmax, step_dt, 100, 1];
%     
%     plot(STDP(:,1), STDP(:,2), '.b')
%     hold on
%     plot(data1Hz(:,2),data1Hz(:,3)./100,'xr')
%     neutral_hline = refline([0 1]);
%     neutral_hline.Color = 'b';   
%     ('Plasticity as a function of pre-post spike delay');
%     xlabel('Pre-post spike delay (ms)');
%     ylabel('Relative change in synaptic strength');
% end
   