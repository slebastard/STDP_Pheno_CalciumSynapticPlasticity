function STDP = poissonCut( params, simu )
%CORRPOISSON Generates uncorrelated Poisson processes
% All rates in Hz
    corr = simu.corr;
    
    T = simu.T;
    dir = simu.dir;
    val = simu.val;
    nuMin = simu.nu.min;
    nuMax = simu.nu.max;
    nuStep = simu.nu.step;
    nTry = simu.nTry;
    rho_max = params.rho_max;
    rho0_step = 5;
    
    nNu = 1 + floor((nuMax-nuMin)/nuStep);
    nRho = 1+floor(rho_max/rho0_step);
    map_cache = zeros(nNu, nRho, nTry);
    map_mean = zeros(nNu, nRho);
    map_corr = zeros(nNu, nRho, nRho);
    
    for nu=nuMin:nuStep:nuMax
        id = 1 + floor((nu-nuMin)/nuStep);
        for tryID = 1:nTry
            if strcmp(dir,'pre')
                if strcmp(corr.type, 'none')
                    t = indPoisson( 2, [1000/val; 1000/nu], simu.T);
                    pre_spikes_hist = t(1,:);
                    post_spikes_hist = t(2,:);
                elseif strcmp(corr.type, 'exp')
                    c12 = corr.c12;
                    tc = corr.tc;
                    C = [c12*1.05*nu_pre/nu_post c12; c12 c12*1.05*nu_post/nu_pre];
                    [t, I] = corrPoisson( 2, [val; nu], C, simu.T, tc);
                    preIds = find(I(1,:));
                    pre_spikes_hist = t(preIds);
                    postIds = find(I(2,:));
                    post_spikes_hist = t(postIds);
                end
            else
                if strcmp(corr.type, 'none')
                    t = indPoisson( 2, [1000/nu; 1000/val], simu.T);
                    pre_spikes_hist = t(1,:);
                    post_spikes_hist = t(2,:);
                elseif strcmp(corr.type, 'exp')
                    c12 = corr.c12;
                    tc = corr.tc;
                    C = [c12*1.05*nu_pre/nu_post c12; c12 c12*1.05*nu_post/nu_pre];
                    [t, I] = corrPoisson( 2, [nu; val], C, simu.T, tc);
                    preIds = find(I(1,:));
                    pre_spikes_hist = t(preIds);
                    postIds = find(I(2,:));
                    post_spikes_hist = t(postIds);
                end
            end
            for rho_0 = 0:rho0_step:rho_max
                rho0_id = 1 + floor(rho_0/rho0_step);
                w_0 = transfer(rho_0, params);
                params.rho_0 = rho_0;
                params.w_0 = w_0;
                
                if strcmp(simu.model, 'caProd')
                    [~, w_end, ~] = caProd_model_efficient(pre_spikes_hist, post_spikes_hist, params, simu);
                else
                    [~, w_end, ~] = pheno_model_efficient(pre_spikes_hist, post_spikes_hist, params, simu);
                end
                map_cache(id, rho0_id, tryID) = w_end/w_0;

                progressbar((rho0_id + (rho_max+1)*(id-1))/(nRho + (rho_max+1)*(nNu-1)));
            end
        end
    end

    map_mean = mean(map_cache, 3);   
    tmp = transfer_ind(0:rho0_step:rho_max, params);
    muW = 0.5;
    rhoDistr = (1/(sqrt(2*pi*muW*(1-muW)))) .* exp(-(tmp - muW*ones(1, 1+floor(rho_max/rho0_step))).^2 ./(2*muW*(1-muW))); 
    rhoDistr = (1/sum(rhoDistr)).*rhoDistr;
    tMean = map_mean*rhoDistr';
    
    var = zeros(nNu,1);
    for nuID=1:nNu
        A = reshape(map_cache(nuID,:,:), [nRho nTry]);
        map_corr(nuID,:,:) = (1/nTry) .* A * A' - map_mean(nuID,:)' * map_mean(nuID,:);
        
        B = reshape(map_corr(nuID,:,:), [nRho nRho]);
        var(nuID,1) = rhoDistr * B * rhoDistr';
    end
    
    STDP.mean = cat(2, (nuMin:nuStep:nuMax)', tMean);
    STDP.var = cat(2, (nuMin:nuStep:nuMax)', var);
end

