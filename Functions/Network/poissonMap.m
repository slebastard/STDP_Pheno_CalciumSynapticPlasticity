function STDP = poissonMap( params, simu )
%CORRPOISSON Generates uncorrelated Poisson processes
% All rates in Hz
    corr = simu.corr;

    T = simu.T;
    nPreMin = simu.nuPre.min;
    nPreMax = simu.nuPre.max;
    nPreStep = simu.nuPre.step;
    nPostMin = simu.nuPost.min;
    nPostMax = simu.nuPost.max;
    nPostStep = simu.nuPost.step;
    nTry = simu.nTry;
    rho_max = params.rho_max;
    rho0_step = 5;
    
    nPre = 1 + floor((nPreMax-nPreMin)/nPreStep);
    nPost = 1 + floor((nPostMax-nPostMin)/nPostStep);
    nRho = 1+floor(rho_max/rho0_step);
    map = zeros(nPre, nPost, nRho);
    
    for nu_pre=nPreMin:nPreStep:nPreMax
        for nu_post=nPostMin:nPostStep:nPostMax
            pre_id = 1 + floor((nu_pre-nPreMin)/nPreStep);
            post_id = 1 + floor((nu_post-nPostMin)/nPostStep);
            cache_qw = zeros(1,nRho);
            for tryID = 1:nTry
                if strcmp(corr.type, 'none')
                    t = indPoisson( 2, [1000/nu_pre; 1000/nu_post], simu.T);
                    pre_spikes_hist = t(1,:);
                    post_spikes_hist = t(2,:);
                elseif strcmp(corr.type, 'exp')
                    c12 = corr.c12;
                    tc = corr.tc;
                    C = [c12*1.05*nu_pre/nu_post c12; c12 c12*1.05*nu_post/nu_pre];
                    [t, I] = corrPoisson( 2, [nu_pre; nu_post], C, simu.T, tc);
                    preIds = find(I(1,:));
                    pre_spikes_hist = t(preIds);
                    postIds = find(I(2,:));
                    post_spikes_hist = t(postIds);
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
                    cache_qw(rho0_id) = cache_qw(rho0_id) + w_end/w_0;
                    progressbar((rho0_id + (rho_max+1)*(post_id-1) + (rho_max*nPost+nPost)*(pre_id-1))/(1 + floor(rho_max/rho0_step) + (rho_max+1)*(nPost-1) + (rho_max*nPost+nPost)*(nPre-1)));

                end
                map(pre_id, post_id, :) = (1/nTry).*cache_qw;
            end
        end
    end
    
    tmp = transfer_ind(0:rho0_step:rho_max, params);
    muW = 0.5;
    rhoDistr = (1/(sqrt(2*pi*muW*(1-muW)))) .* exp(-(tmp - muW*ones(1, 1+floor(rho_max/rho0_step))).^2 ./(2*muW*(1-muW))); 
    rhoDistr = (1/sum(rhoDistr)).*rhoDistr;
    rhoConv = permute(repmat(rhoDistr',1,nPre,nPost),[2 3 1]);
    STDP = sum(map .* rhoConv,3);

end

