function t = poissonMap( pSTDP )
%CORRPOISSON Generates uncorrelated Poisson processes
    
    pSTDP.T;
    pSTDP.nuPre.min;
    pSTDP.nuPre.max;
    pSTDP.nuPre.step;
    pSTDP.nuPost.min;
    pSTDP.nuPost.max;
    pSTDP.nuPost.step;
    rho_max = pSTDP.rho.max;
    rho0_step = pSTDP.rho.step;
    
    ?? prot
    ?? params
    ?? STDP
    ?? no model distinction
    
    map = zeros(1+floor((pSTDP.nuPre.max - pSTDP.nuPre.min)/pSTDP.nuPre.step), 1+floor((pSTDP.nuPost.max - pSTDP.nuPost.min)/pSTDP.nuPost.step), 1+floor(rho_max/rho0_step));
    
    for nu_pre=linspace(pSTDP.nuPre.min, pSTDP.nuPre.max, pSTDP.nuPre.step)
        for nu_post=linspace(pSTDP.nuPost.min, pSTDP.nuPost.max, pSTDP.nuPost.step)
            pre_id = ;
            post_id = ;
            
            t = indPoisson( 2, [nu_pre; nu_post], pSTDP.T);
            pre_spikes_hist = t(1,:);
            post_spikes_hist = t(2,:);
            for rho_0 = 0:rho0_step:rho_max
                rho0_id = 1 + floor(rho_0/rho0_step);
                w_0 = transfer(rho_0, prot);
                params.rho_0 = rho_0;
                params.w_0 = w_0;
                [~, w_end, ~] = pheno_model_efficient(pre_spikes_hist, post_spikes_hist, params, STDP);
                q_w = w_end/w_0;
                map(pre_id, post_id, rho0_id) = q_w;
            end
            
        end
    end
    
end

