function Langevin_GB2007_Stabilization(n_iter)
    
    t0=0; tfinal=100; dt=0.05;
    Sp_a = []; Sp_infty = [];
    for c=0.25:0.05:5
        sa = []; si=[];
        for iter=1:n_iter
            [a, b] = Langevin_GB2007(t0, tfinal, dt, c);
            sa = [sa; a];
            si = [si; b];
        end
        mu_alpha = mean(sa); sig_alpha = std(sa); Sp_a = [Sp_a; sa];
        mu_infty = mean(si); sig_infty = std(si); Sp_infty = [Sp_infty; si];
    end
    
    figure()
    plot(mu_alpha, mu_infty);
    title('Time-limit mean as a function of mean at threshold')
    
    figure()
    plot(mu_alpha,sig_alpha);
    title('Std deviation as a function of mean at threshold')
    
    figure()
    plot(mu_infty, sig_infty);
    title('Time-limit std deviation as a function of time-limit mean')
end