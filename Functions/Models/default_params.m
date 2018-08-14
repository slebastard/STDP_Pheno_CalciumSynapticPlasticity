function params = default_params()

    params.rho_0 = 35;
    params.w_0 = transfer(params.rho_0);
    params.rho_max = 200;
    params.C_pre = 1;
    params.C_post = 2;
    params.tau_Ca = 20;
    params.delay_pre = 0;
    params.theta_dep = 1;
    params.gamma_dep = 200;
    params.theta_pot = 1.3;
    params.gamma_pot = 321;
    params.tau_rho = 150000;
    params.theta_act = 0.6;
    params.noise_lvl = 25;
    params.tau_w = 1000;
    params.S_attr = 40;

end