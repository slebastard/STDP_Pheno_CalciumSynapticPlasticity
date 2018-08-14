function simu = default_simu()

    simu.model = 'pheno';
    simu.d_t = 30;
    simu.n_iter = 100;
    simu.frequency = 1;
    simu.int_scheme = 'euler_expl';
    simu.scheme_step = 0.5;
    simu.T = 1000;

end