function sensibilityAnalysis( out_path )
%GETENV Summary of this function goes here
%   Detailed explanation goes here
    
    % Defining stimulation protocol
    phase1.T = 0.5;
    phase1.dt = 5e-4;
    phase1.strap = 0;
    phase1.InE = 0;
    phase1.InI = 0;
    phase1.EE = 0;
    phase1.EI = 0;
    phase1.IE = 0;
    phase1.II = 0;
 
    phase2.T = 1;
    phase2.dt = 5e-4;
    phase2.strap = 0;
    phase2.InE = 0;
    phase2.InI = 0;
    phase2.EE = 1;
    phase2.EI = 1;
    phase2.IE = 0;
    phase2.II = 0;    

    phase3.T = 0.5;
    phase3.dt = 5e-4;
    phase3.strap = 0;
    phase3.InE = 0;
    phase3.InI = 0;
    phase3.EE = 0;
    phase3.EI = 0;
    phase3.IE = 0;
    phase3.II = 0;  
    
    std.simu.phases = {phase1, phase2, phase3};
    std.simu.T = 2;
    std.simu.nPhases = 3;
    
    % Defining default parameter values
    std.syn = getSynapse();
    std.syn.J = 0.2;
    
    std.net.N = 200;
    std.net.NIn = 200;
    std.net.NE = 160;
    std.net.Connectivity = 0.1;
    std.net.D = 3e-3;
    
    std.neu.V_r = 10;
    std.neu.V_t = 20;
    std.neu.t_rp = 2e-3;
    
    std.net.g = 80;
    std.net.rExtRel = 4;
    
    std.init.c = 1;
    std.init.mode = 'rand';

    std.plt.timeSpl.n = 3;
    std.plt.timeSpl.dur = 40;
    
    std.plt.all.raster = 0;
    std.plt.spl.pres = 0;
    std.plt.spl.hist = 0;
    std.plt.spl.ca = 0;
    std.plt.spl.rho = 0;
    std.plt.spl.w = 0;
    std.plt.spl.phase = 0;
    std.plt.progress = 0;
    
    std.gif.graph = 0;
    std.gif.N = 0;
    std.gif.lapl = 0;
    std.gif.K = 0;
    
    % Creating exploration sets
    
    syn.C_pre = [0.1, 0.2, 0.8, 0.84];
    syn.C_post = [0.5, 0.7, 1, 1.5];
    syn.theta_pot = [1.05, 1.15, 1.35, 1.7];
    syn.theta_dep = [0.4, 0.6, 0.8, 1];
    syn.gamma_pot = [20, 60, 100, 200];
    syn.gamma_dep = [60, 100, 250, 400];
    syn.S_attr = [20, 30, 60, 100];    
    syn.J = [2e-3, 2e-2, 4e-1, 5];
    
    net.rExtRel = [0.7, 1.5, 2.5, 4];
    
    net.N = [180,200,250,300];
    net.Connectivity = [0.05,0.1,0.3,0.7];
    net.NE = [100,140,160,190];
    net.D = [5e-4,1e-3,5e-3,1e-2];
    net.g = [0.5, 2, 4, 8];
    
    % Exploring sets
    fn_syn = fieldnames(syn);
    fn_net = fieldnames(net);
    
    for k=1:numel(fn_syn)
        for lvl=1:4
            run = std;
            fieldList = syn.(fn_syn{k});
            run.syn.(fn_syn{k}) = fieldList(lvl);
            run.out = simuNetwork(run.syn, run.simu, run.net, run.neu, run.init, run.plt, run.gif);
            
            run.out.ID = getSimuID(out_path);
            struct2csv(run.out, out_path);
        end
    end
    
    for k=1:numel(fn_net)
        for lvl=1:4
            run = std;
            fieldList = net.(fn_net{k});
            run.net.(fn_net{k}) = fieldList(lvl);
            run.out = simuNetwork(run.syn, run.simu, run.net, run.neu, run.init, run.plt, run.gif);
            
            run.out.ID = getSimuID(out_path);
            struct2csv(run.out, out_path);
        end        
    end
end

