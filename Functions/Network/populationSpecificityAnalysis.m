function [outputArg1,outputArg2] = populationSpecificityAnalysis(inputArg1,inputArg2)
%POPULATIONSPECIFICITYANALYSIS Summary of this function goes here
%   Detailed explanation goes here

    % Defining stimulation protocol
    phase1.T = 0.3;
    phase1.dt = 5e-4;
    phase1.strap = 0;
    phase1.InE = 0;
    phase1.InI = 0;
    phase1.EE = 0;
    phase1.EI = 0;
    phase1.IE = 0;
    phase1.II = 0;
 
    phase2.T = 0.7;
    phase2.dt = 5e-4;
    phase2.strap = 0;
    phase2.InE = 0;
    phase2.InI = 0;
    phase2.EE = 1;
    phase2.EI = 1;
    phase2.IE = 0;
    phase2.II = 0;    

    phase3.T = 0.3;
    phase3.dt = 5e-4;
    phase3.strap = 0;
    phase3.InE = 0;
    phase3.InI = 0;
    phase3.EE = 0;
    phase3.EI = 0;
    phase3.IE = 0;
    phase3.II = 0;  
    
    std.simu.phases = {phase1, phase2, phase3};
    std.simu.T = 1.3;
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
    
    std.net.g = 4;
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
    
    % Run simulation with plasticity & specific input units
    run.out = simuNetwork(std.syn, std.simu, std.net, std.neu, std.init, std.plt, std.gif);
    
    
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

