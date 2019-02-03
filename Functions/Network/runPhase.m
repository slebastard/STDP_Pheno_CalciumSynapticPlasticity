function [net, plt, splSyn, tuning] = runPhase(net, neu, syn, input, prot, phase, plt, splSyn, gif, nIterTot)
% runPhase Summary of this function goes here
%   Detailed explanation goes here

bSize = floor(input.T / phase.dt);
stims = find(phase.stim);
nStims = length(stims);

phase.nIter = phase.T / phase.dt;

ExInput = zeros(net.NIn,bSize);
for i=phase.firstIter:phase.lastIter
    if (1+mod(i-1,bSize))==1 && ~(input.structured)
%         inSpikeTimes = indPoisson( net.NIn, net.C_ext*net.nu_ext, bSize*phase.dt );
%         inSpikeIter = ceil(inSpikeTimes./phase.dt);
%         inSpikeIter = inSpikeIter(:,2:end);
%         inSpikeIter = (inSpikeIter<bSize).*inSpikeIter  + (inSpikeIter>=bSize).*bSize;
%         ExInput = zeros(net.NIn,bSize);
%         
%         for k=1:bSize
%             ExInput((k-1)*net.NIn*ones(1,net.NIn) + (1:net.NIn)) = transpose(sum(inSpikeIter == k, 2));
%         end
%         
%         ExInput(net.NIn*(bSize-1):net.NIn*bSize - 1) = 0;
          ExInput = transpose(indPoisson2(net.C_ext*net.nu_ext, bSize, net.NIn, phase.dt));
    end
    
    if (1+mod(i-1,bSize))==1 && input.structured
        % Make sure bSize is greater than the input signal duration
        nReps = ceil(bSize*phase.dt/input.T);
        ExInput = zeros(net.NIn,bSize);
        for j=1:net.NIn
            inSpikeIter = [];
            for stID=1:nStims
                inSpikeIter = [inSpikeIter; max(1,floor(input.spikes{stims(stID)}{j}/phase.dt))];
            end
            for k=2:nReps
                inSpikeIter = [inSpikeIter, max(1,floor((k-1)*floor(input.T/phase.dt) + inSpikeIter))];
            end
            inSpikeIter = inSpikeIter(inSpikeIter<=bSize);
            ExInput(j,inSpikeIter) = 1;
        end
        ExInput = ExInput + transpose(indPoisson2(net.nu_inNoise, bSize, net.NIn, phase.dt));
    end

    % Who (input) is spiking?
    spikeIn = ExInput(:,1+mod(i-1,bSize));
    spikeIDIn = find(spikeIn);
    
    % Input passive evolution
    net.caIn = net.caIn.*exp(-phase.dt/syn.tau_Ca);
    net.xpreIn = 1.*net.xpreIn;% - exp(-phase.dt/syn.tau_pre).*(1-net.xpreIn);
    net.xpostIn = 1.*net.xpostIn;% - exp(-phase.dt/syn.tau_post).*(1-net.xpostIn);        
    
    % Recurrent passive evolution
    tIn = net.WIn*spikeIn;
    tRec = net.RI(:,1+mod(i-1,neu.N_del));
    net.V = (1-phase.dt/neu.tau)*net.V + net.WIn*spikeIn + net.RI(:,1+mod(i-1,neu.N_del));
    net.V(net.LS>i-neu.N_rp)=neu.V_r;
    plt.Vspl(i) = net.V(40);
    net.ca = net.ca.*exp(-phase.dt/syn.tau_Ca);
    net.xpre = 1.*net.xpre;% - exp(-phase.dt/syn.tau_pre).*(1-net.xpre);
    net.xpost = 1.*net.xpost;% - exp(-phase.dt/syn.tau_post).*(1-net.xpost);        

    % Who (recurrent) is spiking?
    spikeRec = (net.V>=neu.V_t);
    spikeIDRec = find(spikeRec);
    
    % Input synapses update
    net.caIn(spikeIDRec,:) = net.caIn(spikeIDRec,:) + syn.C_post.*net.xpostIn(spikeIDRec,:);
    net.caIn(:,spikeIDIn) = net.caIn(:,spikeIDIn) + syn.C_pre.*net.xpreIn(:,spikeIDIn);
    net.actPotIn = (net.caIn >= syn.theta_pot).*(net.WIn~=0);
    net.actDepIn = (net.caIn >= syn.theta_dep).*(net.WIn~=0);
    net.xpostIn(spikeIDRec,:) = net.xpostIn(spikeIDRec,:)*(1 - syn.dampFactor);
    net.xpreIn(:,spikeIDIn) = net.xpreIn(:,spikeIDIn)*(1 - syn.dampFactor);
    net.rhoIn = net.rhoIn + phase.PlastInON .* phase.dt./syn.tau_rho .* (syn.gamma_pot.*(syn.rho_max - net.rhoIn).*net.actPotIn - syn.gamma_dep.*net.rhoIn.*net.actDepIn);
    %net.WIn = 2.*syn.J.*transfer(net.rhoIn, prot)./net.NIn;
    net.WIn = net.synSignIn.*transfer(net.rhoIn, prot);
    
    % Recurrent synapses update
    net.ca(spikeIDRec,:) = net.ca(spikeIDRec,:) + syn.C_post.*net.xpost(spikeIDRec,:);
    net.ca(:,spikeIDRec) = net.ca(:,spikeIDRec) + syn.C_pre.*net.xpre(:,spikeIDRec);
    net.actPot = (net.ca >= syn.theta_pot).*(net.W~=0);
    net.actDep = (net.ca >= syn.theta_dep).*(net.W~=0);
    net.xpost(spikeIDRec,:) = net.xpost(spikeIDRec,:)*(1 - syn.dampFactor);
    net.xpre(:,spikeIDRec) = net.xpre(:,spikeIDRec)*(1 - syn.dampFactor);
    net.rho = net.rho + phase.PlastON .* phase.dt./syn.tau_rho .* (syn.gamma_pot.*(syn.rho_max - net.rho).*net.actPot - syn.gamma_dep.*net.rho.*net.actDep);

    % Creation of new synapses
    t = rand(net.N,net.N);
    new_exc = (t < phase.birth.exc*phase.dt/syn.tau_rho);
    new_inh = (t >= phase.birth.exc*phase.dt/syn.tau_rho & t < (phase.birth.exc + phase.birth.inh)*phase.dt/syn.tau_rho);
    net.synSign(new_exc) = 2*syn.J;
    net.synSign(new_inh) = -2*syn.J.*net.g;

    % Weight updates
    net.W = net.synSign.*transfer(net.rho, prot);
    
    % Recurrent reset
    net.RI(:,1+mod(i,neu.N_del))=net.W*spikeRec;
    net.V(spikeRec)=neu.V_r;
    
    % Stats & recordings
    net.meanWexc(i+1,1) = mean(net.W(net.W>0));
    net.meanWinh(i+1,1) = mean(net.W(net.W<0));
    
    plt.histRho(:,i) = (1/net.N^2).*histcounts(net.rho(net.rho>0),plt.edgesRho);
    plt.histW_exc(:,i) = (1/net.N^2).*histcounts((0.5/syn.J).*net.W((0.5/syn.J).*net.W>0),plt.edgesW_exc);
    plt.histW_inh(:,i) = (1/net.N^2).*histcounts((0.5/syn.J).*net.W((0.5/syn.J).*net.W<0),plt.edgesW_inh);
    net.LS(spikeRec)=i;                % Time of last spike
    
    plt.Rasterplot(:,i)=spikeRec;
    plt.RasterplotIn(:,i)=(spikeIn>0);
    
    % Coloring neurons currently in refractory
    % plt.Rasterplot(net.LS<i & net.LS>i-neu.N_del,i) = 0.1;
    
    %
    %plt.Rasterplot(net.WIn*spikeIn~=0,i) = 0.1;
    
    % Updating synapse sample history
    splSyn.ca(:,i) =  2*(net.ca(splSyn.IDs)>=syn.theta_pot) - (net.ca(splSyn.IDs)>=syn.theta_dep);
    %splSyn.ca(:,i) = net.ca(splSyn.IDs);
    splSyn.xpre(:,i) =  net.xpre(splSyn.IDs);
    splSyn.xpost(:,i) =  net.xpost(splSyn.IDs);
    splSyn.rho(:,i) =  net.rho(splSyn.IDs);
    splSyn.w(:,i) = net.W(splSyn.IDs);
    
%     if mod(i,gif.updIter)==0
%         % PrintiIDRec,:) + syn.C_post.*net.xpostIn(spikeIDRec,:);
%         
%         if gif.graph
%             G = digraph(net.W([subIDsExc;subIDsInh],[subIDsExc;subIDsInh])');
%             LWiphase.dths = (5.*(G.Edges.Weight>0) + 0.8.*(G.Edges.Weight<0)).*abs(G.Edges.Weight);
%             EColors = [0 1 0].*(G.Edges.Weight>0) + [1 0 0].*(G.Edges.Weight<0);
%             figure(1)
%             plot(G,'LineWidth',LWiphase.dths,'EdgeColor',EColors)
%             f1 = getframe;
%             im1(:,:,1,floor(i/gif.updIter)) = rgb2ind(f1.cdata,map1,'nodither');
%         end
%         
%         % Eigendecomposition to find clusters
%         if gif.lapl
%             A = abs(net.W + net.W');
%             D = diag(sum(A,1));
%             L_rw = eye(net.N) - D^(-1)*A;
%             [eVals, ~] = eig(L_rw);
%             figure(2)
%             spc = sort(diag(eVals));
%             plot(1:20, spc(1:20), '+r')
%             xlabel('Index')
%             ylabel('Eigenvalue')
%             title(strcat('Eigendecomposition at time ', num2str(phase.dt*i,3),'s'))
%             f2 = getframe;
%             im2(:,:,1,floor(i/gif.updIter)) = rgb2ind(f2.cdata,map2,'nodither');
%         end
%     end
   
    if plt.progress
        progressbar(i/nIterTot)
    end
end

if phase.specAnalysis
    figName = strcat('Specificity analysis - Phase ', num2str(phase.ID));
    figure('Name',figName,'NumberTitle','off')
    
    % Compute firing rate
    filter = gausswin(100,3);
    filter = (1/sum(filter)).*filter;
    firing.rate = zeros(net.N, phase.lastIter - phase.firstIter + 1);
    for n=1:net.N
       firing.rate(n,:) = (1/phase.dt).*conv(plt.Rasterplot(n,phase.firstIter:phase.lastIter), filter, 'same');
    end
    
    % Superimpose firing rate and stimulus value
    tuning.nBins = 500;
    tuning.input = input.stim{stims(1)};
    tuning.curve = zeros(tuning.nBins, net.N);
    
    for b=1:tuning.nBins
        inst = (tuning.input >= (b-1)*pi/tuning.nBins & tuning.input < b*pi/tuning.nBins);
        for n=1:net.N
            tuning.curve(b,n) = mean(firing.rate(n,inst));
        end
    end
    
    tuning.curve(isnan(tuning.curve)) = 0;
    plot(tuning.curve(:,1));
    plot(tuning.curve(:,9));
    
end
end

