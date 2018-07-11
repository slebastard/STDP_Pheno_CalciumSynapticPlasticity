% function [paramCourse,yf] = FullDynamics_ParamAnalysis(model, inputs, paramName, paramSettings)
% Travels back and forth in a given range for a parameter, and gets the
% stable states of the system for those values
% - model is the name of the model that we will look at
% - inputs just contains stimulation descriptors
% -- (1): initial calcium value
% - paramSettings is an array containing
% -- (1): the name of parameter to vary, as declared in the model used
% -- (2): min value taken by this parameter
% -- (3): max value taken by this parameter
% -- (4): step size

model = 'NMDA';
inputs = 0.1;
paramName = 'CaBas';
paramSettings = [0.05, 1.50, 0.05];

paramSetName = 'Graupner';
varInd = 19; % ID of variable to plot as fct of param; set to Sp
CaBas = strcmp(paramName,'CaBas'); 

nSteps = floor((paramSettings(2)-paramSettings(1))/paramSettings(3));
initVals = [];
yf = zeros(2*nSteps,1);

if strcmp(model,'NMDA')
    for step=1:nSteps
        paramVal = paramSettings(1)+(step-1)*paramSettings(3);
        inputs(1) = paramVal*CaBas + inputs(1)*(1-CaBas);
        [~, ~, y] = FullDynamics_NMDA(inputs(1), initVals, paramSetName, false, false, paramName, paramVal);
        yf(step,1) = y(end,varInd);
        initVals = y(end,1:32);
    end 
    
    for step=nSteps:-1:1
        paramVal = paramSettings(1)+(step-1)*paramSettings(3);
        inputs(1) = paramVal*CaBas + inputs(1)*(1-CaBas);
        [~, ~, y] = FullDynamics_NMDA(inputs(1), initVals, paramSetName, false, false, paramName, paramVal);
        yf(1+2*nSteps-step,1) = y(end,varInd);
        initVals = y(end,1:32);
    end 
else
    error('Unknown model name, or model not supported for param analysis')
end

%%
figure()
paramRange = linspace(paramSettings(1),paramSettings(2),nSteps);
paramCourse = cat(2, paramRange, fliplr(paramRange));
plot(paramCourse, yf, 'x')
title(sprintf('Stable state of %s model as a function of %s', model, paramName));

% end