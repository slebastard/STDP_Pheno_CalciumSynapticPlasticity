function env = getEnv(  )
%GETENV Summary of this function goes here
%   Detailed explanation goes here
    env.functionsRoot = 'Functions/';
    env.dataRoot = '../Data/';
    env.outputRoot = '../Output/';
    
    set(groot,'defaultfigureposition',[400 250 900 750])
    set(0,'DefaultFigureWindowStyle','docked')
end

