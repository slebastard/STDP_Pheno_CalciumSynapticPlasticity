function stampParams( params, pos )
%STAMPPARAMS Summary of this function goes here
% params is an object that contains the parameters to print
% pos is an object contating x, y, z, colSepY, colSepZ and ftSize
%   Detailed explanation goes here
    text(pos.x, pos.y + 0.5*pos.colSepY, pos.z + 2*pos.colSepZ, strcat('\tau_{ca} = ', num2str(params.tau_Ca, 3), 'ms'), 'FontSize', pos.ftSize);
    text(pos.x, pos.y + 0.5*pos.colSepY, pos.z + pos.colSepZ, strcat('\tau_{damp} = ', num2str(params.tau_x, 3), 'ms'), 'FontSize', pos.ftSize);
    text(pos.x, pos.y + 0.5*pos.colSepY, pos.z, strcat('C_{pre} = ', num2str(params.C_pre, 3)), 'FontSize', pos.ftSize);
    text(pos.x, pos.y + 0.5*pos.colSepY, pos.z - pos.colSepZ, strcat('C_{post} = ', num2str(params.C_post, 3)), 'FontSize', pos.ftSize);

    text(pos.x, pos.y + 0.5*pos.colSepY, pos.z - 2*pos.colSepZ, strcat('\theta_{dep} = ', num2str(params.theta_dep, 3)), 'FontSize', pos.ftSize);
    text(pos.x, pos.y - 0.5*pos.colSepY, pos.z + 2*pos.colSepZ, strcat('\theta_{pot} = ', num2str(params.theta_pot, 3)), 'FontSize', pos.ftSize);

    text(pos.x, pos.y - 0.5*pos.colSepY, pos.z + pos.colSepZ, strcat('\gamma_{dep} = ', num2str(params.gamma_dep, 3)), 'FontSize', pos.ftSize);
    text(pos.x, pos.y - 0.5*pos.colSepY, pos.z, strcat('\gamma_{pot} = ', num2str(params.gamma_pot, 3)), 'FontSize', pos.ftSize);

    text(pos.x, pos.y - 0.5*pos.colSepY, pos.z - pos.colSepZ, strcat('\sigma = ', num2str(params.noise_lvl, 3)), 'FontSize', pos.ftSize); % 12 factor for effective noise correction - 1/sqrt(N_A*V);
    text(pos.x, pos.y - 0.5*pos.colSepY, pos.z - 2*pos.colSepZ, strcat('dampFactor = ', num2str(params.dampFactor, 3)), 'FontSize', pos.ftSize);
end

