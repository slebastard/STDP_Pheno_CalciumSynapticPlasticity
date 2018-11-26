function stampNetworkParams( fig, names, values, pos)
%STAMPPARAMS Summary of this function goes here
% fig is a handle to the figure to stamp parameters to
% names contains the names of the parameters as an array
% values contains the values in the same order
% pos is an object: 'x' is the start x value
%                   'y' is the start y value
%                   'ftSize'
%                   'colSep'
%                   'rowSep'

set(0, 'CurrentFigure', fig)
l = length(names);

for paramID=1:l
    text(pos.x + pos.colSep*floor(paramID/5), pos.y + (pos.rowSep-1)*rem(paramID,5), char(strcat(names(paramID) + " = " + num2str(values(paramID)))), 'FontSize', pos.ftSize);
end

end

