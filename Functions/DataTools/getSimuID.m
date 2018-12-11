function ID = getSimuID(fn)
%GETSIMUID Summary of this function goes here
%   Detailed explanation goes here
ID = 0; stop = 0; c = 0;

if exist(fn, 'file') ~= 2
    ID = 1;
    return
end
f = fopen(fn,'r');

while ~stop
    L = fgetl(f);
    if ~ischar(L)
        stop = 1;
    else
        spl = strsplit(L,',');
        if c > 0
            ID = str2num(spl{1});
        end
    end
     c = c + 1;
end
ID = ID + 1;

end

