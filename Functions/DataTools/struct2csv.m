function struct2csv(s,fn)
% STRUCT2CSV(s,fn)
%
% Output a structure to a comma delimited file with column headers
%
%       s : any structure composed of one or more matrices and cell arrays
%      fn : file name
%
%      Given s:
%
%          s.Alpha = { 'First', 'Second';
%                      'Third', 'Fourth'};
%
%          s.Beta  = [[      1,       2;
%                            3,       4]];
%          
%          s.Gamma = {       1,       2;
%                            3,       4};
%
%          s.Epsln = [     abc;
%                          def;
%                          ghi];
% 
%      STRUCT2CSV(s,'any.csv') will produce a file 'any.csv' containing:
%
%         "Alpha",        , "Beta",   ,"Gamma",   , "Epsln",
%         "First","Second",      1,  2,      1,  2,   "abc",
%         "Third","Fourth",      3,  4,      3,  4,   "def",
%                ,        ,       ,   ,       ,   ,   "ghi",
%    
%      v.0.9 - Rewrote most of the code, now accommodates a wider variety
%              of structure children
%
% Written by James Slegers, james.slegers_at_gmail.com
% Covered by the BSD License
%

if exist(fn, 'file') == 2
    head = 0;
else
    head = 1;
end

FID = fopen(fn,'a');

headers = fieldnames(s);
m = numel(headers);


if head == 1
    l = '';
    if m>1
        for ii = 1:m-1  
            l = [l,'"',headers{ii},'",'];
        end
    end
    l = [l,'"',headers{m},'"'];
    l = [l,'\n'];
    fprintf(FID,l);
end

l = '';
if m>1
    for j = 1:m-1
        c = s.(headers{j});
        str = '';

        if isnumeric(c)
            str = [str,num2str(c),','];
        elseif islogical(c)
            str = [str,num2str(double(c)),','];
        elseif ischar(c)
            str = ['"',c,'",'];
        elseif iscell(c)
            if isnumeric(c{1,1})
                str = [str,num2str(c),','];
            elseif islogical(c{1,1})
                str = [str,num2str(double(c)),','];
            elseif ischar(c{1,1})
                str = [str,'"',c,'",'];
            end
        end
        l = [l,str];
    end
end

c = s.(headers{m});
str = '';

if isnumeric(c)
    str = [str,num2str(c)];
elseif islogical(c)
    str = [str,num2str(double(c))];
elseif ischar(c)
    str = ['"',c,'"'];
elseif iscell(c)
    if isnumeric(c{1,1})
        str = [str,num2str(c)];
    elseif islogical(c{1,1})
        str = [str,num2str(double(c))];
    elseif ischar(c{1,1})
        str = [str,'"',c];
    end
end
l = [l,str];

l = [l,'\n'];
fprintf(FID,l);
fclose(FID);
