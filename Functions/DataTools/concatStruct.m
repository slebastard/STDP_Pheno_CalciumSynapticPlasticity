function A = concatStruct( A, B )
%CONCATSTRUCTURE Summary of this function goes here
%   Detailed explanation goes here
fields = fieldnames(B);
for i = 1:length(fields)
   A.(fields{i}) = B.(fields{i});
end

end

