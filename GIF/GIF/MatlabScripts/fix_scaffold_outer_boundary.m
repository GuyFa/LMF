function [ x, FixedIndices, FixedValues ] = fix_scaffold_outer_boundary( initialValue, outer_border_size )
%fix the first cone - translation degree of freedom
    n=length(initialValue)/2; 
    FixedIndices = [(1:outer_border_size).'; (n+1:n+outer_border_size).']; 
    FixedValues = [initialValue(1:outer_border_size);initialValue(n+1:n+outer_border_size)]; 
    x=initialValue(setdiff(1:end,FixedIndices));
end

