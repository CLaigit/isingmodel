function [ E ] = energy( net )
%ENERGY Summary of this function goes here
%   Detailed explanation goes here
    [column, row] = size(net);
    Ecolumn = zeros(column, 1);
    E = 0;
    for i = 1:column - 1
        Ecolumn = Ecolumn + net(:,i).*net(:,i+1);
    end
    Ecolumn = Ecolumn + net(:,end).*net(:,1);
    Ecolumn = 2*Ecolumn;
    
    for j = 1:column-1
        E = E + Ecolumn(i,1) * Ecolumn(i+1,1);
    end
    E = E + Ecolumn(end,1) * Ecolumn(1,1);
    E = 2*E;
    
    E = -E / column/row;

end

