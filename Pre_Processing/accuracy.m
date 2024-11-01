% -- ACCURACY CALCULATOR --------------------------------------------------
% Written by Cristian Ciobanu
%
% Compares true ('real') class values to predicted ('pred') class values to
% report accuracy of predictions. Inputs are row vectors.

function acc = accuracy(real, pred)
    acc = real == pred;
    acc = sum(acc)/size(acc,1);
end