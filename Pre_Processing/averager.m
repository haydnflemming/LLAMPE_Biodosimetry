% -- AVERAGER -------------------------------------------------------------
% Written by Cristian Ciobanu
%
% This function averages every n spectra (sequentially) in the 'data' set,
% while discaring any leftovers

function averaged = averager(data, n)
    if n == 1
        averaged = data;
    else
        numRows = size(data,1);
        remainder = mod(numRows,n);
        data(numRows-remainder+1:end,:) = [];

        averaged = [];

        for i = 1:n:numRows-remainder
            temp = data(i:i+n-1,:);
            temp = mean(temp);
            averaged = [averaged; temp];
        end
    end
end