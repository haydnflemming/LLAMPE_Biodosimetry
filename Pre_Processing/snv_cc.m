% -- STANDARD NORMAL VARIATE ----------------------------------------------
% Written by Cristian Ciobanu
%
% Performs standard normal variate (SNV) correction on data.

function snvCorrected = snv(data)
    meanValues = mean(data, 2);
    stdDev = zeros(size(data, 1), 1);
    
    for i = 1:size(data, 1)
        stdDev(i, 1) = std(data(i, :));
    end
    
    snvCorrected = (data - meanValues) ./ stdDev;
end