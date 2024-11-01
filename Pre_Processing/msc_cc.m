% -- MULTIPLICATIVE SCATTER CORRECTION ------------------------------------
% Written by Cristian Ciobanu
%
% Performs multiplicative scatter correction (MSC) on data. Mean spectrum
% is taken as the reference.

function [mscCorrected, coeff] = msc_cc(data)
    meanSpectrum = mean(data);
    coeff = zeros(size(data, 1), 2);
    
    for i = 1:size(data, 1)
        coeff(i, :) = polyfit(meanSpectrum, data(i, :), 1);
    end
    
    mscCorrected = (data - coeff(:, 2)) ./ coeff(:, 1);
end