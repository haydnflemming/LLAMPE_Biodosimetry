% -- SNR FUNCTION ---------------------------------------------------------
% Written by Cristian Ciobanu
%
% SNR calculation. Assumes data is background subtracted. (I recommend
% removing outliers, cosmic rays, and background.) 'data' contains the
% Raman spectra. 'peaks' is a row vector indicating the column numbers of
% the peaks you wish to calculate SNR for. the 'SNR' output is a row vector
% with the same size as 'peaks'; each entry corresponds to the SNR at the
% respective peak.

function SNR = SNR_calc(data, peaks)
    data_std = std(data);
    data_mean = mean(data);
    
    SNR = data_mean(peaks) ./ data_std(peaks);
end