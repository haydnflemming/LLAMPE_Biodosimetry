% -- BASELINE CORRECTION --------------------------------------------------
% Written by Cristian Ciobanu
%
% Performs a baseline correction (corrects for hot pixels). 'data' contains
% all the Raman spectra to be corrected. 'baselineData1' and
% 'baselineData2' are the mean baseline measurements before and after Raman
% spectrum collection respectively (though order does not really matter).
% They are both row vectors.

function baselineCorrected = baselineRemove(data, baselineData1, baselineData2)
    baseline = mean([mean(baselineData1); mean(baselineData2)]);
    baselineCorrected = data - baseline;
end