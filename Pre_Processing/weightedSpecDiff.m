% -- WEIGHTED SPECTRAL DIFFERENCE -----------------------------------------
% Written by Cristian Ciobanu
%
% This function calculates the weighted spectral difference as defined by
% Dinh et al. (2014) between two spectra (a reference and sample).

function wsd = weightedSpecDiff(ref, sample)
    wsd = sqrt( sum( ( abs(ref) / abs(mean(ref)) ) .* ( ( ref-sample ).^2 ) ) / size(sample,2) );
end