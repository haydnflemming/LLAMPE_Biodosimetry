% -- MEAN WEIGHTED SPECTRAL DIFFERENCE ------------------------------------
% Written by Cristian Ciobanu
%
% The weighted spectral difference (WSD) formula on its own only calculates
% the difference between 2 spectra. For a set of spectra, this function
% uses the mean spectrum as the reference, calculates the WSD of each
% spectrum in relation to the mean spectrum, then outputs the mean WSD of
% the set.

function meanWSD = wsdMean(samples)
    num = size(samples, 1);
    meanSample = mean(samples);
    meanWSD = 0;
    
    for i = 1:num
        meanWSD = meanWSD + weightedSpecDiff(meanSample, samples(i, :));
    end
    
    meanWSD = meanWSD / num;
end