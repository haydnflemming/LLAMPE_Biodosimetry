% -- RAMAN SHIFT CONVERTER ------------------------------------------------
% Written by Cristian Ciobanu
%
% Converts wavelengths to Raman shifts. 'specWavelengths' is a row vector
% containing wavelengths for the Raman shifts.
%
% NOTE: INPUTS ARE IN NANOMETERS!

function converted = shiftConvert(specWavelengths, laserWavelength) % inputs are in nm
    converted = (1e7 / laserWavelength) - (1e7 ./ specWavelengths);
end