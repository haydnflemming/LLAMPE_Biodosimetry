% -- NIST CORRECTOR -------------------------------------------------------
% Written by Cristian Ciobanu
% edited 02/08/23
% 
% Applies NIST correction to data. 'bloodSpectra' contains Raman spectra.
% 'correctionData' contains all NIST spectra. 'shifts' is a row vector
% corresponding to Raman shifts in 'bloodSpectra' and 'correctionData'.

function nistCorrected = nistCorrect(bloodSpectra, correctionData, shifts)
    cosRayCorrected = CRremove(correctionData,6,2,1,1);
    % minimum = min(cosRayCorrected, [], 2);
    % maximum = max(cosRayCorrected, [], 2);

    normalized = mean(cosRayCorrected);
    % normalized = (cosRayCorrected - minimum) ./ (maximum - minimum);
    % normalized = mean(normalized);
    
    a0 = 9.71937E-2;
    a1 = 2.28325E-4;
    a2 = -5.86762E-8;
    a3 = 2.16023E-10;
    a4 = -9.77171E-14;
    a5 = 1.15596E-17;
    
    idealCurve = a0 + (a1 .* shifts) + (a2 .* shifts.^2) + (a3 .* shifts.^3) ...
        + (a4 .* shifts.^4) + (a5 .* shifts.^5);
    
    correctionCurve = idealCurve ./ normalized;
    nistCorrected = bloodSpectra .* correctionCurve ;
end