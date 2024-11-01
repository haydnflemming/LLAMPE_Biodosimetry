% -- HOTELLING T2 VS Q RESIDUAL TEST --------------------------------------
% Written by Cristian Ciobanu
%
% This function calculates the Hotelling T2 and Q residual statistics of
% the input data (from PCA). Spectra with either T2 or Q values above the
% calculated 95% confidence limits are removed.
%
% NOTE: REQUIRES PLS_TOOLBOX!

function [outlierRemoved, indices] = hqTest(data)
    options = pca('options');
    options.display = 'off';
    options.outputversion = 2;
    options.preprocessing = {preprocess('default', 'mean center')};
    options.plots = 'none';
    options.confidencelimit = 0.95;
    [~, ~, ~, res, reslm, tsqlm, tsq] = pca(data, 3, options);
    
    outlierRemoved = data;
    outlierCheck = tsq > tsqlm | res > reslm;
    
    i = 1;
    k = 1;
    indices = [];
    
    while i <= size(outlierCheck, 1)
        if outlierCheck(i) == 1
            outlierRemoved(i,:) = [];
            outlierCheck(i,:) = [];
            i = i - 1;
            indices = [indices k];
        end
        
        i = i + 1;
        k = k + 1;
    end
end