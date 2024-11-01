% -- COVARIATE ADJUSTMENT -------------------------------------------------
% Written by Cristian Ciobanu
%
% Performs covariate adjustment on Raman spectra. 'data' contains Raman
% spectra with first three rows corresponding to donor, dose, and serial
% number. It is important to ensure the first column contains donor number.
% 'covariates' is a matrix with the first column containing donor number,
% and respective rows containing donor information (e.g., age, gender, WBC
% count, RBC count, HGB count, etc.). If a covariate is categorical (e.g.,
% donor either IS or ISN'T a smoker), then assign a 0 and 1 to the two
% outcomes (e.g., IS smoker -> 1, ISN'T smoker -> 0). The function assigns
% donor info from 'covariates' to spectra in 'data' that contain matching
% donor numbers in the first column. It then performs a linear regression
% along all variables/Raman shifts, and outputs the residuals. These
% residuals are the covariate-adjusted data.

function res = covAdj(data, covariates)
    sizeData = size(data);
    sizeCov = size(covariates);
    
    if isa(covariates,'table')
        covariates = covariates{:,:};
    end

    if sizeData(1) ~= sizeCov(1)
        cov = zeros(sizeData(1),sizeCov(2)-1);
        for i = 1:sizeData(1)
            donor = data(i,1);
            for j = 1:sizeCov(1)
                if covariates(j,1) == donor
                    cov(i,:) = covariates(j,2:end);
                    break;
                end
            end
        end
    else
        cov = covariates;
    end
    
    temp = data;
    data = data(:,4:end);
    sizeData = size(data);
    res = zeros(sizeData);
    
    for i = 1:sizeData(2)
        [~,~,r] = regress(data(:,i),[ones(size(cov,1),1) cov]); % column of ones is necessary to add a constant term to the linear fit
        res(:,i) = r;
    end
    
    res = [temp(:,1:3) res];
end
