% Linear mixed-effect (LME) model covariate adjustment code written by Cristian
% Ciobanu.
% edited 31/07/23
%
% 'res' - Covariate-adjusted data (fit residuals).
%
% 'lmes' - Cell array containing linear mixed effects model for each Raman
% shift.
%
% 'data' - The Raman intensity data (2D matrix), with first three columns containing
% donor ID, dose, and serial number.
%
% 'covariates' - A table (not matrix) containing donor IDs and respective covariate
% data. Note that the columns must be labeled as: Donor, Gender, Age,
% Smoker, WBC, RBC, HGB. The data must be categorical (e.g., for smoker
% status, use a 1 or 0 to indicate if they are smoker or not, instead of
% writing something like "true" or "false").
%
% 'type' - The algorithm used to perform the linear mixed effects model
% fitting. Can either be 'ML' or 'REML', see the MATLAB 'fitlme' page for
% more information.
%
% 'smoker' - 1 if you want to consider smoker as a covariate, and 0 if not.
%
% NOTE: To use multiple linear regression models instead, you just have to
% adjust the code a little bit at the end to use the 'fitlm' function
% instead of 'fitlme'.

function [res, lmes] = covAdjTable(data, covariates, type, dataset)
    lmes = {}; % Will store LME models in here later.

    sizeData = size(data); % Getting size of data. Output format: [rows, columns]
    sizeCov = size(covariates); % Getting size of covariate table. Output format: [rows, columns]
    
    % Here we intialize an empty table that will hold all intensities from
    % 'data' in the last column, and the respective covariate information in
    % the first seven columns.

    if sizeCov(2) == 8 % 'types' should be a 1xN cell array where N is sizeCov(2) + 1.
        types = {'double', 'double','double','double','double','double','double', 'double', 'double'};
    elseif sizeCov(2) == 7
        types = {'double','double','double','double','double','double','double', 'double'};
    elseif sizeCov(2) == 2
        types = {'double', 'double', 'double'};
    elseif sizeCov(2) == 9
        types = {'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double'};
    end
    cov = table('Size',[sizeData(1) sizeCov(2)+1],'VariableTypes',types,'VariableNames',[covariates.Properties.VariableNames 'Intensity']); % Create the table.
    
    % Here we fill the table in with the covariates and intensity data.
    % Finally, the LME models are build with the 'fitlme' function.
    
    for i = 1:sizeData(1)
        donor = data(i,1);
        for j = 1:sizeCov(1)
            if covariates{j,1} == donor
                cov{i,1:sizeCov(2)} = covariates{j,:};
                break;
            end
        end
    end
    
    for i = 4:sizeData(2)-2
        cov{:,sizeCov(2)+1} = data(:,i);
        if dataset == 0
            lme = fitlme(cov, "Intensity~Age+WBC+RBC+HGB+(1|Donor)+(Gender-1|Donor)+(Colour-1|Donor)+(Day-1|Donor)");
        % if smoker == 0 % Include smoker in covariates if desired ('smoker' = 1)
        %     lme = fitlme(cov,"Intensity~Gender+Age+WBC+RBC+HGB+Day+Colour+(1|Donor)", "FitMethod", type);
        % else % Or exclude them with 'smoker' = 0.
        %     lme = fitlme(cov,"Intensity~Gender+Age+Smoker+WBC+RBC+HGB+(1|Donor)", "FitMethod", type);
        elseif sizeCov(2) == 2
            lme = fitlme(cov,"Intensity~Date+(1|Tube)", "FitMethod", type);
        elseif dataset == 16
            % november 2016 set of ocean optics data
            lme = fitlme(cov, "Intensity~Age+WBC+RBC+HGB+(1|Donor)+(Gender-1|Donor)+(Day-1|Donor)+(Medication-1|Donor)+(Smoker-1|Donor)");
        elseif dataset == 17
            % february 2017 set of ocean optics data
            lme = fitlme(cov, "Intensity~Age+WBC+RBC+HGB+(1|Donor)+(Gender-1|Donor)+(Day-1|Donor)+(Smoker-1|Donor)");
        elseif dataset == 23
            lme = fitlme(cov, "Intensity~Age+(Gender-1|Donor)+(Day-1|Donor)+(Colour-1|Donor)+(1|Donor)");
        end
        data(:,i) = residuals(lme); % Replacing the intensity data with the residuals for each column.
        lmes{end + 1} = lme; % Add the LME model to the 'lmes' variable (just there for reference).
    end
    
    res = data; % Output the data in the 'res' variable.
end