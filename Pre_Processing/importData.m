% -- IMPORT FUNCTION ------------------------------------------------------
% Written by Cristian Ciobanu
%
% 'Database' of datasets for easy access in workstation.
%
% If you'd like to use this function to import data, you will need to
% replace all the if/elseif statements to point to your own data. I'm
% leaving all mine here to serve as a template.

function [shifts, data] = importData(nickname)
    if strcmp(nickname, 'lysate_feb2017_OO')
        cd('C:\Users\cristi\Desktop\LLAMPE\Data\Blood Lysate\Nov24_2016-Feb16_2017-OceanOptics');
        data = load('112416_021617_LysedBlood_Raw_OceanOptics_AllSpectra.mat');
        [shifts, data] = formatData(cell2mat(struct2cell(data)), 201:967);
    elseif strcmp(nickname, 'lysate_nov2017_T')
        cd('C:\Users\cristi\Desktop\LLAMPE\Data\Blood Lysate\Oct-Nov_2017-Tornado');
        data = load('Oct_Nov_2017_LysedBlood_Raw_Tornado_AllSpectra.mat');
        [shifts, data] = formatData(cell2mat(struct2cell(data)), 401:1601);
    elseif strcmp(nickname, 'plasma_31jul2019_T')
        cd('C:\Users\cristi\Desktop\LLAMPE\Data\Blood Plasma\Blood_190731');
        data = load('190731_BloodPlasma_Raw_Tornado_AllSpectra.mat');
        [shifts, data] = formatData(cell2mat(struct2cell(data)), 671:1732);
    elseif strcmp(nickname, 'plasma_15may2021_T')
        cd('C:\Users\cristi\Desktop\LLAMPE\Data\Blood Plasma\Blood_2021_05_15');
        data = load('210515_BloodPlasma_BSL+NIST_Tornado_AllSpectra.mat');
        [shifts, data] = formatData(cell2mat(struct2cell(data)), 671:1732);
    elseif strcmp(nickname, 'lysate_oct2017_T')
        cd('C:\Users\cristi\Desktop\LLAMPE\Data\Blood Lysate\Oct_2017_0Gy_vs_50Gy_Tornado');
        data = load('Oct_2017_LysedBlood_0Gy_vs_50Gy_Raw_Tornado_AllSpectra.mat');
        [shifts, data] = formatData(cell2mat(struct2cell(data)), 401:1601);
    elseif strcmp(nickname, 'lysate_may2018_T')
        cd('C:\Users\cristi\Desktop\LLAMPE\Data\Blood Lysate\May30_2018_Mock_Tornado');
        data = load('May_30_2018_LysedBlood_Mock_Raw_Tornado_AllSpectra.mat');
        [shifts, data] = formatData(cell2mat(struct2cell(data)), 401:1601);
    elseif strcmp(nickname, 'plasma_18jun2021_T')
        cd('C:\Users\cristi\Desktop\LLAMPE\Data\Blood Plasma\Blood_2021_06_18');
        data = load('210618_BloodPlasma_BSL+NIST_Tornado_AllSpectra.mat');
        [shifts, data] = formatData(cell2mat(struct2cell(data)), 671:1732);
    elseif strcmp(nickname, 'plasma_30jul2019_T')
        cd('C:\Users\cristi\Desktop\LLAMPE\Data\Blood Plasma\Blood_190730');
        data = load('190730_BloodPlasma_Raw_Tornado_AllSpectra.mat');
        [shifts, data] = formatData(cell2mat(struct2cell(data)), 671:1732);
    end
end

function [shifts, data] = formatData(raw, fingerprintRegion)
    shifts = raw(1,4:end);
    shifts = shifts(1,fingerprintRegion);
    data = raw(2:end,:);
    class = data(:,1:3);
    data = data(:,4:end);
    data = data(:,fingerprintRegion);
    data = [class, data];
end