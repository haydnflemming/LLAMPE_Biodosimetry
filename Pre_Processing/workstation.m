% -- SPECTRAL PREPROCESSING WORKSTATION -----------------------------------
% Written by Cristian Ciobanu
%
% This is a workstation for preprocessing and analyzing Raman spectra.

close all;
clear all;
clc;

% -- GLOBAL INTERMEDIATE PREPROCESSED SPECTRA -----------------------------
global csr sgf MSC SNV bgr bge hqt nrm avg hqr mdn

% -- GLOBAL PREPROCESSING SETTINGS ----------------------------------------
global avgNum sgfOrder sgfWindow snipOrder snipWindow

sgfOrder = 2;
sgfWindow = 11;

avgNum = 12;

snipWindow = 111;
snipOrder = 2;

% -- NICKNAMES ------------------------------------------------------------
% lysate_feb2017_OO (ocean optics), lysate_nov2017_T (tornado data), plasma_31jul2019_T (old 2s plasma),
% plasma_15may2021_T (10s plasma), lysate_oct2017_T (0v50Gy), lysate_may2018_T (mock)
% plasma_18jun2021_T (15s plasma), plasma_30jul2019_T (old 2s plasma)

[shifts, data] = importData('lysate_nov2017_T');

% -- PREPROCESSING OPTIONS -------------------------------- 6,1,8,2,4,5,6,7
% 1 - Cosmic Ray Removal         2 - Savitzky-Golay Smoothing
% 3 - MSC Correction             4 - SNV Correction
% 5 - Background Removal         6 - Outlier Removal
% 7 - Normalization              8 - Shuffle Averaging
% 9 - Median

[all, train, test] = prepSet(data, 17, [0,5], -1, [6,1,8,2,4,5,6,7]);

function [all, train, test] = prepSet(data, donor, doses, ser, preprocessing)
    all = 0;
    
    if donor == -1
        donor = unique(data(:,1))';
    end
    
    if doses == -1
        doses = unique(data(:,2))';
        if ismember(999,doses)
            doses = doses(1:end-1);
        end
    end
    
    compiled_to_shuffle = cell(size(doses));
    
    for i = 1:size(donor,2)
        for j = 1:size(doses,2)
            [info, crnt] = dataPicker(data, donor(i), doses(j), ser);
            [crnt, removed] = apply(crnt, preprocessing);
            
            if ismember(8,preprocessing) || ismember(9,preprocessing)
                numRows = ones(size(crnt,1),1);
                info = [donor(i)*numRows doses(j)*numRows 999*numRows];
            else
                for k = 1:size(removed,2)
                    info(removed{k},:) = [];
                end
            end
            
            if all == 0
                all = [info crnt];
            else
                all = [all; info crnt];
            end
            
            if compiled_to_shuffle{j} == 0
                compiled_to_shuffle{j} = [info crnt];
            else
                compiled_to_shuffle{j} = [compiled_to_shuffle{j}; info crnt];
            end
        end
    end
    
    all(:,2) = all(:,2) + 1;
    [train, test] = combine(compiled_to_shuffle);
    train(:,2) = train(:,2) + 1;
    test(:,2) = test(:,2) + 1;
end

function [crnt, removed] = apply(data, operations)
    global csr sgf MSC SNV bgr bge hqt nrm avg hqr mdn
    
    global avgNum sgfOrder sgfWindow snipOrder snipWindow
    crnt = data;
    
    removed = {};
    
    for i = operations
        if i == 1
            [crnt, ~] = removeCosmicRays(crnt,1.59,2,15);
            csr = crnt;
        elseif i == 2
            crnt = sgolayfilt(crnt',sgfOrder,sgfWindow)';
            sgf = crnt;
        elseif i == 3
            [crnt, ~] = msc(crnt);
            MSC = crnt;
        elseif i == 4
            crnt = snv(crnt);
            SNV = crnt;
        elseif i == 5
            [crnt, bge, ~] = bgrem(crnt, snipWindow, 2, 15, snipOrder);
            bgr = crnt;
        elseif i == 6
            [crnt, hqr] = hqTest(crnt);
            removed{1,end+1} = hqr;
            hqt = crnt;
        elseif i == 7
            crnt = normr(crnt);
            nrm = crnt;
        elseif i == 8
            crnt = shuffle(crnt);
            crnt = averager(crnt, avgNum);
            avg = crnt;
        elseif i == 9
            crnt = median(crnt);
            mdn = crnt;
        end
    end
end

