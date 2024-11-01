% -- DATA PICKER ----------------------------------------------------------
% Written by Cristian Ciobanu
%
% Given a set of data 'original' containing donor, dose, and serial number
% information in its first 3 columns, the function will return all data
% that satisfies the 'donor', 'dose', and 'serial' specified.

function [class, data] = dataPicker(original, donor, dose, serial)
    picked = zeros(1,size(original,2));
    for i = 1:size(original,1)
        if (ismember(original(i,1), donor) || ismember(-1, donor)) && (ismember(original(i,2), dose) || ismember(-1, dose)) && (ismember(original(i,3), serial) || ismember(-1, serial))
            picked = [picked; original(i,:)];
        end
    end
    picked = picked(2:end,:);
    class = picked(:,1:3);
    data = picked(:,4:end);
end