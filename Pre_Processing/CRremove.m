%% ---- COSMIC RAY REMOVER ---- %%
% Written by Amiel Beausoleil-Morrison (09/06/2023)
% based on method described in the paper: 
% Barton, S.J., Henneley, B.M. "An algorythm for the removal of cosmic ray 
% artifacts in spectral data sets," Applied Spectroscopy, 73(8) 893-901 (2019).

function[fixed, cosmic_ray1, cosmic_ray2, cind, cray_num] = CRremove(data, threshold1,...
    threshold2, range, remType)

global meds neighbours stdev

% -------------------------------------------------------------------------
% INPUTS: 
% 1. data - Raman intensity data
% 2. threshold1 - threshold to check if a point is far enough from the
%    median to be identified as a CR (if threshold1*stdev <= absdiff)
% 3. threshold2 - threshold to check neighbouring points around CR's (if
%    prefered removal method is to check each point of the CR)
% 4. range - range to remove around each CR if this is the prefered method
% 5. remType - 0 if one wants to remove the CR width using a threshold2
%    check and 1 if one wants to remove a range around each CR

% OUTPUTS:
% 1. fixed - spectra following CR removal 
% 2. cosmic_rays - logical matrix containing the CR's
% 3. cind - spectra with CR's indices
% 4. cray_num - number of CR's found 

% -------------------------------------------------------------------------

spectraNum = size(data, 1);
binNum = size(data, 2);

fixed = data;

% calculate the normalized covariance between each spectrum with the others
% this formula is from the paper
for i = 1:spectraNum
    for j = 1:spectraNum
        covariance(i, j) = (dot(fixed(i, :),fixed(j, :)).^2)./...
            (dot(fixed(i, :), fixed(i, :)).*dot(fixed(j, :), fixed(j, :)));
    end
end

% order the covariance (col by col) such that the 
% row 1: 1 2 3 4 ... (spectral indices with cov = 1)
% first largest covariance
%       .
%       .
%       .
[~, inds] = sort(covariance, "descend");

% chooses 4 spectra with the highest covariance as nearest neighbours
neighbours = transpose(inds(1:5, :));

% finding the median and standard deviation of each set of nearest neighbours 
% which the test spectrum will be compared to (median of the 4 nearest neighbours)
meds = [];
for i = 1:spectraNum
    ners = fixed(neighbours(i, :), :);
    std_ner(i, :) = std(ners);
    meds = [meds ; median(ners)];
    % meds(i, :) = median(ners);
end

% mean standard deviation for each set of nearest neighbours
stdev = mean(transpose(std_ner));

% if the abs diff between the test data point and the median of the nearest
% neighbour set, then the point is identified as a cosmic ray
for i = 1:spectraNum 
    for j = 1:binNum
        absdiff(i, j) = abs(fixed(i, j)-meds(i, j));
        if absdiff(i, j) >=  threshold1*stdev(i)
            cosmic_ray1(i, j) = 1;
        else 
            cosmic_ray1(i, j) = 0;
        end
    end
end

[n, m] = find(cosmic_ray1); % indices of the cosmic ray locations
cray_num1 = nnz(cosmic_ray1);

if remType == 0
    % checking the points directly preceding and following the CR's located
    % with a different threshold to remove all artefacts of the CR
    cosmic_ray2 = cosmic_ray1;
    for i = 1:cray_num1 
        if m(i) <= binNum -1 && absdiff(n(i), m(i)+1) >=  threshold2*stdev(n(i))
            cosmic_ray2(n(i), m(i)+1) = 1;
        elseif m(i) >= 2 && absdiff(n(i), m(i)-1) >= threshold2*stdev(n(i))
            cosmic_ray2(n(i), m(i)-1) = 1;
        end
    end
    
    cray_num2 = nnz(cosmic_ray2-cosmic_ray1);
    cray_num = cray_num1+cray_num2;
    cosmic_rays = cosmic_ray2;
    [n2, m2] = find(cosmic_ray2);
    
    % remove entire affected area for each cosmic ray
    for i = 1:size(n2)
        fixed(n2(i), m2(i)) = nan;
    end

elseif remType == 1
    for j = 1:range
        for i = 1:size(n)
            if m(i) <= range
                fixed(n(i), 1:m(i)+j) = nan;
            elseif binNum - m(i) <= range
                fixed(n(i), m(i)-j:binNum) = nan;
            else
                fixed(n(i), m(i)-j:m(i)+j) = nan;
            end
        end
    end
    cosmic_ray2 = [];
    cray_num = cray_num1;
    cosmic_rays = cosmic_ray1;
end

ind = isnan(fixed);
[cind, ~] = find(isnan(fixed));

% set the removed area to the median of the nearest neighbour group
for i = 1:spectraNum
    for j = 1:binNum
        if ind(i, j) == 1
            fixed(i, j) = meds(i, j);
        end
    end
end

end