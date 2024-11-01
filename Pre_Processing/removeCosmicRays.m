% -- COSMIC RAY REMOVER ---------------------------------------------------
% Written by Harry Allen
%
% Detects and removes cosmic rays from spectra by smoothing over them with
% a Savitzky-Golay filter.

function [fixed,crays,count] = removeCosmicRays(spectra,C,order,framelen)

fixed = zeros(size(spectra));
fixed(:,:) = spectra;
weights = ones(1,framelen);

med = median(fixed);
absdiff = abs(fixed-med).*fixed;
%smad = max(absdiff(absdiff < max(absdiff)));
sort_absdiff = sort(absdiff,'descend'); % new
smad = sort_absdiff(3,:); %new
cosmic_rays = absdiff > C*smad;
count = nnz(cosmic_rays);
crays = false(size(cosmic_rays));
crays(:,:) = cosmic_rays;
nf = any(cosmic_rays,2);
midpoint = floor(framelen/2 + 1);
stop = midpoint/2;
w = 1;
while any(nf)
    weights(1,midpoint-w:midpoint+w) = 0.0001; % 0.0001
    sg = sgolayfilt(fixed(nf,:),order,framelen,weights,2);
    inds = find(nf);
    for i = 1:length(inds)
        fixed(inds(i),cosmic_rays(inds(i),:)) = sg(i,cosmic_rays(inds(i),:));
    end
    med(1,:) = median(fixed);
    absdiff(:,:) = abs(fixed-med);
    smad(1,:) = max(absdiff(absdiff < max(absdiff)));
    cosmic_rays(:,:) = absdiff > C*smad;
    nf(:) = any(cosmic_rays,2);
    w = w+1;
    if w > stop
        break
    end  
end
end