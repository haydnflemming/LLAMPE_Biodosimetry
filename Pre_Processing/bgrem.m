% -- SNIP BACKGROUND REMOVAL ----------------------------------------------
% Written by Harry Allen
% 'filterOrder' 4 and 6 added by Cristian Ciobanu
%
% Estimates fluorescence background using a sensitive non-linear iterative
% peak (SNIP) clipping algorithm.

function [bgremd, bg, smooth] = bgrem(spectra,window,order,sgw,filterOrder)
    %window must be odd
    %bgremd is the data with background removed
    %bg is the background fit to each spectrum
    %so both bgremd and bg have the same size as the spectra data
    [rows, cols] = size(spectra);
    %index 1 in spectra will be a = w+1 in padded version
    w = floor(window/2);
    a = w+1;
    %index end in spectra will be b = end + w in padded version
    b = cols+w;
    %create padded version of spectra
    padded = zeros(cols + 2*w,rows);
    padded(a:b,:) = spectra.';
    %padded region is reflected difference
    for i = 1:w
        padded(a-i,:) = 2.*padded(a,:) - padded(a+i,:);
        padded(b+i,:) = 2.*padded(b,:) - padded(b-i,:);
    end
    %smooth padded dataset
    padded(:,:) = sgolayfilt(padded,order,sgw);
    smooth = padded(a:b,:).';
    %now apply SNIP
    if filterOrder == 2
        for j = w:-1:1
            A = padded(a:b,:);
            B = (padded(a-j:b-j,:) + padded(a+j:b+j,:))./2;
            
            padded(a:b,:) = min(A, B);
        end
    elseif filterOrder == 4
        for j = w:-1:1
            A = padded(a:b,:);
            B = (padded(a-j:b-j,:) + padded(a+j:b+j,:))./2;
            C = (-padded(a-j:b-j,:) + 4*padded(a-j/2:b-j/2,:) + 4*padded(a+j/2:b+j/2,:) - padded(a+j:b+j,:))./6;
            
            padded(a:b,:) = min(A, max(B, C));
        end
    elseif filterOrder == 6
        for j = w:-1:1
            A = padded(a:b,:);
            B = (padded(a-j:b-j,:) + padded(a+j:b+j,:))./2;
            C = (-padded(a-j:b-j,:) + 4*padded(a-j/2:b-j/2,:) + 4*padded(a+j/2:b+j/2,:) - padded(a+j:b+j,:))./6;
            D = (padded(a-j:b-j,:) - 6*padded(a-2*j/3:b-2*j/3,:) + 15*padded(a-j/3:b-j/3,:))./20 ...
                + (padded(a+j:b+j,:) - 6*padded(a+2*j/3:b+2*j/3,:) + 15*padded(a+j/3:b+j/3,:))./20;
            
            padded(a:b,:) = min(A, max(B, max(C, D)));
        end
    end
    
    bg = padded(a:b,:).';
    smooth = smooth - bg;
    bgremd = spectra - bg;
end