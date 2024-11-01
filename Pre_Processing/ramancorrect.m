function [correctedshifts,interpolatedInt] = ramancorrect(ramanshift,standardInt,standardSamples,sampleInt)
    %Written by Connor McNairn
    %
    %Corrects Raman Shifts of measured data to specified list of known peaks in
    %standard samples. Function finds peaks in each standard sample
    %corresponding to known peaks and then fits a Lorentzion function to
    %each peak to determine the true peak location. True peak location is
    %used with knownpeaks to fit a 3rd order polynomial that corrects for
    %shifts in the Raman spectrum
    %
    %Inputs:
    %
    %   raman shift - Vector of Raman shifts include in measurement. 
    %
    %   int - Matrix of intensity for each standard sample included in
    %      measurement. Each row corresponds to one standard sample
    %
    %   standardSample - Names of standard samples used. Each entry
    %       corresponds to a row in int
    %
    %Outputs:
    %   correctedshifts - Raman spectrum corrected to known peaks of specified
    %      standard samples
    %
    %   interpolatedShifts - Interpolated intensity values after shift correction
    %       at original Raman shifts
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Standard Sample names and known peaks taken from:
    %https://www.chem.ualberta.ca/~mccreery/ramanmaterials.html
    standardPeaksStrings = ["Silicon","Polystyrene","Tylenol"];
    standardPeaks = {520.2,...
        [620.9,1001.4,1031.8,1155.3,1583.1,1602.3],...
        [651.6,710.8,797.2,834.5,857.9,968.7,1105.5,1168.5,1236.8,1278.5,1371.5,1561.6,1648.4]};
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if nargin < 3
        error('Error: ramanshift, standardInt, standardSamples are required inputs')
    end

    if nargin < 4
        sampleInt = standardInt;
    end


    if size(standardInt,2) ~= size(ramanshift,2) 
        error('Raman shift vector must be same length as intensity matrix')
    elseif size(standardInt,1) ~= length(standardSamples)
        
        error('Number of rows of int and must match number of standard samples specified by standardSamples')

    elseif isstring(standardSamples) ~=1 

        error("standardSamples must be an array of strings")
       

    elseif size(sampleInt,2) ~= size(ramanshift,2) &&  size(sampleInt,2) ~= size(standardInt,2)
        error('Sample intensity matrix must be same length as raman shift vector and standard sample intensity matrix')
    end

    knownpeaks = cell(1,length(standardSamples));
    for i = 1:length(standardSamples)
        for j = 1:length(standardPeaksStrings)
            if standardSamples(i) == standardPeaksStrings(j)
                knownpeaks{i} = standardPeaks{j};
                break;
            elseif j == length(standardPeaksStrings)
                error(standardSamples(i) + " does not match any of the following standard samples: " + strjoin(standardPeaksStrings + repelem([",",""], [length(standardPeaksStrings)-1, 1])))
            end
        end
    end

    tempsize = cellfun(@size,knownpeaks,'UniformOutput',false);
    peaksize = sum([tempsize{:}]) - length(knownpeaks);
    ft = fittype('y0 + (2*A./pi).*(w./(4*(x - xc).^2 + w.^2))');
    truepeak = zeros(1,peaksize);
    truepeakerror = zeros(1,peaksize);
    peakloc = zeros(1,peaksize);
    knownpeaksvector = zeros(1,peaksize);


    u = 1;
    for i = 1:size(standardInt,1)

        [pks,locs,widths,prominances] = findpeaks(standardInt(i,:),'MinPeakHeight',max(standardInt(i,:))*0.03,'MinPeakWidth',2);
        peakraman = ramanshift(1,locs);
        for j = 1:length(knownpeaks{i})
          
            logicalvect = abs(peakraman-knownpeaks{i}(j))<5;
            if sum(logicalvect) ~= 1
                error(standardSamples(i)+" measured peaks do not match known peaks from standard sample")
            end
            peakloc(u) = locs(logicalvect);
      
            ind_range = peakloc(u)-round(widths(logicalvect))+2:peakloc(u)+round(widths(logicalvect))-2;
           
            y0 = pks(logicalvect) - prominances(logicalvect);
            w = ramanshift(1,ind_range(end))-ramanshift(1,ind_range(1));
            A = sum(standardInt(i,ind_range(1):ind_range(end))-y0);
            xc = ramanshift(1,peakloc(u));

            [ramanfit,gof,fit_output1] = fit(ramanshift(1,ind_range)',standardInt(i,ind_range)',ft,'StartPoint',[A,w,xc,y0]);
            
            truepeak(u) = ramanfit.xc;


            %Calculates jacobian of each fit
            J1 = fit_output1.Jacobian;
            %Epsilon is covariance matrix, middle line of matrix are the variances of
            %each parameter, to get parameter error take square root
            alpha = J1'*J1;
            epsilon = inv(alpha)*(gof.sse/gof.dfe);
            truepeakerror(u) = sqrt(epsilon(3,3));
            knownpeaksvector(u) = knownpeaks{i}(j);
            u = u+1;
        end
    end

    warning('off','all');
    diff = knownpeaksvector-truepeak;
    %f = fit(truepeak',diff','poly3','Weights',truepeakerror2'.^(-2));
    ramanfit = fit(truepeak',diff','poly3');
    warning('on','all');
    correctedshifts = ramanshift(1,:) + ramanfit(ramanshift(1,:))';
    interpolatedInt = spline(correctedshifts,sampleInt,ramanshift(1,:));
    
    
end