function [ intLum, meanLum, errs ] = getImageIntensity( im, cellsize, topt )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%% Inputs and Defaults Options
% 
% if ~exist('dataDir','var') || isempty(filestub)
%     filestub = 'mucin';
% end
% 
% if ~exist('channel','var') || isempty(channel)
%     channel = 1;
% end
% 
% if ~exist('rescale','var') || isempty(rescale)
%     rescale = 1;
% end
% 
% %% Load images
% 
% tifstruct = dir([pwd,'\*.jpg']); %Only want .tifs from FIJI output
% 
% if ~exist('totFrames','var') || isempty(totFrames)
%     totFrames = numel(tifstruct);
% end

% intLum = zeros(totFrames,1);

%% calcs 

lnoise = 2;
kernel_size = 5;
h = fspecial('gaussian', kernel_size, lnoise);

if numel(size(im))>3
    fixedim = squeeze(sum(im,4));
else
    fixedim = im;
end

for j = 1:size(fixedim,3)

    thisim = fixedim(:,:,j);
    thisim = imfilter(thisim,h);
    subx = 1:cellsize:size(im,1);
    suby = 1:cellsize:size(im,2);
    
    subx = subx(1:end-1);
    suby = suby(1:end-1);

    if topt
level = graythresh(thisim); %Threshhold to find pixels that actually contain mucin (Ohtsu)
    else
level = multithresh(thisim,2);
level = min(level);
    end

thisim = thisim-level;

for jj = 1:length(subx)-1
        for kk = 1:length(suby)-1
            tempsub = thisim(subx(jj):subx(jj+1),suby(kk):suby(kk+1));
            courseim(jj,kk) = mean(tempsub(:)>0);
        end
    end

intLum(j) = sum(thisim(thisim(:)>0));

N = sum(thisim(:)>0);

meanLum(j) = mean(thisim(thisim(:)>0));
stdLum(j) = nanstd(thisim(:))/sqrt(numel(thisim(:)));
if numel(size(im))>3
    errs(j) =std(courseim(:));
else
    errs(j) = std(courseim(:));
end
% errs(j) = std(courseim(:));


end



% dataOut(:,:,1) = layerLum;
% dataOut(:,:,2) = intLum;

end

