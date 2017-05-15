% filtthreshimg.m
% 
% filter and threshold (2D) image
% uses same procedure as in f04_rp.m
%
% ** To apply to 3D images we'll need to modify bpass.m ** 
%
% uses bpass.m (filtering)
%
% INPUTS
% img = image (2D) to locate objects within
% objsize : size in pixels of objects to find. Used to determine:
%    ste = structuring element with diameter approximately particle size
%    and image-space bandpass filter
% thresh = number in [0, 1] that sets the local max threshold.  Keep pixels
%    with intensity at the thresh "percentile"
%
% OUTPUTS
% ftimg = filtered, thresholded image (double precision array)
%
% Suggestion: use testthresh to determine good objsize, thresh parameters
%
% Raghuveer Parthasarathy, June 29, 2010
% Last modified: August 8, 2010 (RP)

function [ftimg] = filtthreshimg(img, objsize, thresh)

if length(size(img))>2
    disp('ERROR! works only for 2D images; need to modify bpass.m for 3D images.');
    disp('Press Control-C');
    pause
end

showplots=false;  % for debugging -- plot things.

% make img processable
img = double(img);

% now do bandpass filter (if objsize > 0) -- use Grier et al bpass.m
if objsize>0
    noisesize = 1;  % size of noise, px
    filtimg = bpass(img,noisesize,objsize);
else
    filtimg = img;
end

if showplots
    figure(1)
    imagesc(img); colormap(gray); title('1 original image')
    figure(2)
    imagesc(filtimg); colormap(gray); title('2 bandpass filtered')
end

% now compute noise level and threshold
[hs, bins] = hist(filtimg(:),100);

ch = cumsum(hs);
ch = ch/max(ch);
noiseind = find(ch > thresh); %
noiseind = noiseind(1); % The index value below which "thresh" fraction
% of the pixels lie.  (Originally noisind(2), but this can lead to errors
% if there is only one bin above threshold.)

% keep pixels with >threshold intensity; set all others to zero (black) 
ftimg = filtimg.*(filtimg > bins(noiseind));  
