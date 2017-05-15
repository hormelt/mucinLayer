function [objs_ring HM2] = houghring_rp(img, objsize, thresh, Houghthresh, rrange, fitstr)

% houghring_rp : Uses a Hough transform to find the center of bright rings 
%                (not filled circles) objects in images
%
% Uses Hough Transform method (standard) from (Author?)
% http://basic-eng.blogspot.com/2006/02/hough-transform-for-circle-detection.html
% Corrects error in this of not taking +/- square roots in determination of
% "a" -- see approach in
% http://www.mathworks.com/matlabcentral/fileexchange/9898
% by Tan Chin Luh, 2006
%
% first, filters image, then..
% Thresholds, then looks for objects in binary thresholded image.  Keeps
% only objects with >2*pi*r pixels -- eliminating small noise specks.  RP
% 15Oct09 -- this seems to work!
% Then, implements a Hough transformation, then
%   calls fo4_rp (object finding) to find the centers in the
%   parameter space.
% Also weight "voting" by the brightness of the thresholded image.
%
% Suggestion:  First run viewHoughseries with a known estimate
% of the ring radius to determine good threshold and object size
% parameters.
%
% INPUTS
% img = image to locate objects within
% objsize : size in pixels of objects to find. Used to determine:
%    ste = structuring element with diameter approximately particle size
%    and image-space bandpass filter
% thresh = number in [0, 1] that sets the local max threshold
% Houghthresh = number in [0, 1] that sets the local max threshold for
%    object detection in the Hough parameter space.  (Use radius/2 for
%    object size for detection)
% rrange = [rmin rmax] 2 element array of min and max radii to search (px)
% fitstr : (Optional) string that selects the fitting option (Gaussian, 
%    non-linear Gaussian fit (default) or centroid) for fo4_rp.    
%    For centroid, enter 'centroid'.
%    'none' for no fitting (leave objs empty; just return Hough image HM2)
%
% OUTPUT
% Centers are returned in the matrix 'objs' which has the following form:
%  [x;
%   y;
%   mass;
%   radius;  ** Note -- not particle ID as in the usual object matrix
%   frame; 
%   trackid]
%
% Also returns the Hough parameter space matrix; sometimes useful to look
% at;
%
% frame and trackid are set to zero for all objects by houghring_rp() and must be
% dealt with by the functions calling houghring_rp().
%
% Raghuveer Parthasarathy
% October 13, 2009
% July 13, 2010 -- modifications for speed, and quadratic fit
%    to best radius value
% Dec. 7, 2010 -- Alter object size for fitting Hough-space peak
% Last modified Aug. 8, 2011 

% Fitting option; default is nonlinear Gaussian
if ~exist('fitstr', 'var') || isempty(fitstr)
    fitstr = 'nonlineargauss';
end

showplots=false;  % for debugging -- plot things.

% make img processable
img = double(img);

% now do bandpass filter -- use Grier et al bpass.m
noisesize = 1;  % size of noise, px
filtimg = bpass(img,noisesize,objsize);

if showplots
    figure; imagesc(img); colormap(gray); title('1 original image')
    figure; imagesc(filtimg); colormap(gray); title('2 bandpass filtered')
end

% Threshold image
% compute noise level and threshold
[hs, bins] = hist(filtimg(:),100);
ch = cumsum(hs);
ch = ch/max(ch);
noiseind = find(ch > thresh); %
noiseind = noiseind(2); % The index value below which "thresh" fraction
% of the pixels lie
isthresh = (filtimg > bins(noiseind));  % Boolean array: above threshold?
[L,num] = bwlabel(isthresh,4);  % look for connected objects (4-nn)
if showplots
    figure; imagesc(L)
end
nj = zeros(1,num);
threshobj = 2*pi*min(rrange);  % we'll only keep connected objects 
    % with more than this many pixels -- i.e. the (min) circumference of a
    % full circle 1 px thick
LSTATS = regionprops(L, 'area', 'PixelIdxList'); % region areas and pixelIDs
badpx = [];
for j=1:num
    nj(j) = LSTATS(j).Area; % number of pixels in object j
    % SLOW:  nj(j) = sum(sum(L==j)); % number of pixels in object j
    if (nj(j) < threshobj)
        badpx = [badpx; LSTATS(j).PixelIdxList];  % pixels of small regions
        % SLOW:  L = L.*(L~=j);
    end
end
% remove bad (small) regions
if ~isempty(badpx)
    L(badpx) = 0;
end

threshimg = filtimg.*(L>0);  % image, with the below-threshold pixels black
if showplots
    figure; imagesc(L)
    figure; plot(nj, 'ko')
    h3=figure; imagesc(threshimg); colormap(gray); title('3 thresholded')
end


% ------------------------------------------------------
% Hough transformation
% It should be possible to vectorize this, but it's too painful to figure
% out how to "wrap" the ind variables to account for the different R
% values.  Could adopt the solution at "http://basic-eng.blogspot.com," but
% this involves the "find" function, which I don't like.
% Simply loop.

if length(rrange)>1
    nR = rrange(2)-rrange(1)+1;  % number of radius values to examine
    R = rrange(1):rrange(2);
else
    nR = 1;
    R = rrange;
end

[sy,sx]=size(threshimg);
[y,x]=find(threshimg>0);  % all nonzero points
totalpix = length(x);

b = 1:sy;
y = repmat(y',[sy,1]);
x = repmat(x',[sy,1]);

HM2 = zeros(sy,sx,nR); % allocate memory

for j = 1:nR,
    b1 = repmat(b',[1,totalpix]);  % sloppily re-define this, since altered later (gooda1)
    b2 = b1;  % for other square root
    R2 = R(j)*R(j);
    HM = zeros(sy*sx,1);  % allocate memory for the Hough matrix

    % all values of the 'a' parameter.  Time-consuming
    tempa = R2 - (y - b1).^2;
    gooda = tempa>0;  % the negative a's will get discarded as imaginary
    sq = zeros(size(y));
    sq(gooda) = sqrt(tempa(gooda));
    % sq = sqrt(R2 - (y - b1).^2);
    % a1 = round(x - sq);
    % a2 = round(x + sq);
    a1 = (x - sq);
    a1(~gooda) = -1;
    a2 = (x + sq);
    a2(~gooda) = -1;

    % Removing all the invalid value in matrices a and b
    %    Use 1 and sx-1 as the bounds rather than 0 and sx to avoid
    %    checking for rounding issues
    gooda1 = (imag(a1)==0 & a1>1 & a1<(sx-1));
    gooda2 = (imag(a2)==0 & a2>1 & a2<(sx-1));
    xx1 = x(gooda1);
    yy1 = y(gooda1);
    xx2 = x(gooda2);
    yy2 = y(gooda2);
    b1 = b1(gooda1);
    a1 = round(a1(gooda1));
    b2 = b2(gooda2);
    a2 = round(a2(gooda2));
    ind1 = sub2ind([sy,sx],b1,a1);
    ind2 = sub2ind([sy,sx],b2,a2);
    ind = [ind1; ind2];   % the indices of all the valid points in the ab plane
    
    % For weighting: find the indices in the real images
    wind1 = sub2ind([sy,sx],yy1,xx1);
    wind2 = sub2ind([sy,sx],yy2,xx2);
    wind = [wind1; wind2];   % the indices of all the valid points in the
    %  xy (image) plane, for weighting.

    % Reconstruct the Hough Matrix
    % val = ones(length(ind),1);  % Use this line for no weighting
    val = 1.*(threshimg(wind));  % Weighted vote (weight by intensity)
    data=accumarray(ind,val);
    HM(1:length(data)) = data;
    HM2(:,:,j) = reshape(HM,[sy,sx]);
end

if showplots
    figure; hold on
    for j=1:nR
        imagesc(HM2(:,:,j));
        pause(0.25)
    end
end

% determine circle radius -- only looks for one (max) value!
H = zeros(1,nR);
for cnt = 1:nR
    H(cnt) = max(max(HM2(:,:,cnt)));
end
[maxval, maxind] = max(H);
if ((nR>1) && (maxind>1) && (maxind<nR))
    % find peak in H(R) -- very simple parabolic fit
    m = (maxind-1):1:(maxind+1);
    p = polyfit(R(m), H(m), 2);
    bestR = -p(2)/2/p(1);
else
    bestR = R(maxind);
end
if showplots
    figure; plot(R, H, 'ko-');
    hold on
    plot(bestR*[1 1],[min(H) max(H)],'r:');
end

% [B,A] = find(HM2(:,:,maxind)==maxval);
% if showplots
%     figure; imagesc(double(HM2(:,:,maxind))); shading interp; title('HM2')
%     figure(h3); hold on;
%     plot(mean(A),mean(B),'xm')
%     text(mean(A),mean(B),num2str(bestR),'color','green')
% end

if ~strcmp(lower(fitstr), 'none');
    % At the radius (crudely) determined above, call fo4_rp to find the
    % peaks in the Hough matrix
    % Use fixed object size, equal to the object size used in the original filtering
    %   -- may need to alter / think further
    [objs_ring] = fo4_rp(HM2(:,:,maxind), objsize, Houghthresh, fitstr);
    % HMobjsize = round(bestR/2.0);  %objsize;
    % [objs_ring] = fo4_rp(HM2(:,:,maxind), HMobjsize, Houghthresh, fitstr);
    objs_ring(4,:) = bestR;
else
    objs_ring = [];
end

