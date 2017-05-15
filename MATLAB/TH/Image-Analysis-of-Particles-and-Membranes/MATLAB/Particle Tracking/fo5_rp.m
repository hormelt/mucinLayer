% fo5_rp.m : finds objects in one 2D image
%    Determines neighborhoods in which to apply precise particle 
%    localization algorithms, calls these algorithms, and returns particle
%    positions.
% Finds neighborhoods by processing the image by spatial filtering or 
%    gradient voting (both optional), finding local maxima of the processed 
%    image, and thresholding. 
% Then refines positions (using the unprocessed image) with one of 
%    several possible methods noted below
% Calls bpass.m (Grier et al.) for spatial filtering
% Calls gradientvote.m for gradient line voting
% Calls calcthreshpts.m for local max. calculation, thresholding
%
% Note: Replaces fo4_rp.m � different input parameters, input order!
%
% Inputs:
% img = image to locate objects within
% processopt = option for processing to find neighborhoods
%     processopt = 'spatialfilter': spatial filter [default],
%     processopt = 'gradientvote': gradient line voting
%     processopt = 'none', or empty matrix: no filtering (but still need to
%        indicate neighborhood size using processparam)
% processparam:  Parameters for processing
%    If spatially filtering, 2 elements, both a measure of, roughly, the 
%    size of the object to find:
%      (1) bpfiltsize, used to determine the low-pass filter size sent to 
%        bpass.m.  "0" indicates no filtering.
%      (2) nsize = size of the "neighborhood" around each particle, used to
%        make the structuring element for dilation (local max finding) and 
%        to set the array size for single-particle localization.
%      If one element, use the same value for both
%    If using gradient voting, 1-3 elements:
%      (1) gradobjsize  : object size, px.  (Shouldn't be smaller than the object;
%              could be larger).  Gradient votes are counted along the 
%              gradient line +/- gradobjsize pixels from each point if 
%              graddir=0; only in one direction if graddir = +/- 1.
%              Also uses this value for the neighborhood size.
%      (2) graddir  : option for the relevant direction of the intensity gradient
%              0 : [default] gradient vote in both directions
%              +1 or -1 : count votes only for pixels in the direction of
%              positive or negative intensity gradient, respectively.
%      (3) grthresh : (optional) gradient magnitude threshold [0-1).  If >0, allow
%              only points in the >=grthresh fraction of gradient
%              magnitudes to "vote" (Default 0)
% thresh : intensity threshold 
%    *Various options,* inferred from the form of the input
%    (1) if a number in [0, 1), sets the local max intensity threshold, 
%        keeping all points with intensity above the thresh*100 percentile
%    (2) if a number <0, keeps pixels with intensity > thresh*std.dev. above
%          the median intensity
%    (3) if a number >= 1, keeps brightest "thresh" number of particles
%          by finding the "thresh" brightest regions and allowing only one
%          maximum in each.  (Slow -- due to the enforcement of 1 particle
%          per neighborhood -- see "try1pernhood" below)
% fitstr : String (or array of strings) that selects the fitting option 
%     fitstr can be a character array (i.e. string) such as 'radial,' in which
%        case it indicates the center-finding method.
%        OR fitstr can be a cell array such as {'radial'; 'momentcalc'} --
%        note the curly braces -- in which case the first element is the
%        center-finding method and the second is the orientation-finding
%        method (for elliptical particles}
%      Center-finding:
%      -- [default] 'radial'.  Radial-symmetry based fit.  Fast, accurate -- see
%         notes July-August 2011, and RP's Nature Methods paper, June 2012.
%         If the number of points in the image is <10, call radialcenter.m .
%         Otherwise, call radialcenter_stk.m (stack input, avoiding 
%         redundant grid calculations for speed).
%      --  'gaussmle' , radially symmetric 2D Gaussian, fit by maximum
%         likelihood estimation.  Slow, most accurate.
%      --  'nonlineargauss' , radially symmetric 2D Gaussian nonlinear fit.  
%         Slow, accurate. Fit to intensity = offset  A*exp(-(x^2+y^2)/2/sigma^2);
%        The "intensity" (row 3 of the object matrix) is pi*A*sigma_x*sigma_y
%      -- 'lineargauss' , 2D Gaussian, linear fit (i.e. parabolic 
%         fit to log(intensity).  Fast, moderate accuracy � dangerous
%      -- 'centroid'.  Centroid fit.  Least accurate (heavily biased); 
%         may be necessary for particles with saturated image intensities.
%      -- 'weightedlineargauss' , Linearized Gaussian fit, weighted 
%         for noise (Stephen M. Anthony, Steve Granick -- see Langmuir paper)
%         DON'T USE THIS
%      Orientation-finding:
%      -- [default] 'none'.  Don't assess orientation.  (E.g. usual
%         point-particles or spherical colloids)
%      -- 'momentcalc' call simpleellipsefit.m for simple determination of
%         the ellipsoid corresponding to the intensity covariance matrix
% try1pernhood : (optional; default false).  Input to calcthreshpts.m
%     If true, dilate local maxima to try to have only one local max 
%     per neighborhood.  (This is always done for threshold option 3, 
%     regardless of  this input variable.)  If nhoodctrs is not empty (i.e.
%     if neighborhood locations are input), force this to be false.
%     Note that this operation is slow (see Aug. 5, 2012 notes)
% nhoodctrs : Array of neighborhood centers, from previous iteration of 
%     fo5_rp.m for example.  Optional.  If input, avoids processing and 
%     finding of local maxima.  N x 2 array -- col. 1 = x values, col 2 = y
%     fo5_rp.m discard NaNs.
% lsqoptions : [optional] options structure for nonlinear least-squares
%         fitting, from previously running 
%         "lsqoptions = optimset('lsqnonlin');"
%         Inputting this speeds up the function by avoiding redundant 
%         calls to optimset.
%         This variable is only used for non-linear Gaussian fitting.
%
% Output
%   objs : array of object properties; one column per object.  
%          Rows (if only finding centers):
%            [x;
%             y;
%             mass;  (i.e. brightness)
%             particleid;
%             frame;
%             trackid;
%             sigma]
%          sigma is the 2D Gaussian width (std)('nonlineargauss'), the square
%          root of the second moment ('radial') or zero ('centroid')
%          frame and trackid are set to zero for all objects by fo5() and must be
%          dealt with by the functions calling fo5().
%          Rows (if finding centers and orientation):
%            [x;
%             y;
%             mass;  (i.e. brightness)
%             particleid;
%             frame;
%             trackid;
%             theta; (ellipse orientation angle, radians)
%             ra; (ellipse semimajor axis)
%             rb; (ellipse semimajor axis)]
%
% Raghuveer Parthasarathy, (original fo_rp begun April 2007)
% Modifications -- 
%   March 24, 2011 (option for nonlinear 2D Gaussian fit)
%   March 31, 2011 (returns the Gaussian width, sigma)
%   June 8, 2011: Allow threshold to set a max no. of particles
%   August 8, 2011: Allowing thresholding to find centers of only 
%      the 'n' brightest objects
%   August 31, 2011: Allow radial symmetry based fitting
%   October 24, 2011: std. dev. based thresholding
%   Feb. 11, 2012: add max. likelihood Gaussian fit
%   May 2, 2012 : move thresholding to a function: calcthreshimg.m
%   June 27, 2012: allow two-element object size array 
%   July 3, 2012: Calls radialcenter_stk.m if there are >10 pts for radial
%     symmetry fitting, for greater speed.
%   August 2012: fo5_rp.m
%   August 23, 2012: ellipse orientation
%   April 2014: size finding alrorithms for phase separated domains (TH)
% Last modified April 8, 2014 (TH)

function [objs] = fo5_rp(img, processopt, processparam, thresh, fitstr, ...
    try1pernhood, nhoodctrs, lsqoptions)


showplots=false;  % for debugging -- plot things.


%% Defaults
if ~exist('processopt', 'var') || isempty(processopt)
    processopt = 'spatialfilter';
end
% Fitting option; default is radial symmetry (RP 2012)
if ~exist('fitstr', 'var') || isempty(fitstr)
    fitstr = 'radial';
end
if ~exist('try1pernhood', 'var') || isempty(try1pernhood)
    try1pernhood = false;
end
if ~exist('nhoodctrs', 'var')
    nhoodctrs = [];
end
if ~isempty(nhoodctrs)
    try1pernhood = false;
end

% Finding centers, or centers and orientation?
if ischar(fitstr)
    % character array -> only one string, so center-finding only
    ctrfitstr = fitstr;
    orientationstr = [];
    sizestr = [];
elseif iscell(fitstr)
    ctrfitstr = char(fitstr(1));
    orientationstr = char(fitstr(2));
    sizestr = char(fitstr(3));
else
    errordlg('fo5_rp.m: invalid "fitstr"!');
end
% 'none' as an orientation string is equivalent to empty array -- don't
% calculate orientations
if strcmpi(orientationstr, 'none')
    orientationstr = [];
end
if strcmpi(sizestr, 'none')
    sizestr = [];
end
        
if ~exist('lsqoptions', 'var') || isempty(lsqoptions)
    % This variable is only used for non-linear Gaussian fitting.
    if strcmpi(ctrfitstr, 'nonlineargauss')
        lsqoptions = optimset('lsqnonlin');
    else
        lsqoptions = [];
    end
end

% Size parameters -- examine inputs
switch lower(processopt)
    case {'spatialfilter'}
        if length(processparam)==1
            bpfiltsize = processparam;
            nsize = processparam;
        else
            bpfiltsize = processparam(1);
            nsize = processparam(2);
        end
        if bpfiltsize<1
            processopt = 'none';
        end
    case {'gradientvote'}
        gradobjsize = processparam(1);
        nsize = gradobjsize;
        if length(processparam) < 3
            grthresh = 0;
        else
            grthresh = processparam(3);
        end
        if length(processparam) < 2
            graddir = 0;
        else
            graddir = processparam(2);
        end
    case {'none'}
        nsize = processparam(1);
    otherwise
        errordlg('Error in fo5_rp.m: invalid processing option')
end

% Determine thresholding option -- see header comments 
if thresh >= 1.0
    threshopt = 3;
elseif thresh >= 0.0
    threshopt = 1;
else
    threshopt = 2;
    thresh = -thresh;  % note that negative thresh is the indicator of this option;
                       % flip so it's positive.
end

% make img processable
img = double(img);

if showplots
    h1 = figure('Name', 'Original image');
    imshow(img,[]); title('1 Original image')
end

%% Determine neighborhoods for fine localization (unless input)

% Process image
if isempty(nhoodctrs)
    switch lower(processopt)
        case {'spatialfilter'}
            % Bandpass filter -- use Grier et al bpass.m
            noisesize = 1;  % size of noise, px
            processedimg = bpass(img,noisesize,bpfiltsize);
        case {'gradientvote'}
            % gradient vote -- use gradientvote.m
            processedimg = gradientvote(img, gradobjsize, graddir, grthresh);
        case {'none'}
            % do nothing
            processedimg = img;
    end
    if showplots
        h2 = figure('Name', 'Processed image');
        imshow(processedimg,[]); title('2 Processed image')
    end
    % Three options for thresholding -- see above
    % For the chosen option, determine the points that "pass" the threshold.
    % Move to a separate function, so it can be called by the GUI
    [y, x] = calcthreshpts(processedimg, threshopt, thresh, nsize, try1pernhood);
else
    x = nhoodctrs(:,1);
    y = nhoodctrs(:,2);
end

% Get rid of maxima too close to the edge
lenx = 2*floor(nsize/2) + 1;  % 'floor' isn't really necessary, but this
% is the size of "nhood = getnhood(ste);" for a disk structuring
% element of that size.  Note that lenx is forced to be odd.
leny = lenx;  % in principle, could make different
edgeind = ((x < lenx/2) | (x > (size(img,2) - lenx/2))) | ...
    ((y < leny/2) | (y > (size(img,1) - leny/2)));
x(edgeind) = [];
y(edgeind) = [];

% get rid of bad neighborhood centers -- these won't arise from local max
% finding, but could arise if bad nhoodctrs are input
wNan = find(isnan(x) | isnan(y));
x(wNan) = [];
y(wNan) = [];


%% Refine positions

% Compute "masses"
savemass = zeros(1, length(x));
rect = zeros(length(x), 4);
% Compute the first local neighborhood to know the image size, and then all
% the rest.  A bit inelegant; could calculate all rect's at once...
% Skip if there are no objects to find
if ~isempty(x)
    rect(1,:) = [(round(x(1)) - floor(lenx/2)) (round(y(1)) - floor(leny/2)) (lenx-1) (leny-1)];
    cropimg1 = imcrop(img, rect(1,:));
    % all the other neighborhoods
    cropimg = repmat(cropimg1, [1 1 length(x)]); % to allocate memory
    for k = 2:length(x)
        rect(k,:) = [(round(x(k)) - floor(lenx/2)) (round(y(k)) - floor(leny/2)) (lenx-1) (leny-1)];
        cropimg(:,:,k) = imcrop(img, rect(k,:));
    end
end

% Calculate "mass" (intensity) in each neighborhood
nhood = getnhood(strel('disk', floor(nsize/2),0));  % somewhat silly; 
        % could avoid call to image processing toolbox if needed
for k = 1:length(x)
    tempreg = cropimg(:,:,k);
    cropreg = tempreg(nhood);
    savemass(k) = sum(cropreg(:));
end

% Do refinement (find center) around each local max
% If there are many points, use radialcenter_stk.m rather than radialcenter.m
% for radial symmetry method -- even faster!
xcent = zeros(1,length(x));
ycent = zeros(1,length(x));
sigma = zeros(1,length(x));
lsumx = 1:size(cropimg,2);
lsumy = 1:size(cropimg,1);
Lx = lsumx(end);
Ly = lsumy(end);
switch lower(ctrfitstr)
    case {'radial'}
        % Radial-symmetry based fit -- fast, accurate
        % If <10 points, use radialcenter_stk.m for extra speed (avoid
        % redundant grid calculations); else radialcenter.m
        if length(x) < 10
            for j = 1:length(x)
                [xcent(j), ycent(j), sigma(j)] = radialcenter(cropimg(:,:,j));
            end
        else
            [xcent ycent sigma] = radialcenter_stk(cropimg) ;
        end
        % Is the center within reasonable bounds?
        % If not, replace with centroid
        % frequency of bad cases ~ 1/100,000 !  (Extremely rare)
        % See notes Oct. 26, 2011
        % This conditional statement can slow things; Delete?
        badcase = find(abs(xcent - Lx/2)>1.5*Lx | abs(ycent - Ly/2)>1.5*Ly);
        for j = badcase
            ci = cropimg(:,:,j);
            xcent(j) = sum(sum(ci) .* lsumx) / sum(sum(ci));
            ycent(j) = sum(sum(ci,2) .* lsumy') / sum(sum(ci));
        end
    case {'gaussmle'}
        % Gaussian fit via maximum likelihood estimmation -- most accurate
        for j = 1:length(x)
            [A, xcent(j), ycent(j), sigma(j), offset] = gaussfit2DMLE(cropimg(:,:,j));
        end
    case {'nonlineargauss'}
        % Gaussian fit via nonlinear least squares
        for j = 1:length(x)
            [A, xcent(j), ycent(j), sigma(j), offset] = gaussfit2Dnonlin(cropimg(:,:,j), [], [], [], [], lsqoptions);
            % savemass(j) = A*2*pi*sigma(j)*sigma(j);  % Area under a 2D gaussian
            %        Don't use this, since gives large values for weak Gaussians
        end
    case {'lineargauss'}
        % Linear gaussian fit
        % Linear regression -- fit of cropimg intensity to a 2D Gaussian,
        % via polynomial fit of log(intensity),
        % using gaussfit2D.m with a 0.2 threshold
        for j = 1:length(x)
            [A, x0, sigma_x, y0, sigma_y] = gaussfit2D(lsumx, lsumy, cropimg(:,:,j), 0.2);
            if imag(x0)>0.0
                xcent(j) = 0.0;  % return zero -- nonsense
            else
                xcent(j) = x0;
            end
            if imag(y0)>0.0
                ycent(j) = 0.0;  % return zero -- nonsense
            else
                ycent(j) = y0;
            end
            sigma(j) = 0.5*(sigma_x+sigma_y);  % mean Gaussian width
        end
    case {'weightedlineargauss'}
        % Linearized Gaussian fit, weighted for noise (Stephen M.
        % Anthony, Steve Granick -- see Langmuir paper)
        % Need "noiselevel"  std dev of background noise
        % Function gauss2dcirc.m written by Stephen M.
        % Anthony.
        noiselevel = 220;
        disp('hardwiring noise level!!!  -- re-write this')
        for j = 1:length(x)
            [xcent(j),ycent(j),A,sigma(j)] = ...
                gauss2dcirc(cropimg(:,:,j),repmat(lsumx,size(cropimg,1),1),...
                repmat(lsumy',1,size(cropimg,2)),noiselevel);
        end
    case {'centroid'}
        % centroid (center of mass) fit
        % don't subtract background
        % consider all points at once (not looping)
        sumcropimg = sum(sum(cropimg,1),2); % length = length(x)
        xcent = sum(sum(cropimg,1).*repmat(lsumx, [1,1,length(x)]),2) ./ sumcropimg;
        xcent = reshape(xcent, [1 length(x)]);
        ycent = sum(sum(cropimg,2).*repmat(lsumy', [1,1,length(x)]),1) ./ sumcropimg;
        ycent = reshape(ycent, [1 length(x)]);
    otherwise
        errordlg('Unknown method! [fo5_rp.m]');
end
    
% center position relative to image boundary
xn = xcent + rect(:,1)' - 1; % -1 is to correct for matlab indexing
yn = ycent + rect(:,2)' - 1;

% If desired, calculate the ellipsoid orientation
if ~isempty(orientationstr)
    theta = zeros(1,length(x));
    ra = zeros(1,length(x));
    rb = zeros(1,length(x));
    switch lower(orientationstr)
        % Note that 'none' has already been turned into an empty
        % orientation string
        case {'momentcalc'}
            % simple moment calculation
            % Redundant grid calculations, etc., that could be sped up...
            for j = 1:length(x)
                [~, theta(j) r] = simpleellipsefit(cropimg(:,:,j), [xcent(j) ycent(j)], false, true);
                ra(j) = r(1);
                rb(j) = r(2);
            end
        otherwise
            errordlg('Unknown orientation method! [fo5_rp.m]');
    end
end

%equivdiameter and watershed don't work particularly well. Use the
%bilateral filter method -- TH

if ~isempty(sizestr)
    switch lower(sizestr)
        case {'equivdiameter'}
            for j = 1:length(x)
                cropimgscaled(:,:,j) = cropimg(:,:,j)./max(max(cropimg(:,:,j)));
                threshold(j) = graythresh(cropimgscaled(:,:,j));
                bw(:,:,j) = im2bw(cropimgscaled(:,:,j),threshold(j)); %so, could process more, but
                s = regionprops(bw(:,:,j),'EquivDiameter');
                alldiameters = cat(1, s.EquivDiameter);
                radii(j) = max(alldiameters)/2;
            end
        case {'watershed'}
            %for j = 1:length(x)
            radii = zeros(1,length(x));
            [ L, bgm, fgm4 ] = DomainWatershed( img, [xn; yn] );
            region_IDs = diag(L(round(yn),round(xn)));
            for j = 1:length(x)
                radii(j) = (sum(L(:)==region_IDs(j))/pi)^(1/2);
            end
            %                 Lrgb = label2rgb(L, 'jet', 'w', 'shuffle');
            %                 figure, imshow(Lrgb);
            %                 hold on
            %                 plot(xn,yn,'ro','LineWidth',2);
            %end
        case{'bilateral_filter'}
            sigspace = 3;
            sigintense = 1;
            radii = zeros(1,length(x));
            for j = 1:length(x)
                cropimg(:,:,j) = cropimg(:,:,j)/max(max(cropimg(:,:,j)));
                B = bfilter2(cropimg(:,:,j),5,[sigspace, sigintense]);
                bw = im2bw(B,graythresh(B));
                %                 figure; imshow(bw); %Useful for debugging
                s = regionprops(bw,'area');
                areas = cat(1,s.Area);
                if length(areas)==1;
                    radii(j) = sqrt(areas/pi);
                    xup = 0; %will make this larger to make window larger, if needed
                    xdown = 0;
                    yup = 0;
                    ydown = 0;
         
                    while sum([sum(bw(1,:)), sum(bw(:,1)), sum(bw(end,:)), sum(bw(:,end))])>0 % domain overlaps edge of window
                        if xup<30
                            if yup<30
                                if xdown<30
                                    if ydown<30
                                        if sum(bw(1,:))>0
                                            ydown = ydown+1;
                                        elseif sum(bw(end,:))>0
                                            yup = yup+1;
                                        elseif sum(bw(:,1))>0
                                            xdown = xdown+1;
                                        else
                                            xup = xup+1;
                                        end
                                        rerect = [(round(x(j)) - floor(lenx/2)-xdown) (round(y(j)) - floor(leny/2)-ydown) (lenx+xup+xdown-1) (leny+yup+ydown-1)];
                                        recropimg = imcrop(img, rerect);
                                        recropimg = recropimg/max(max(recropimg));
                                        B = bfilter2(recropimg,5,[sigspace, sigintense]);
                                        %figure; imshow(B); %debug
                                        bw = im2bw(B,graythresh(B));
                                        %figure; imshow(bw); %for debugging
                                        s = regionprops(bw,'area');
                                        areas = cat(1,s.Area);
                                        if length(areas)==1
                                            radii(j) = sqrt(areas/pi);
                                        else
                                            radii(j) = NaN;
                                            bw = 0;
                                        end
                                    else
                                        break
                                    end
                                else
                                    break
                                end
                            else
                                break
                            end
                        else
                            break
                        end
                    end
                else
                    radii(j) = NaN;
                end
            end
            
%             img = img/max(max(img));
%             B = bfilter2(img, 5, [3, .1]);
%             bw = im2bw(B,graythresh(B));
%             s = regionprops(bw,'area','centroid');
%             areas = cat(1,s.Area);
%             centroids = cat(1,s.Centroid);
%             if length(areas)>1
%                 separations = dist(centroids,[xn;yn]);
%                 [mindists, mininds] = min(separations);
%                 ordered_areas = areas(mininds);
%                 radii = (ordered_areas/pi).^(1/2);
%             else
%                 radii = NaN;
%             end
        otherwise
            errordlg('Unknown size method! [fo5_rp.m]');
    end
end
    
if showplots
    h6 = figure('Name', 'Particle centers');
    imshow(zeros(size(img))); title('6 particle centers')
    for j=1:length(xn)
        rectangle('Position', [xn(j)-1 yn(j)-1 2 2], 'Curvature', [1 1], ...
            'Linewidth', 2.0, 'EdgeColor', [1.0 1.0 0.3]);
    end
end

if showplots
    figure; imagesc(img); colormap('gray');
    hold on
    for j=1:length(xn)
        rectangle('Position', [xn(j)-1 yn(j)-1 2 2], 'Curvature', [1 1], ...
            'Linewidth', 2.0, 'EdgeColor', [1.0 1.0 0.3]);
        plot([xn(j)-diameters(j)/2 xn(j)+diameters(j)/2],[yn(j) yn(j)],'y','LineWidth',2);
    end
end

nrows = 7;
objs = zeros(nrows, length(xn));
objs(1,:) = xn;
objs(2,:) = yn;
objs(3,:) = savemass;
objs(4,:) = 1:length(x);
if isempty(orientationstr) && isempty(sizestr)
    objs(7,:) = sigma;
elseif isempty(sizestr)
    % ellipse properties
    objs(7,:) = theta;
    objs(8,:) = ra;
    objs(9,:) = rb;
else
    objs(7,:) = radii;
end