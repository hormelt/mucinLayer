% FRAPevolve2D.m
% A program to analyze two or more FRAP images to determine the molecular 
%    diffusion coefficient.
%
% USES: Image Processing Toolbox
% CALLS:  (imagebinRP.m), subtractplane.m, (fitplane.m), binavg.m
%
% Brief description:
%   Given at least two images, the program returns the each of the 
%   best-fit diffusion coefficients calculated from each pair 
%   {Image1, Image2}, {Image1,Image3}, {Image1, Image4}, ... 
%   {Image1, ImageN} -- i.e. always using the first image as the "initial"
%   image. These are calculated by assuming the image is "blurred" over 
%   time by a simple Fickian diffusion process (i.e. random walks).
%   The program can handle background subtraction, plane subtraction,
%   cropping, and excluded regions.
%
% Function Inputs:
%
% RawImageArray : Two or more (observed) fluorescence images, taken at times
%   "t1," "t2," ... "tN". Input as a 3D array, or leave empty to choose
%   from a dialog box.
%  usetInf : if true, use the last images as "t = Infinity" to calculate an
%   immobile fraction. Default: false. Can also be input from a dialog box.
%   Immobile fraction calculation.  The last (input) fluorescence 
%   image, if at least 3 images are input, can be identified as
%   an image taken at "t -> Infinity" to be used to determine the mobile fraction
%   of fluorophores.  See notes 5 July 2007.
%   * WARNING:*  I am sceptical of the accuracy of this procedure,
%   especially if the initial image isn't close enough to t=0 and the final
%   isn't a good approximation to "t=Infinity"  -- RP 17 August 2007
% backgroundLevel : Background intensity to subtract from each image. Enter
%   a negative number to determine this by selecting a region from the 
%   first image. Can also be input from a dialog box.
%   Can be either a single number, to be subtracted from all images, or 
%   an array of numbers.
%   An accurate background value is important for accurate FRAP
%   measurements! -- (measure e.g. by taking images of a blank slide) 
% tdelay : The time delay between image 1 and each of the other images
%    (except for the t=Infinity image, see below).  E.g. if FRAP images
%    were taken at t = 10, 30, and 65 seconds, the delay should be entered
%    as "20 55". Can also be input from a dialog box.
%
% Dialog Box Inputs (partial list):
%
% The scale of the images (microns/px)
% The background fluorescence intensity to be subtracted -- 
% Various analysis parameters and options.
%
% Options: 
%
% Cropping. The user can crop the images, to consider only the region
%   near the bleached spot.
% Re-scaling the image size.  NOTE: This is probably unnecessary on 
%   contemporary computers (2016).
%   A typical 20X image has a needlessly high pixel
%   density.  Due to the large number of pixels, the analysis will be slow.
%   The program allows the images to be re-scaled, e.g. by a factor of 4
% Region of interest. 'Speckles' and other junk do not diffuse, and so
%   disturb the analysis.  The user can specify a "region of interest", a
%   binary map.  Points not in these regions are not included in the
%   calculation of chi^2.  The region of interest can either be a
%   user-input image file, or (better) can be created within the program 
%   by specifying polygonal regions.  In addition to any user-specified region of
%   interest, the program creates a circular region of interest around the center
%   of the bleached spot, reaching the image boundaries, since a spherical region
%   seems to give slightly more accurate results than a rectangular image array.
%   The bleached spot center is determined as the centroid of the inverted first
%   image.
%  Radial intensity profile.  (Optional, recommended)
%   The program can calculate the 
%   azimuthally-averaged radial intensity profile, based on the centroid 
%   position determined as above.  This takes several seconds to calculate
%   and is not used further in the analysis.  Still, it is illustrative.
% Option: Scaling image intensities.  To deal with overall bleaching, the mean 
%   intensity of the images can be scaled to equal that of the first
%   (Recommended)
% Option: Subtracting a plane.  To deal with non-uniform illumination, can 
%   fit a plane to the last image, and subtract this from all images. 
%   Not recommended, unless the illumination is significantly non-uniform 
%   *and* the FRAP spot is very centered.
% 
% Procedure:  
% The initial (t1) image is 'evolved' for j = 1 to N time
%   steps, where N is the maximum possible given the image scale and the
%   final time (t2).  At each step, the image is compared to the final (t2)
%   image.  The timestep j for which the evolved image best matches the
%   final observed image (minimal chi^2) gives the Diffusion coefficient:
%   D = dx*dx*j/(2.0*(t2-t1)), where dx is the pixel size.
%   (See "Random Walks in Biology" by H. Berg for details, Chapter 1 and
%   Appendix B)  Implementing the diffusion equation (D Del^2(c) = dc/dt)
%   as a 2D difference equation is equivalent to replacing each point,
%   at each time step, with the average of its neighbors.
%   This is done by convolution with a nearest-neighbor matrix, 
%   including first-, second-, and third-nearest neighbors with 
%   appropriate weights.
% Each 'timestep' corresponds to time dx*dx / (2D).  (In 2D or in any  
%   dimensions, this gives the correct macroscopic behavior.
% Noise: To account for non-ideal noise in the image, the program adds
%   noise to the evolved image.  This noise has the same standard deviation
%   as that of a user-selected 'uniform' region.  The region is determined
%   based on the first image, but the noise is determined from the last
%   (presumably most uniform) image.
% Display images are re-scaled to enhance the contrast.  This does not
%   alter the image arrays or any FRAP calculations.
%
% TESTING:  Running FRAPevolve2D on simulated images of Gaussian intensity
%   profiles (ideal diffusion from a point source, simdiff2D.m) shows agreement 
%   within about 10% between the calculated and true diffusion coefficients
%   for images calculated over a wide range of diffusion times
%   (For e.g. D~3 um^2/s, 600x600 px, 0.11 um/px, times t= 5:10:45s).  
%   For both simulated and real data, the best-fit D value has a small, 
%   systematically incorrect dependence on the time between image, decreasing
%   as this time increases.  See the note below, on the chi^2 curve, for further
%   discussion.
% NOTE regarding chi^2(D).  For each time delay examined, the program
%   displays a chi^2 vs D curve, the minimum of which indicates the
%   best-fit D value.  The shape of this curve is important: it will drop
%   to its minimal value, rise, and plateau.  The plateau (d(chi^2)/dD=0)
%   at large D indicates that continued "evolution" of the 
%   initial image does nothing to its
%   agreement with the observed final image -- the evolved image is already
%   flat, and the disagreement is due to the non-ideal properties of either
%   the data or the analysis.  *Analyses with large plateaus give best-fit
%   D values that are lower than the true value.*  I.e. the data appear
%   "slower" than they are, since they never match the calculated image.
%   *Analyses with no plateau give best-fit D values that are higher than
%   the true value,* perhaps because too little image
%   space or temporal space is sampled by the evolution.  The most accurate
%   values are for data that just indicate an incipient plateau at the edge
%   of the chi^2(D) curve.
%
% Outputs
%
% best_smoothed_D : Diffusion coefficient values (micron^2/s) for each
%   image pair (Image_1-Image_j), after smoothing the chi^2 curve
% sigma_D : uncertainty in D for each image pair
%
%
% Raghuveer Parthasarathy
%   25 April 2004 -- first version
%   22 Dec., 2005  Binary region of interest map
%   5 July, 2007 -- Determines mobile fraction, 
%       comparing first and last images (t-> inf).
%   18 August 18, 2007 -- Significant revisions to convolution,
%       image resizing, etc.  Most important: Increasing number of neighbors
%       included in the convolution matrix from 1 to 3 (4 also tested)
%   29 May, 2016 -- rewrite to be a function, rather than a program
% Last modified May 29, 2016 (RP)

function [best_smoothed_D, sigma_D, RawImageArray, backgroundLevel, tdelay, ...
    FileNames, PathName] = ...
    FRAPevolve2D(RawImageArray, usetInf, backgroundLevel, tdelay)

% turn off warning that image may be re-scaled to fit the screen
iptsetpref('ImshowBorder', 'tight');
warning off Images:initSize:adjustingMag

% Randomize (for noise calculation.)
rng('shuffle'); % initialize random numbers


% -------------------------------------------------------------
% 
disp(' ');  disp(' ');  disp(' '); 
disp('********************************************');
disp('****                                    ****');
disp('****          FRAPevolve2D.m            ****');
disp('****                                    ****');
disp('********************************************');
disp('   ');
disp('** Use gray-scale .tif images, 8- or 16-bit');
disp('** Note instructions in command window (this window) as well as the dialog boxes..');
disp(' ');


presentDir = pwd;  

% Load images
if ~exist('RawImageArray', 'var') || isempty(RawImageArray)
    [RawImageArray, Nimages, usetInf, FileNames, PathName] = LoadImages(presentDir);
end

if ~exist('usetInf', 'var') || isempty(usetInf)
    % if the image array was input without specifying this, assume the last
    % image is NOT t=Inf
    usetInf = false;
    Nimages = size(RawImageArray,3);
end
if ~exist('Nimages', 'var') || isempty(Nimages)
    % Number of images to use -- an inelegant way to deal with this
    switch usetInf
        case false
            Nimages = size(RawImageArray,3);
        case true
            Nimages = size(RawImageArray,3)-1;
    end
end


% ----------------------------------------------------------------------
%  Various calculation parameters
prompt = {'Enter image scale (microns/pixel):', ...
          'Enter the time between image #1 and each subsequent image (not Inf.) (s)',...
          'Enter background value(s) to subtract from each image (0.0 for none; -1 to choose from Image 1)', ...
          'Crop images?  (1==yes):', ...
          'Scale mean intensity of all images to match image 1?  (1==yes):', ...
          'Subtract a plane from all images?  (1==yes):', ...
          'Mobile fraction  (Overridden if calculated from last image)', ...
          'Re-scale image pixels?  Enter factor (1==no change):', ...
          'Enter the resolution in D desired for the smoothed chi^2 curve (um^2/s):'};
dlg_title = 'Fluorophore parameters'; num_lines= 1;
if exist('tdelay', 'var') 
    default_tdelay = num2str(tdelay*(1:(Nimages-1)));
else
    default_tdelay = num2str(15*(1:(Nimages-1)));
end
if exist('backgroundLevel', 'var')
    default_bkgLevel = num2str(backgroundLevel);
else
    default_bkgLevel = '-1';
end
def     = {'0.11', default_tdelay, default_bkgLevel, '1', '1', '0', '1.0', '1',  '0.1'};  % default values
answer  = inputdlg(prompt,dlg_title,num_lines,def);
scale = str2double(answer(1));
tdelay = str2num(char(answer(2))); %#ok<ST2NM>
if ~(length(tdelay) == (Nimages-1))
    disp('ERROR!  number of time delays must equal Nimages-1.  Press Control-C');
    pause
end
backgroundLevel = str2double(answer(3));
crop = logical(str2double(answer(4)));
normint = logical(str2double(answer(5)));
subpln = logical(str2double(answer(6)));
fm = str2double(answer(7));
resize = str2double(answer(8));
Dbox = str2double(answer(9));

% ---------------------------------------------------------------------
% Subtract background, and display

if backgroundLevel < 0
    % Choose a region in Image 1 in which to calculate the background value
    h_sub = figure('Name', 'Choose Background (to subtract)');
    imagesc(RawImageArray(:,:,1));
    colormap('gray'); title('Image 1: for background measurement');
    disp(' '); disp('* Determine Bacground *');
    disp('   Select a region of fairly uniform intensity -- will define background here.');
    disp('   Double-click or Right-click ends polygon definition.');
    bkgreg = roipoly;
    close(h_sub)    
    backgroundLevel = sum(sum(bkgreg.*double(RawImageArray(:,:,1))))/sum(sum(bkgreg));
    fs = sprintf('Background value: %.1f', backgroundLevel); disp(fs);
end

if max(backgroundLevel)>0
    Asub = subtractbackground(RawImageArray, backgroundLevel);
    % h_sub_image = figure;
    %for j=1:Nimages+usetInf
    %    subplot(1,Nimages+usetInf,j);
    %    imshow(scaletomax(Asub(:,:,j),class(RawImageArray)));
    %    colormap('gray'); title('Sub');
    %end
else
    Asub = RawImageArray;
end


% --------------------------------------------------------------
% Calc min. and max. intensity value for each image (for display)
h_allimages = figure('Name', 'All Images');
% display images rescaled so maximum = maximum of range
for j=1:(Nimages+usetInf)
    subplot(1,Nimages+usetInf,j);
    imshow(scaletomax(Asub(:,:,j),class(RawImageArray))); colormap('gray'); title('Sub');
end


% --------------------------------------------------------------
% Crop image, scale intensities for display

% Calc min. and max. intensity value for each image (for display)
maxAsub = zeros(1,Nimages+usetInf);
minAsub = zeros(1,Nimages+usetInf);
for j=1:Nimages+usetInf
    maxAsub(j) = max(max(Asub(:,:,j)));
    minAsub(j) = min(min(Asub(:,:,j)));
end

hcrop = figure('name', 'Cropping');
if (crop==1)
    disp(' '); disp('* Cropping image *');
    disp('   Select the cropping rectangle from image 1 -- will apply to all images.');
    scaleA1sub = double(intmax(class(RawImageArray)))*(Asub(:,:,1) - minAsub(1))./...
        (maxAsub(1) - minAsub(1)); % scale for display
    scaleA1sub = cast(scaleA1sub, 'like', RawImageArray);
    [tempA1,rect] = imcrop(scaleA1sub); % scale for display
    A = zeros(size(tempA1,1), size(tempA1,2), Nimages);
    for j=1:Nimages+usetInf
        A(:,:,j) = imcrop(Asub(:,:,j), floor(rect));
    end
else
    % Don't crop
    A = Asub;
end
maxA = zeros(1,Nimages+usetInf);
minA = zeros(1,Nimages+usetInf);
maxA(1) = max(max(A(:,:,1)));
minA(1) = min(min(A(:,:,1)));

close(h_allimages);  % earlier display of all images
close(hcrop);

% ------------------------------------------------------------
% Normalize intensities, so all images have the same value (optional)

mA1 = mean(mean(A(:,:,1)));
if normint
    for j=2:Nimages+usetInf
        A(:,:,j) = A(:,:,j)*mA1/mean(mean(A(:,:,j)));  % Normalize
    end
end

% ------------------------------------------------------------
% Subtract a plane from each (normalized) image (optional)
%    For plane subtraction, use the fit plane coefficients from the *last*
%    image, since it is presumably the most uniform
% For all, subtract a plane, keeping the mean value constant
planecoeff = [];  % define to avoid errors in .mat file saving
if subpln
    [A(:,:,Nimages+usetInf), planecoeff] =...
        subtractplane(A(:,:,Nimages+usetInf), true);
    for j=Nimages+usetInf-1:-1:1
        A(:,:,j) = subtractplane(A(:,:,j), true, planecoeff);
    end
end

for j=1:Nimages+usetInf
    maxA(j) = max(max(A(:,:,j)));
    minA(j) = min(min(A(:,:,j)));
end

% ----------------------------------------------------------------
% Determine the region (in Image 1) in which to measure noise
%   (Will measure noise after image rescaling)

hnoise = figure('name', 'Noise region');
imagesc(A(:,:,1));
colormap('gray'); title('Image 1: for noise measurement');
disp(' '); disp('* Determine Noise *');
disp('   Select a region of fairly uniform intensity -- will define "noise" here.');
disp('   Double-click or Right-click ends polygon definition.');
noisereg = roipoly;
close(hnoise)


% ---------------------------------------------------------------
% Region of Interest
% Only use a 'region of interest' to compare images (partially optional)
% Force use of a circular ROI, centered, of maximal size, since this seems
%    aid the accuracy of fitting the circularly symmetric FRAP spot.  In
%    addition, user can define an ROI to ignore various regions
% ROI should be a binary TIFF image, prepared separately, or defined by the
%   user by selecting polygons to ignore

% Find the center of the bleached spot by finding the centroid of the
% dimmest 25%
disp(' ');
NA1 = [size(A,1) size(A,2)];
xarray = repmat(1:NA1(2),NA1(1),1);     % array of x-distances, px
yarray = repmat((1:NA1(1))',1,NA1(2));  % y-distance, px
Ainv = maxA(1)-A(:,:,1);  % inverse brightness of image #1
maxAinv = max(Ainv(:));

dlg_title = 'Center'; num_lines= 1;
prompt = {'Calculate bleached spot centroid (1==yes) or use image center (0)?'};
def     = {'1'};  % default values
answer  = inputdlg(prompt,dlg_title,num_lines,def);
centroidcalc = logical(str2double(answer(1)));
if centroidcalc
    cmx = sum(sum(xarray.*Ainv.*(Ainv>(0.75*maxAinv)))) / ...
        sum(sum(Ainv.*(Ainv>(0.75*maxAinv))));
    cmy = sum(sum(yarray.*Ainv.*(Ainv>(0.75*maxAinv)))) / ...
        sum(sum(Ainv.*(Ainv>(0.75*maxAinv))));
else
    cmx = NA1(2)/2.0;
    cmy = NA1(1)/2.0;
end

% Circular ROI
rcircle = min([cmx (NA1(2)-cmx) cmy (NA1(1)-cmy)]); 
   % radius of the largest circle totally contained in the image
dx = xarray - cmx;   % x-distance from center, px
dy = yarray - cmy;  % y-distance from center, px
r = sqrt(dx.*dx + dy.*dy);
ROIcircle = (r < rcircle);  % circular ROI around the center

hroi = figure('name', 'ROI region');
imagesc(A(:,:,1).*ROIcircle);
colormap('gray'); title('Image for ROI');
hold on
plot(cmx, cmy, 'yo');

% Additional ROI selection
prompt = {'Use a binary region-of-interest map (defined here, or from a separate file)?   (1==yes):  ', ...
        'If yes, Region of interest: (1) Define by selecting polygons to ignore; (2) Load binary image that defines ROI', ...
        '   If selecting polygons, number of polygonal regions to IGNORE', ...
        '   If loading image, (1) to choose filename from a dialog box, (2) to type it manually', ...
        '      If typing manually, ROI image filename (including extension)'};
dlg_title = 'Region of interest parameters'; num_lines= 1;
def     = {'1', '1', '1', '1', 'ROIfile.tif'};  % default values
answer  = inputdlg(prompt,dlg_title,num_lines,def);
useroi = logical(str2double(answer(1)));
defROIopt   = uint8(str2double(answer(2)));
numROI      = uint8(str2double(answer(3)));
loadROIopt  = uint8(str2double(answer(4)));
rFileName   = char(answer(5));
if useroi
    switch defROIopt;
        case 1
            disp(' '); disp('* Define Region of Interest *');
            reg = ones(size(A(:,:,1)));
            for j=1:double(numROI),
                % Select a polygon to IGNORE
                figure(hroi);
                fs = sprintf('   Click corners of ignored polygon no. %d; double-click to end polygon definition', j);
                disp(fs);
                partialreg = roipoly;  % should use presently open Fig. 3
                reg = reg.*not(partialreg);  % multiply all selected regions to get their union
            end
        case 2
            if (loadROIopt==1)
                % dialog box for filename
                [rFileName,rPathName] = uigetfile('*.*', 'Image to load...'); 
                roiimage = imread(strcat(rPathName, rFileName));
            else
                roiimage = imread(rFileName);
            end
            fs = sprintf('Image file %s.', rFileName); disp(fs);
            if (size(roiimage,3) > 1)
                roiimage = mean(roiimage,3); % make gray
            end
            if (max(double(roiimage(:))) > min(double(roiimage(:))))
                % making sure that the image isn't all one value
                mr = 0.5*(max(double(roiimage(:))) + min(double(roiimage(:))));
            else
                mr = 0.5*max(double(roiimage(:)));
            end
            reg = (roiimage > mr);  % make ones and zeros; 
        otherwise
            disp('Error!  Bad region of interest paramters; using entire image!');
            msgbox('Error!  Bad region of interest paramters; using entire image!','Error', 'error') ;
            reg = ones(size(A1)); % Region of interest is the entire image
    end
    close(hroi);
else
    reg = ones(size(A(:,:,1))); % Region of interest is the entire image
end

% Combine circular and selected ROI (intersection)
reg = ROIcircle.*reg;
Nreg = sum(reg(:));  % number of points in region of interest; 

% ----------------------------------------------------------------------
% Since we have the radial position information, plot it (optional)
%  averaging over r using binavg.m.  (Slow)

dlg_title = 'Radial calculation'; num_lines= 1;
prompt = {'Calculate & plot radial intensity profile? (1==yes)'};
def     = {'1'};  % default values
answer  = inputdlg(prompt,dlg_title,num_lines,def);
radcalc = logical(str2double(answer(1)));
mxA = [];

if radcalc
    disp(' ');
    disp('Calculating radial-averaged brightness...');
    binr = 0.0:1.0:rcircle;
    A1d = zeros(1+Nimages+usetInf, size(A,1)*size(A,2));
    A1d(1,:) = r(:);
    for k=1:Nimages+usetInf,
        A1d(k+1,:) = reshape(A(:,:,k),1,size(A,1)*size(A,2));
    end
    [mxA, ~, ~, ~, ~] = binavg(A1d, binr);
    hradial = figure; hold on
    for k=1:Nimages+usetInf,
        plot(scale*mxA(1,:), mxA(k+1,:), '-', 'Color', ...
            [0.7 0.7 0.7]*k/(Nimages+usetInf));
    end
    plot(scale*mxA(1,:), mean(mean(A(:,:,1))), 'b:');
    xlabel('Radial position from centroid, \mum');
    ylabel('Intensity, a.u.');
    title('Radial intensity function');
end

% ----------------------------------------------------------------------
% re-sample image (re-scaling it)
% Also re-scale region of interest -- and re-make it as binary

if (resize > 1)
    % First image, to get size
    %tempA1 = imresize(A(:,:,1), 1.0/resize);
    tempA1 = imagebinRP(A(:,:,1), resize);
    newA = zeros(size(tempA1,1), size(tempA1,2), Nimages+usetInf);
    newA(:,:,1) = tempA1;
    for j=2:Nimages+usetInf
        %newA(:,:,j) = imresize(A(:,:,j), 1.0/resize);
        newA(:,:,j) = imagebinRP(A(:,:,j), resize);
    end
    A = newA;
    % Region of interest
    % regtemp = imresize(reg, 1.0/resize);
    regtemp = imagebinRP(reg, resize);
    reg = (regtemp > 0.5);  % re-make binary, in case re-scaling led to values not 0 or 1
    Nreg = sum(sum(reg));  % number of points in region of interest;
    scale = scale*resize;
    
    % for noise calculation:
    % noiseregtemp = imresize(noisereg, 1.0/resize);
    noiseregtemp = imagebinRP(noisereg, resize);
    noisereg = (noiseregtemp > 0.5);   % re-make binary, in case re-scaling led to values not 0 or 1
    
end

% ----------------------------------------------------------------------
% Noise calculation -- region determined and rescaled earlier.
%   Get noise from last (most uniform) image
Nnoisereg = sum(sum(noisereg));  % number of points in region of interest for noise calc.
meannoisereg = sum(sum(noisereg.*A(:,:,Nimages+usetInf)))/Nnoisereg;  
     % mean intensity in region of interest for noise calc.
stdnoisereg = sqrt(sum(sum(noisereg.*(A(:,:,Nimages+usetInf)-meannoisereg).* ...
    (A(:,:,Nimages+usetInf)-meannoisereg)))/(Nnoisereg-1));
disp(' ');
disp('Noise Region');
fs = sprintf('Mean: %.2f', meannoisereg); disp(fs)
fs = sprintf('Standard deviation: %.2f', stdnoisereg); disp(fs)

% ---------------------------------------------------------
% Determine immobile fraction, comparing first and last images (optional)
% Note that the user was asked for the mobile fraction -- use this (fm) if
% not calculating from the images
if usetInf
    meanA1 = sum(sum(reg.*A(:,:,1)))/Nreg;  % mean of first image, in ROI
    d1Inf = reg.*(A(:,:,1) - A(:,:,Nimages+1)) ./ (A(:,:,1) - meanA1);
    d1Inf = d1Inf(d1Inf>0);
    fm = median(d1Inf(:));     % median over ROI -- don't use mean, since skewed
    disp(' ');   fs = sprintf('Mobile fraction: %.2f', fm); disp(fs);
end

if or((fm<0.0),(fm>1.0))
    fs = sprintf('WARNING! \n  Mobile fraction is calculated to be %.2f! \n  (Will set to 0 or 1)', fm);
    warndlg(fs)
    if (fm > 1.0)
        fm = 1.0;
    elseif (fm < 0.0)
        fm = 0.0;
    end
    fs = sprintf('Mobile fraction set to: %.2f', fm); disp(fs);
end

% --------------------------------------------------------------------
% Determine Nmax, the maximum number of steps

dx = scale;  % microns per pixel
tmax = max(tdelay);  % max delay to consider
Nmax = round(2*tmax / (dx*dx/(2*10.0)));  % maximal number of steps possible to run
Dmin = dx*dx/(2.0*tmax);  % minimal D than can be tested
Dmax = dx*dx*Nmax/(2.0*tmax); % maximal D that will be tested
disp('  ');
fs1 = sprintf('* Calculation will take Nmax = %d steps, spanning a range of D %.2f to %.2f um^2/s', ...
    Nmax, Dmin, Dmax); 
fs2 = sprintf('  Enter a new, lower, maximal D range to examine.');
fsDmax = strcat(fs1, fs2);

prompt = {fsDmax};
dlg_title = 'D range and Nmax'; num_lines= 1;
Dmaxstr = sprintf('%.2f', Dmax);
def     = {Dmaxstr};  % default values
answer  = inputdlg(prompt,dlg_title,num_lines,def);
newDmax = str2double(answer(1));

if (newDmax < Dmax)
    Nmax = round(newDmax*2.0*tmax / dx / dx);  % new Nmax
    Dmax = newDmax;
end


% ---------------------------------------------------------------------
% Evolve first image (simulated diffusion) and compare to *each* of the
%    input images (except t=Infinity, if used)

minchi2=repmat(9e99,Nimages-1,1);
c1evolve = A(:,:,1);  % starting image
N = size(c1evolve);
progbar = waitbar(0, 'evolving image...');  % will display progress
chi2 = zeros(Nimages-1,Nmax);  % preallocate memory
Deff = zeros(Nimages-1,Nmax);  % preallocate memory
% in case no good fit is found, return minimal D and original image
bestN = repmat(1,Nimages-1,1);  
bestevolve = repmat(A(:,:,1), [1 1 Nimages-1]);
bestevolvemob = bestevolve;
bestevolvenoise = bestevolve;
% nn = [0 1 0; 1 0 1; 0 1 0]/4.0;  % nearest neighbor matrix
% nn2 = ([0 1 0; 1 0 1; 0 1 0] + [1 0 1; 0 0 0 ; 1 0 1]/2)/6.0;
%    nearest + 2nd-nearest neighbor matrix.  Note 2nd nearest neighbors are
%    sqrt(2) px away, so their contribution to the Diff. Eq. is divided by 2
nn3 = [0 0 0.25 0 0; 0 0.5 1 0.5 0; 0.25 1 0 1 0.25; ...
       0 0.5 1 0.5 0; 0 0 0.25 0 0]/7;  % to third-neighbors
%nn4 = [0 0.2 0.25 0.2 0; 0.2 0.5 1 0.5 0.2; 0.25 1 0 1 0.25; ...
%       0.2 0.5 1 0.5 0.2; 0 0.2 0.25 0.2 0]/8.6; % to fourth-neighbors
%nn4 = [0 0.1 0.25 0.1 0; 0.1 0.5 1 0.5 0.1; 0.25 1 0 1 0.25; ...
%       0.1 0.5 1 0.5 0.1; 0 0.1 0.25 0.1 0]/7.8; % to fourth-neighbors,
%       counter-balancing the larger (8) number of fourth neighbors
%      Neither nn4 seems as accurate as nn3 when run on simulated data.
%  I am unsure of either\92s correctness as a model of the continuum Laplacian. 
halfnn = (size(nn3,1)-1)/2;
for j=1:Nmax,
    c1evolveold = c1evolve;
    % 2D convolution (finite difference method) of the image matrix with
    % the nearest neighbor matrix, to simulate diffusion.
    % The array edges are problematic, especially because convolution must
    % not be allowed to change the mean of the image values.  Use "mirror
    % image" boundaries, so points on the edges are convolved with the
    % usual array values and their reflections for points outside the
    % boundaries.
    mirrorc1 = zeros(N(1)+2*halfnn, N(2)+2*halfnn);
    mirrorc1(halfnn+1:halfnn+N(1), halfnn+1:halfnn+N(2)) = c1evolveold;
    mirrorc1(halfnn+1:halfnn+N(1), 1:halfnn) = fliplr(c1evolveold(:,1:halfnn));
    mirrorc1(halfnn+1:halfnn+N(1), halfnn+N(2)+1:N(2)+2*halfnn) = ...
        fliplr(c1evolveold(:,N(2)-halfnn+1:N(2)));
    mirrorc1(1:halfnn, halfnn+1:halfnn+N(2)) = flipud(c1evolveold(1:halfnn,:));
    mirrorc1(halfnn+N(1)+1:N(1)+2*halfnn, halfnn+1:halfnn+N(2)) = ...
        flipud(c1evolveold(N(1)-halfnn+1:N(1),:));
    mirrorc1(1:halfnn,1:halfnn) = fliplr(flipud(c1evolveold(1:halfnn,1:halfnn)));
    mirrorc1(1:halfnn,halfnn+N(2)+1:N(2)+2*halfnn) = ...
        fliplr(flipud(c1evolveold(1:halfnn,N(2)-halfnn+1:N(2))));
    mirrorc1(halfnn+N(1)+1:N(1)+2*halfnn,1:halfnn) = ...
        fliplr(flipud(c1evolveold(N(1)-halfnn+1:N(1),1:halfnn)));
    mirrorc1(halfnn+N(1)+1:N(1)+2*halfnn,halfnn+N(2)+1:N(2)+2*halfnn) = ...
        fliplr(flipud(c1evolveold(N(1)-halfnn+1:N(1),N(2)-halfnn+1:N(2))));
    % Convolution with the first-third nearest neighbor matrix
    c1evolve = conv2(mirrorc1, nn3, 'valid');
    %  Enforce the invariance of the mean image intensity, though a quick
    %  examination confirms that it isn't necessary
    c1evolve = c1evolve*(mean(c1evolveold(:))/mean(c1evolve(:)));
    c1evmob = fm*c1evolve + (1.0-fm)*A(:,:,1);  % account for immobile fraction
    c1evnoise = c1evmob + stdnoisereg*randn(N(1),N(2));  
        % adds noise of same std as the roi.
        % previously I had calculated all random numbers at the start to
        % save time, but this requires too much memory
    % compare the evolved image, with random noise,
    %     with each subsequent image, ONLY in the region of interest!
    for k=1:Nimages-1
        % Calculate chi^2.  As usual, chi^2 is the mean square deviation
        % between the calc. and measured images.  However, since we're
        % adding random noise of std=stdnoisereg, even perfect agreement
        % will not give chi^2=0, but rather chi^2 = 2*stdnoisereg^2.
        % Subtract this as an "offset" for uncertainty calculations.
        % Doesn't affect the fit value.
        dev = A(:,:,k+1)-c1evnoise;
        dev2 = dev.*dev.*reg;
        chi2(k,j) = sum(dev2(:))/Nreg - 2*stdnoisereg*stdnoisereg;
        Deff(k,j) = dx*dx*j/(2.0*tdelay(k));  % Diffusion coefficient,
        %  if 'time delay' j were the best match to the data of image k
        if (chi2(k,j) < minchi2(k))
            minchi2(k) = chi2(k,j);
            bestN(k) = j;
            bestevolve(:,:,k) = c1evolve;
            bestevolvemob(:,:,k) = c1evmob;
            bestevolvenoise(:,:,k) = c1evnoise;
        end
    end
    if (mod(j,100)==0)
       waitbar(j/Nmax, progbar, 'evolving image...');
    end
end
close(progbar);

% best-fit diffusion coefficient (one for each image pair)
% Will report the smoothed value, below
bestD = dx*dx*bestN./(2.0*tdelay');


% ---------------------------------------------------------------------
% Calculate boxcar smoothed chi^2 and print results
disp('  '); disp('  ');
disp('**  RESULTS:'); disp('  ');
Dsize = dx*dx./(2.0*tdelay);  % D resolution per point in above scan
best_smoothed_D = zeros(1,Nimages-1);
Dbin = 0.0:Dbox:Dmax;  % bins of D for smoothing
smDeff = zeros(Nimages-1,length(Dbin)-1);
smchi2 = zeros(Nimages-1,length(Dbin)-1);
for k=1:Nimages-1,
    if (Dbox > max(Dsize))
        % boxcar averaging
        Darray = zeros(2,size(Deff,2));
        Darray(1,:) = Deff(k,:);
        Darray(2,:) = chi2(k,:);
        [mxD, ~, ~, ~, ~] = binavg(Darray, Dbin);
        smDeff(k,:) = mxD(1,:);
        smchi2(k,:) = mxD(2,:);
    else
        disp('   Resolution in D too poor for desired smoothing.');
        % used un-smoothed values;
        smDeff = Deff;
        smchi2 = chi2;
    end
end
[msmc, ci] = min(smchi2,[],2);  % msmc(k) is the minimal smoothed
% chi^2, which occurs at column ci(k) of row k
for k=1:Nimages-1
    best_smoothed_D(k) = smDeff(k,ci(k));  % After smoothing: D value at minimal chi^2
end


% Estimate the uncertainty in the diffusion coefficient values by finding
% the D at which chi^2 = 2*minchi^2.  If minchi2 is negative, instead use 
% 2*stdnoisereg as the "confidence" D value -- somewhat crude.
Dhigh  = zeros(1,Nimages-1);
Dlow  = zeros(1,Nimages-1);
cols = 1:size(smchi2,2);
for k=1:Nimages-1
    D2 = 2*msmc(k);
    if (D2 < 0.0)
        D2 = 2*stdnoisereg;
    end
    highchi2 = (smchi2(k,:)>D2) & (cols>ci(k));
    lowchi2 = (smchi2(k,:)>D2) & (cols<ci(k));
    if sum(highchi2)>0
        highci = find(highchi2,1,'first');
        Dhigh(k) = smDeff(k,highci);  % D above optimal D at which chi^2=double
    else
        Dhigh(k) = NaN;
    end
    if sum(lowchi2)>0
        lowci  = find(lowchi2==1,1,'last');
        Dlow(k) = smDeff(k,lowci); % D below optimal D at which chi^2=double
    else
        Dlow(k) = NaN;
    end
end
sigma_D = 0.5*(Dhigh-Dlow);  % Estimate of sigma_D

for k=1:Nimages-1
    fs = sprintf('Images 1 & %d: ', k+1);  disp(fs);
    fs = sprintf('   D = %.2f um^2/s.',  best_smoothed_D(k));  disp(fs);
    fs = sprintf('   sigma_D = %.2f (half range of D= %.2f to %.2f um^2/s)', ...
        sigma_D(k), Dlow(k), Dhigh(k)); disp(fs)
    fs = sprintf('   Minimal chi^2 = %.2e.', msmc(k)); disp(fs);
    disp('  ');
end

% -------------------------------------------------------------------------
% Save MAT file (optional)

dlg_title = 'Save MAT file'; num_lines= 1;
prompt = {'Save a MAT file with important variables? (1==yes)',...
          'Output filename (will add .mat)'};
def     = {'1', 'FRAPout'};  % default values
answer  = inputdlg(prompt,dlg_title,num_lines,def);
outmat = logical(str2double(answer(1)));
Moutfiletemp = char(answer(2));

if outmat
    Moutfile = strcat(Moutfiletemp, '.mat');  % output filename
    save(Moutfile, 'A',  'bestevolve',  'bestevolvemob',  ...
        'bestevolvenoise',  'fm', 'minchi2',  'best_smoothed_D', 'msmc',  'mxA', ...
        'planecoeff', 'reg', 'sigma_D', ...  
        'smDeff', 'smchi2',  'stdnoisereg', 'tdelay');
end


% -------------------------------------------------------------------------
% Plots of chi^2 and smoothed chi^2
hchi2 = figure; hold on
for k=1:Nimages-1
    subplot(1,Nimages-1,k); hold on
    plot(Deff(k,:), chi2(k,:), '-', 'LineWidth', 2.0, 'Color', [0.8 0.8 0.8]);
    plot(smDeff(k,:), smchi2(k,:), '-', 'Linewidth', 2.0, 'Color', [0.0 0.3 0.3]);
    plot([min(smDeff(k,:)) max(smDeff(k,:))], 2*msmc(k)*[1 1], 'b:');
    a = axis;
    axis([0 Dmax a(3) a(4)])
    xlabel('D, \mum^2/s','FontWeight', 'bold');
    ylabel('\chi^2','FontWeight', 'bold');
    title(strcat('\chi^2: ', num2str(1), '-', num2str(k+1)),'FontWeight', 'bold');
end

% Plot mean diffusion coefficient for each pair of images, if Nimages>2
if Nimages>2
    hDplot = figure;
    errorbar(tdelay, best_smoothed_D, sigma_D, 'ro');
    a = axis; axis([0 1.1*max(tdelay) a(3) a(4)]);
    xlabel('Time delay, s'); ylabel('Diffusion coefficient, \mum^2/s');
    title('Best-fit (smoothed) diffusion coefficients');
end

% -------------------------------------------------------------------------
% Display images

dispscale = intmax(class(RawImageArray))/max(mean(A(:,:,1))); % value by which to scale image intensities for display
maxA1 = max(max(A(:,:,1)));
hfinalimages = figure; 


for j=1:Nimages-1

        subplot(Nimages-1,3,(j-1)*3+1);
        imshow(scaletomax(A(:,:,1), class(RawImageArray), maxA1));
        colormap('gray'); title('Image 1');
        subplot(Nimages-1,3,(j-1)*3+2);
        imshow(scaletomax(A(:,:,j+1), class(RawImageArray), maxA1));
        colormap('gray'); title(sprintf('Image %d', j+1));
        subplot(Nimages-1,3,(j-1)*3+3);
        imshow(scaletomax(bestevolvenoise(:,:,j), class(RawImageArray), maxA1));
        %imshow(uint8(bestevolvenoise(:,:,j)*dispscale)); 
        colormap('gray'); title('Best Evolved Image 1');
end

if 1==0
for j=1:Nimages-1
    if isa(A, 'uint8')
        subplot(Nimages-1,3,(j-1)*3+1);
        imshow(scaletomax(A(:,:,1), 'uint8', maxA1));
        % uint8(A(:,:,1)*dispscale)); 
        colormap('gray'); title('Image 1');
        subplot(Nimages-1,3,(j-1)*3+2);
        imshow(uint8(A(:,:,j+1)*dispscale)); colormap('gray'); title(sprintf('Image %d', j+1));
        subplot(Nimages-1,3,(j-1)*3+3);
        imshow(uint8(bestevolvenoise(:,:,j)*dispscale)); colormap('gray'); title('Best Evolved Image 1');
    else
        subplot(Nimages-1,3,(j-1)*3+1);
        imshow(uint16(A(:,:,1)*dispscale)); colormap('gray'); title('Image 1');
        subplot(Nimages-1,3,(j-1)*3+2);
        imshow(uint16(A(:,:,j+1)*dispscale)); colormap('gray'); title(sprintf('Image %d', j+1));
        subplot(Nimages-1,3,(j-1)*3+3);
        imshow(uint16(bestevolvenoise(:,:,j)*dispscale)); colormap('gray'); title('Best Evolved Image 1');
    end
end
end


end

function [RawImageArray, Nimages, usetInf, FileNames, PathName] = LoadImages(presentDir)
%   
% Load images into an array

    % Dialog box for file loading options
    % Note that Nimages is the number of images *not counting* the last if the
    % last is to be used for the immobile fraction calculation
    dlg_title = 'Number of images (options)'; num_lines= 1;
    prompt = {'Total number of images to load (2 or 3)', ...
        'Use the last image as t=Infinity, for immobile fraction calc.? (1==yes)'};
    def     = {'2', '0'};  % default values
    answer  = inputdlg(prompt,dlg_title,num_lines,def);
    Nimages = round(str2double(answer(1)));
    usetInf = logical(str2double(answer(2)));
    % The total number of images is "Nimages" if no
    %    "t=infinity" image is used (see below), or "Nimages+1" otherwise.
    if usetInf
        Nimages = Nimages - 1;
        fs = sprintf('WARNING! Immobile fraction calculation \n   is unreliable -- RP ');
        warndlg(fs)
    end
    if (Nimages<2)
        disp('Number of images must be at least 2!  Press Control-C');
        pause;
    end

    % Load the first image
    % dialog box for filename\E5
    [tempFileName,PathName] = uigetfile('*.*', 'Image to load...');
    tempImage1 = imread(strcat(PathName, tempFileName));
    % switch working directory to image directory
    cd(PathName);
    fs = sprintf('Image file no. 1: %s', tempFileName); disp(fs);

    % Load the other images, including the t=Infinity image
    % Store all as layers in a 3D array
    FileNames = repmat({tempFileName},3,1); % allocates memory, and puts the first file name in Cell 1
    RawImageArray = repmat(tempImage1, [1, 1, Nimages+usetInf]);
    % RawImageArray(:,:,1) = tempImage1;
    for j=2:Nimages+usetInf,
        [tempFileName] = uigetfile('*.*', 'Image to load...');
        FileNames(j) = {tempFileName};
        RawImageArray(:,:,j) = imread(strcat(PathName, char(FileNames(j))));
        fs = sprintf('Image file no. %d: %s', j, char(FileNames(j))); disp(fs);
    end
    
    cd(presentDir)
end

function Asub = subtractbackground(A, bkg)
% subtract a background value from all images of the array "A"
Asub = double(A);
if length(bkg)==1
    bkg = repmat(bkg, [1 size(Asub,3)]);
end
for j=1:size(A, 3)
    Asub(:,:,j) = Asub(:,:,j)-bkg(j);
end

end


function scaleA = scaletomax(A, imclass, maxA)
   % scales 2D array "A" so that the max value (or maxA) is the max of the desired class
   % imclass: desired class ('uint8' or 'uint16'); default: use A's class
   % maxA : use this as the max value of A; leave empty for max(A(:))
   if ~exist('imclass', 'var') || isempty(imclass)
       imclass = class(A);
       disp('here')
   end
   if ~exist('maxA', 'var') || isempty(maxA)
       maxA = max(A(:));
   end
   
   scaleA = double(intmax(imclass))*double(A)/double(maxA);
   scaleA = cast(round(scaleA), imclass);
end