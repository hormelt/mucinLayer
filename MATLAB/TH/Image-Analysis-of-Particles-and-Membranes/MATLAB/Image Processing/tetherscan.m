% tetherscan.m
%
% program to determine the coordinates of a tether axis and return the
% intensity as a function of distance from the axis, 
% averaged over tether length.
% Works by using Gaussian line fits to find tether axis, then determining
% the radial distance from the tether axis of all image points, and
% binning.  Lots of simple geometry!
%
% See notes Sept. 16-18, 2008, and instructions below
%
% Allows bypassing of algorithm's determination of tether axis (useful for
% noisy images) and instead using user's "line 2"
%
% Calls gaussfit.m
% Calls fitline.m
% 
% Raghuveer Parthasarathy
% September 17-18, 2008
% Last modified May 12, 2009


clear all
close all

%--------------------------------------------------

disp(' '); disp(' '); disp(' '); disp(' ');
disp('***********************************');
disp(' ');
disp('tetherscan.m');
disp('  (1) Crop image -- drag mouse to select a tether region');
disp('  (2) Select region to determine background intensity ');
disp('      -- select a polygon, rt. click or dlb. click at end');
disp('  (3) Select a line crossing the tether, roughly perpendicular to it');
disp('      and roughly in the center of the region of interest.');
disp('      Left click for endpoint #1; rt. or double click for endpoint #2.');
disp('  (4) Select a line along the tether, spanning the region of interest.');
disp('      Left click for endpoint #1; rt. or double click for endpoint #2.');
disp('***********************************');

% -------------------------------------------------------------

% Load image

% File loading options
prompt = {'Enter 1 to choose filenames from a dialog box, 0 to type them manually '};
dlg_title = 'Loading option'; num_lines= 1;
def     = {'1'};  % default values
answer  = inputdlg(prompt,dlg_title,num_lines,def);
loadopt = logical(str2double(answer(1)));
% Load image
if (loadopt)
    [filename,pathname] = uigetfile('*.*', 'First image of the sequence...');
    cd(pathname);
    %disp('Changing directory')
    %pwd
    %disp(filename);
else
    fn = input('Enter first image filename (assumes current directory) -- Don"t forget the extension!:  ', 's');
    filename = strcat(pwd, '\', fn);
end
A = imread(filename); 

imclass = class(A);  % Determine the class -- uint8 or uint16
fs = sprintf('Image class: %s', imclass); disp(fs);
if ~or(strcmp(imclass, 'uint8'), strcmp(imclass, 'uint16'))
    disp('Image must be 8 or 16 bit!  Press Control-C');
    pause
end


% -------------------------------------------------------------

% Crop

dA = double(A);
scaleA = uint8((dA - min(dA(:)))*255.0 / (max(dA(:)) - min(dA(:))));
disp('Crop the image to the *desired tether anslysis region.*');
[A,rect] = imcrop(scaleA);
  
% -------------------------------------------------------------
% Check format

A = double(A); % double to allow processing.
if (size(A,3) > 1)
    A = mean(A,3); 
    disp('Color image:  averaging over third (color) dimension!');
end

h=figure;
imagesc(A);
colormap('gray')

% -------------------------------------------------------------
% Subtract background intensity
%
title('Select polygon for background subtraction')
disp(' ');
disp('Select a polygon (rt-click at end) for *Background* intensity');
reg = roipoly;  % a binary mask of the region of interest
Nreg = sum(reg(:));  % number of pixels in region
backint = sum(sum(A.*reg)) / Nreg;

A = A - backint;
        
% -------------------------------------------------------------------------

% Select points for line #1, drawn across the tether, to find tether axis

figure(h)
title('Select line #1: across tether')
disp(' ');
disp('Click the endpoints of line #1 crossing the tether.');
disp('    (left click to start, right click to end) -- ONE LINE ONLY.  ');

% User-selected line; get line coordinates. Also gets the profile, but we'll ignore this.
[cx1, cy1, c1] = improfile;
% cx and cy are arrays that give the x and y co-ordinates of the line
n1 = length(c1);  % the number of points in the line
rectangle('Position', [cx1(1)-0.5 cy1(1)-0.5 1 1], 'Facecolor', 'y');
rectangle('Position', [cx1(n1)-0.5 cy1(n1)-0.5 1 1], 'Facecolor', 'y');

% Select points for line #2, drawn along the tether, to set the range over
% which we find the tether axis

title('Select line #2: range along tether')
disp(' ');
disp('Click the endpoints of line #2 along the tether.');
disp('    (left click to start, right click to end) -- ONE LINE ONLY.  ');

% User-selected line; get line coordinates. Also gets the profile, but we'll ignore this.
[cx2, cy2, c2] = improfile;

% cx and cy are arrays that give the x and y co-ordinates of the line
n2 = length(c2);  % the number of points in the line


% ---------------------------------------------------------------------
% Locate tether axis using Gaussian fit of the intensity along 
% lines parallel to line #1 drawn above

% same calculation method as calc_avgprofile.m

thetaperp = atan2(-1.0*(cx1(n1)-cx1(1)), cy1(n1)-cy1(1));  
    % angle of the perpendicular to line #1

% Use line #2 to determine how many parallel lines to make: parallel lines
% spaced roughly by sp pixels
sp = 1.0;
l2 = sqrt((cx2(n2)-cx2(1))*(cx2(n2)-cx2(1)) + ...
    (cy2(n2)-cy2(1))*(cy2(n2)-cy2(1)));   % length of line 2, px.
npm = round(l2/2.0/sp);  % half-length, line2

xa = zeros(1, 2*npm+1);  % initialize array of tether axis x points
ya = zeros(1, 2*npm+1);  % initialize array of tether axis y points
threshold = 0.25;  % threshold for Gaussian fit
fs = sprintf('Gaussian threshold: %.2f', threshold); disp(fs);
for j=-npm:npm,
    xshift = cx1 + sp*j*cos(thetaperp);
    yshift = cy1 + sp*j*sin(thetaperp);
    cpar = improfile(A, xshift, yshift, n1);
    % Gaussian fit
    [Ag, xa(j+npm+1)] = gaussfit(xshift, cpar, threshold, false);
    [Ag, ya(j+npm+1)] = gaussfit(yshift, cpar, threshold, false);
    % Draw line
    if or(abs(j)==npm,j==0)
        figure(h);
        LSpar = line(xshift, yshift);
        colorarray = [j==0 0.6 1.0];
        set(LSpar, 'Color', colorarray);  
    end
end
% Get rid of any bad points (NaN xa or ya)]
isbadxy = or(isnan(xa), isnan(ya));
xa = xa(~isbadxy);
ya = ya(~isbadxy);
% Draw all the determined axis points
for j=1:length(xa),
    rectangle('Position', [xa(j)-0.5 ya(j) 1 1], 'Facecolor', 'g');
end
xamin = min(xa); xamax = max(xa);
yamin = min(ya); yamax = max(ya);
rangex = abs(xamax-xamin);
rangey = abs(yamax-yamin);

% tether axis equation y = mx+b
% Fit a line -- consider "horizontal" and "vertical" separately, since line
% fitting is sensitive to outliers
% Also, due to tsensitivity, do twice to filter out noisy points (> 2 stdev)
if abs(abs(thetaperp)-(pi/2.0)) < (pi/4.0)
    % "vertical" line
    [bb, sigbb, mm, sigmm] = fitline(ya, xa, false);
    dev = xa-bb-mm*ya;
    v = std(dev);
    [bb, sigbb, mm, sigmm] = fitline(ya(abs(dev)<(2*v)), xa(abs(dev)<(2*v)), false);
    m = 1.0/mm;
    b = -1.0*bb/mm;
else
    % horizontal line
    [b, sigb, m, sigm] = fitline(xa, ya, false);
    dev = ya-b-m*xa;
    v = std(dev);
    [b, sigb, m, sigm] = fitline(xa(abs(dev)<(2*v)), ya(abs(dev)<(2*v)), false);
end    

% Bypass determination of tether axis
prompt = {'Enter 1 bypass the finding of the tether axis and instead use line 2 '};
dlg_title = 'Bypass axis finding'; num_lines= 1;
def     = {'0'};  % default values
answer  = inputdlg(prompt,dlg_title,num_lines,def);
bypassopt = logical(str2double(answer(1)));
if bypassopt
    m = (cy2(n2)-cy2(1))/(cx2(n2)-cx2(1));  % slope
    b = 0.5*(cy2(n2)+cy2(1)-m*(cx2(n2)+cx2(1)));  % intercept
end

% Draw the line
figure(h)
LSpar = line(xa, m*xa+b);
axiscolor = [1.0 0.3 0];
set(LSpar, 'Color', axiscolor);  


% ---------------------------------------------------------------------
% determine r (distance to tether axis) for all points in cropped image

s = size(A);
x = repmat(1:s(2),s(1),1);  % all x positions
y = repmat((1:s(1))',1,s(2));  % all y positions
x2 = (x + m*y - m*b)/(1+m*m);  % all tether axis intersect x values
y2 = (m*(x + m*y)+b)/(1+m*m);  % all tether axis intersect y values
r = sqrt((x2-x).*(x2-x)+(y2-y).*(y2-y));  % distance values (abs)
% Positive and negative radial positions
if abs(m)>1.0
    r = r.*sign(x-x2);
else
    r = r.*sign(y-y2);
end
% figure; surf(r); shading interp

% ---------------------------------------------------------------------
% bin data to get average intensity function

% bin and scale options
prompt = {'Bin size (pixels)', 'Enter the image scale, microns/px ', ...
    'Threshold for radial Gaussian fit'};
dlg_title = 'Bin and Scale'; num_lines= 1;
def     = {'0.5', '0.11', '0.5'};  % default values
answer  = inputdlg(prompt,dlg_title,num_lines,def);
dr = str2double(answer(1)); % resolution, pixels
scale = str2double(answer(2));
gthresh = str2double(answer(3));

rall = r(:);
intall = A(:);
binlb = min(rall):dr:max(rall);
[mx, stdx, nx, ind, outN] = binavg([rall'; intall'], binlb);
rbin = mx(1,:);  % binned distances
intbin = mx(2,:);  % binned intensities
sigintbin = stdx(2,:);  % standard deviation of binned intensities


% --------------------------------------------------------
% determine width (Gaussian fit)

% Use only binned points with relative std < 3
rbinum = rbin*scale;
goodrbinum = rbinum(abs(sigintbin./intbin)<3.0);
goodintbin = intbin(abs(sigintbin./intbin)<3.0);

% Gaussian fit
[Ag, r0, sigma] = gaussfit(goodrbinum, goodintbin, gthresh, false);

disp(' ');
disp(' ');
disp('Gaussian fit to radial intensity profile: ');
fs = sprintf('   sigma = %.2f um', sigma); disp(fs);
fs = sprintf('   HWHM from sigma = %.2f um', sqrt(log(2)*2*sigma*sigma)); disp(fs);

figure; 
xlabel('Radial position, \mum');
ylabel('Intensity');
hold on
box on
plot(rbinum, Ag*exp(-(rbinum-r0).*(rbinum-r0)/2/sigma/sigma), '-',...
    'color', [0 0.25 1]);
plot(rbinum, intbin, 'ko');

% Also plot interpolated positions
interppos = min(rbinum):scale/10:max(rbinum);
goodbin = ~isnan(rbinum);
plot(interppos, interp1(rbinum(goodbin), intbin(goodbin), ...
    interppos, 'spline'), 'Color', 0.7*[0.8 1 1]);


% -------------------------------------
% output data: radial position (px), radial position (um), intensity
prompt = {'Output line intensity data to file? 1==yes.', 'Output Filename'};
dlg_title = 'Output'; num_lines= 1;
def     = {'0', 'tetherscan_out.dat'};  % default values
answer  = inputdlg(prompt,dlg_title,num_lines,def);
outp   = logical(str2double(answer(1)));
fnout   = char(answer(2));
if outp
    fout = fopen(fnout, 'w');
    for j=1:length(rbinum),
        fprintf(fout, '%.2f\t %.3f\t %.1f\n', rbin(j), rbinum(j), intbin(j));
    end
    fclose(fout);
end

