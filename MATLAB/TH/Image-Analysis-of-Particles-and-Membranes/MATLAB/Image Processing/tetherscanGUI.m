% tetherscanGUI.m
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
% Interactive GUI for gaussian fitting. 
% Displays Gaussian amplitude (Ag), width (sigma), and half width at half max (HWHM)
%
% Calls binavg.m, gaussfit.m, fitline.m, fitplane.m
%
% Raghuveer Parthasarathy
% based on tetherscan.m (begun Sept. 2008)
% June 25, 2009
% Last modified June 25, 2009.  Could use lots of cleaning up!

% clear all
% close all

function  tetherscanGUI

%--------------------------------------------------

disp(' '); disp(' '); disp(' '); disp(' ');
disp('*******************************************************');
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

% ------------------------------------------------------------------------

% Load image

% File loading options
prompt = {'Enter 1 to choose filenames from a dialog box, 0 to type them manually ', ...
    'Enter the image scale, microns/px'};
dlg_title = 'Loading option and scale'; num_lines= 1;
def     = {'1', '0.11'};  % default values
answer  = inputdlg(prompt,dlg_title,num_lines,def);
loadopt = logical(str2double(answer(1)));
scale = str2double(answer(2));
% Load image
if (loadopt)
    [filename,pathname] = uigetfile('*.*', 'First image of the sequence...');
    cd(pathname);
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
A = imcrop(scaleA);

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
% Polygon with which to subtract background intensity; optional plane
% subtraction
%
title('Select polygon for background subtraction (& opt. plane sub.)')
disp(' ');
disp('Select a polygon (rt-click at end) for *Background* intensity');
reg = roipoly;  % a binary mask of the region of interest
Nreg = sum(reg(:));  % number of pixels in region
backint = sum(sum(A.*reg)) / Nreg;  % background value

A = A - backint;  % subtract background

% Subtract plane -- optional whether we'll use this plane subtracted image
% Calls fitplane, and uses a procedure that's essentially the same as
% subtractplane.m (simpler syntax)
sA = size(A);
x = repmat((1:sA(1))',1,sA(2));
y = repmat(1:sA(2), sA(1), 1);
[a1, a2, c] = fitplane(A(reg), x(reg), y(reg));  % fit plane, using only the polygon region
Asubplane = A - a1*x - a2*y - c;


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
fs = sprintf('Gaussian threshold for axis finding: %.2f', threshold); disp(fs);
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

% for later use
rall = r(:);
intall = A(:);
intallsub = Asubplane(:);

% -------------------------------------------------------------
% GUI -- user interface for determining line properties

%  Initialize and hide the GUI as it is being constructed.
fGUI = figure('Name','Tetherscan Fitting', 'Menubar','none', ...
    'Visible','off','Position',[100,100,1000,800], 'Color', [0.8 0.4 0.1]);

% Construct the components.
% Create axes, for the images
buttonheight = 30;
hsig = uicontrol('Style','text',...
    'String','tetherscanGUI.m -- 2009, R. Parthasarathy',...
    'FontWeight', 'bold', 'Position',[875,50,115,30]);
% Initialization
hupdate    = uicontrol('Style','checkbox',...
    'String','Update','FontWeight', 'bold',...
    'Position',[460,765,95,buttonheight],'value', 0, ...
    'Callback',{@update_Callback});
hsubplane    = uicontrol('Style','checkbox',...
    'String','Subtr.Plane','FontWeight', 'bold',...
    'Position',[460,715,95,buttonheight],'value', 0, ...
    'Callback',{@subplane_Callback});
hcl    = uicontrol('Style','checkbox',...
    'String','Clear plot','FontWeight', 'bold',...
    'Position',[460,665,95,buttonheight],'value', 0, ...
    'Callback',{@clearplot_Callback});
hthreshtextinit = uicontrol('Style','text',...
    'String','threshold (0-1)', 'FontWeight', 'bold', ...
    'Position',[585,780,85,15]);
hthreshtext  = uicontrol('Style','edit',...
    'String','thresh', 'Position',[675,780,40,15], ...
    'Callback',{@threshtext_Callback});
hthresh = uicontrol('Style','slider', ...
    'Max', 0.9999, 'Min', 0.01, 'Value', 0.5, ...
    'SliderStep', [0.02 0.1], 'Position', [720,780,140,15], ...
    'Callback',{@thresh_Callback});
hbinsizetextinit = uicontrol('Style','text',...
    'String','binsize (px)', 'FontWeight', 'bold', ...
    'Position',[585,755,85,15]);
hbinsizetext  = uicontrol('Style','edit',...
    'String','binsize', 'Position',[675,755,40,15], ...
    'Callback',{@binsizetext_Callback});
hbinsize = uicontrol('Style','slider', ...
    'Max', 2, 'Min', 0.25, 'Value', 0.5, ...
    'SliderStep', [0.05 0.25], 'Position', [720,755,140,15], ...
    'Callback',{@binsize_Callback});
hexit    = uicontrol('Style','pushbutton',...
    'String','Exit',...
    'Position',[875,765,115,30],...
    'Callback',{@exit_Callback});
% Messages
hmsg    = uicontrol('Style','text',...
    'String','Message: ', 'FontWeight', 'bold', ...
    'HorizontalAlignment', 'left', 'Position',[50,780,400,15]);


% Initialize the GUI.
% Change units to normalized so components resize automatically.
set([fGUI, hupdate, hsubplane, hcl, hsig, hthreshtextinit, hthresh, hthreshtext, ...
    hbinsizetextinit, hbinsizetext, hbinsize, hexit, hmsg], ...
    'Units','normalized');

% Create all variables here
gthresh = 0.25;  % Gaussian threshold
dr = 0.5;  % bin size, pixels
showlive = 0;
subplaneopt = false;
clearplotopt = false;
Ag = 0.0;
sigma = 0.0;

% Move the GUI to the center of the screen.
movegui(fGUI,'center')
% Make the GUI visible.
set(fGUI,'Visible','on');
figure(fGUI);
set(hthreshtext, 'String', sprintf('%.2f', gthresh));
set(hthresh, 'Value', gthresh);
set(hbinsizetext, 'String', sprintf('%.2f', dr));
set(hbinsize, 'Value', dr);
subplot(1,2,1)
imagesc(A);
colormap('gray');
title('Tether image');


% Callback functions

    function subplane_Callback(source,eventdata)
        % If checked, subtract plane
        subplaneopt = get(source,'Value');
    end

    function clearplot_Callback(source,eventdata)
        % If checked, clear the plot window
        clearplotopt = get(source,'Value');
    end

    function update_Callback(source,eventdata)
        % Show image and Gaussian fit to line profile
        showlive = get(source,'Value');
        while showlive==1
            subplot(1,2,1);
            if subplaneopt
                imagesc(Asubplane); colormap('gray'); title('Plane sub.')
                useint = intallsub;  % use the plane-subtracted image
            else
                imagesc(A); colormap('gray');
                useint = intall;
            end

            % ---------------------------------------------------------------------
            % bin data to get average intensity function
            binlb = min(rall):dr:max(rall);
            [mx, stdx, nx, ind, outN] = binavg([rall'; useint'], binlb);
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

            HWHM = sqrt(log(2)*2*sigma*sigma);  % half width at half max
            fmsg = sprintf('   Amplt. = %.2f, sigma = %.2f um, HWHM from sigma = %.2f um', Ag, sigma, HWHM);
            set(hmsg, 'String', fmsg);

            h2 = subplot(1,2,2);
            if clearplotopt
                cla(h2)
                set(hcl, 'Value', false);
            end
            xlabel('Radial position, \mum');
            ylabel('Intensity');
            hold on
            box on
            plot(rbinum, Ag*exp(-(rbinum-r0).*(rbinum-r0)/2/sigma/sigma), '-',...
                'color', [0.75*subplaneopt 0.25 (1-0.25*subplaneopt)]);
            plot(rbinum, intbin, 'ko');

            % Also plot interpolated positions
            interppos = min(rbinum):scale/10:max(rbinum);
            goodbin = ~isnan(rbinum);
            plot(interppos, interp1(rbinum(goodbin), intbin(goodbin), ...
                interppos, 'spline'), 'Color', 0.7*[0.8 1 1]);

            pause(0.2)
        end
    end


    function threshtext_Callback(source,eventdata)
        % set thresh (text entry)
        gthresh = str2double(get(source,'string'));
        set(hthresh, 'Value', gthresh);
    end

    function thresh_Callback(hObject, source,eventdata)
        gthresh = get(hObject,'Value');
        % Also update text entry box:
        set(hthreshtext, 'String', sprintf('%.2f', gthresh));
    end

    function binsizetext_Callback(source,eventdata)
        % set objsize (text entry)
        dr = str2double(get(source,'string'));
        set(hbinsize, 'Value', dr);
    end

    function binsize_Callback(hObject, source,eventdata)
        dr = get(hObject,'Value');
        % Also update text entry box:
        set(hbinsizetext, 'String', sprintf('%.2f', dr));
    end

    function exit_Callback(source,eventdata)
        % Exit
        if get(hupdate, 'Value')==1
            set(hmsg, 'String', 'UNCHECK the update box before exiting!');
            pause(0.6)
        else
            close all;
        end
    end


end