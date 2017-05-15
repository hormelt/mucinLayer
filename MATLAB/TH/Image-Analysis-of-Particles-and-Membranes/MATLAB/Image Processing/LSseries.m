% LSseries.m
% 
% program to generate intensity profiles along the same line for a series
% of images; averages along pixels parallel to the line.
% Can also be used for a single image.
% Also returns (displays as text) mean and total intensity along the line.
%
% Also returns the first line scan as a plot with 10X points (interpolated)
%
% If images are color (RGB), averages to return only the mean (grayscale)
% intensity.
%
% calls function calc_avgprofile.m.
%
% Raghuveer Parthasarathy
% 11-12 March, 2004
% Sept. 15, 2010 -- use getnumfilelist.m to get an array of file names
% Jan. 18, 2011 (moved line drawing from calc_avgprofile.m)
% July 10, 2011: Allow multipage TIFF input (file info from getnumfilelist)
% last modified Sept. 12, 2011 

clear all
close all

disp('Ideally, images should be grayscale uint8 .tif; ');
disp('   Color will be averaged -> gray.  Other formats also allowed.')

% -----------------------------------------------------------------------
% Load images: Get File Names
[fbase, frmin, frmax, formatstr, FileName1, FileName2, PathName1 ext ismultipage] = ...
    getnumfilelist;
if (frmax-frmin) > 1
    prompt = {'Analyze every "k-th" image; k = '};
    dlg_title = 'Series option'; num_lines= 1;
    def     = {'1'};  % default values
    answer  = inputdlg(prompt,dlg_title,num_lines,def);
    Nskip = str2double(answer(1));
else
    Nskip = 1;  % not really necessary
end


% -----------------------------------------------------------------------

if ismultipage
    A1  = double(imread(strcat(fbase, ext), 1));
else
    A1 = double(imread(FileName1));
end
if (size(A1,3) > 1)
    iscolor = true;
    A1 = mean(A1,3); 
    disp('Color image:  averaging over third (color) dimension!');
else
    iscolor = false;
end

imgfig = 1;  % the number of the figure with the image in it
figure(imgfig);
imagesc(A1);
colormap('gray')
title('first image; scaled intensity');

% Determine line endpoints
disp('Zoom in if desired; Press Enter to continue.')
pause

disp('Click the endpoints of the linescan line (left click to start, right click to end) -- ONE LINE ONLY.  ');

% User-selected line; get line coordinates. Also gets the profile, but we'll ignore this.
[cx, cy, c] = improfile;
% cx and cy are arrays that give the x and y co-ordinates of the line
n = max(size(c));  % the number of points in the line
fs = sprintf('* Range of line (x,y) = (%.1f,%.1f) to (%.1f,%.1f)', cx(1), cy(1), cx(n), cy(n));
disp(fs);

prompt = {'Enter the half-width (+/-) of pixels over which to average.', ...
        'Enter the image scale (microns/px), or -1 to ignore.'};
dlg_title = 'Line properties'; num_lines= 1;
def     = {'3', '0.11'};  % default values
answer  = inputdlg(prompt,dlg_title,num_lines,def);
navg   = round(str2double(answer(1)));
scale  = str2double(answer(2));


% Plot the peripheral and central lines (removed from calc_avgprofile; 
% somewhat redundant code)
figure(imgfig);
thetaperp = atan2(-1.0*(cx(n)-cx(1)), cy(n)-cy(1));  % angle of the perpendicular to the linescan line
if navg > 0
    jarray = -navg:navg:navg;
else
    jarray = 0;
end
for j=jarray,  % peripheral and central lines
    xshift = cx + j*cos(thetaperp);
    yshift = cy + j*sin(thetaperp);
    LSpar = line(xshift, yshift);
    colorarray = [(j==0) 0 1*(j~=0)];
    set(LSpar, 'Color', colorarray);  % make the peripheral lines blue, central line red
end

% convert pixels to distance (or not)
pos = sqrt((cx-cx(1)).*(cx-cx(1)) + (cy-cy(1)).*(cy-cy(1)));
if (scale > 0.0),
    pos = scale*pos;
end

% total number of frames to consider
if isempty(frmin)
    % only one image
    nimages = 1;
else
    frarray = frmin:Nskip:frmax;
    nimages = length(frarray);
end
avgc = zeros(n,nimages);


for j=1:nimages
    if (j==1)  % first frame
        A = A1;
    else
        if ismultipage
            A  = imread(strcat(fbase, ext), frarray(j));
        else
            framestr = sprintf(formatstr, frarray(j));
            FileName = strcat(fbase, framestr, ext);
            A  = imread(FileName);  % image
        end
    end
    if iscolor
        A = mean(A,3);
    end
    avgc(:,j) = calc_avgprofile(A, cx, cy, n, navg);
end


% ----------------------------------

% subtract a background, if desired -- useful since main application is FRAP.
prompt = {'Enter a background intensity to subtract, or 0 for none.'};
dlg_title = 'Background'; num_lines= 1;
def     = {'0.0'};  % default values
answer  = inputdlg(prompt,dlg_title,num_lines,def);
background   = str2double(answer(1));
avgc = avgc - background;

% ----------------------------------

% Note mean and total intensity of line 

meanc = zeros(1, nimages);
sumc = zeros(1, nimages);
for j=1:nimages,
   meanc(j) = mean(avgc(:,j));
   sumc(j) = sum(avgc(:,j));
end

% ----------------------------------

% plot subtracted intensities
figure(2)
hold on
for j=1:nimages,
    plot(pos, avgc(:,j), 'color', [0.8 0.3*(mod(j+1,3)) 0.7*mod(j,2)]);
end
title('Averaged intensities; background subtracted', 'FontWeight', 'bold');
ylabel('Intensity', 'FontWeight', 'bold');
if (scale > 0.0)
    xlabel('Position (microns)', 'FontWeight', 'bold');
else
    xlabel('Position (pixels)', 'FontWeight', 'bold');
end

% figure(3); 
% interppos = min(pos):scale/10:max(pos);
% plot(interppos, interp1(pos, avgc(:,1),interppos, 'spline'));
% title('INTERPOLATED intensity #1; background subtracted', 'FontWeight', 'bold');
% ylabel('Intensity', 'FontWeight', 'bold');
% if (scale > 0.0)
%     xlabel('Position (microns)', 'FontWeight', 'bold');
% else
%     xlabel('Position (pixels)', 'FontWeight', 'bold');
% end
% 

% display mean and sum intensities
fs = sprintf('File \t mean \t sum\n'); disp(fs)
for j=1:nimages
    fs = sprintf('%d \t %.2f \t %.2f\n', j, meanc(j), sumc(j)); disp(fs)
end

% outp = 1;  % output data
prompt = {'Output line intensity data to file? 1==yes.', 'Output Filename'};
dlg_title = 'Output'; num_lines= 1;
def     = {'1', 'lscan.dat'};  % default values
answer  = inputdlg(prompt,dlg_title,num_lines,def);
outp   = logical(str2double(answer(1)));
fnout   = char(answer(2));
if (outp == true)
    fout = fopen(fnout, 'w');
    for j=1:n,
        fprintf(fout, '%.3f\t ', pos(j));
        for k=1:nimages,
            fprintf(fout, '%.2f\t ', avgc(j,k));
        end
        fprintf(fout, '\n');
    end
    fclose(fout);
end

