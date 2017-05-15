% RegionIntensities.m
% 
% Program to calculate the mean intensity from several (user-defined) regions 
%    of an image.
%
% April 26, 2007:  Allows use of a *second* image for "calibration" --
% adjusting, e.g., for non-uniform fluorescence illumination.  If the
% calibration option is chosen, the image intensities are measured for the
% second image in exactly the same regions as those chosen for the first
% image.  The "raw" values are first reported, and then the intensities for
% the first image regions divided by those of the second image (multiplied
% by the mean value, to keep the intensity range unchanged).  In other
% words, for each region, Intensity = I1 / (I2 / mean(I2)).  Neglect the
% standard deviation of the calibration image intensity values.
%
% June 3, 2008
% Calculates the mean intensity in all the regions and uses this to
% (optionally) subttract from the image the best-fit intensity plane 
% and save the resulting image.  Useful for correcting for
% non-uniform illumination.  Does not weight by the region sizes.
% Calls fitplane.m and subtractplane.m.
% June 18, 2008:  Allows flipping of image
% August 28, 2008: Verifies that user wants to stop selecting regions
% Raghuveer Parthasarathy
% begun 4 April, 2007
% last modified August 28, 2008

clear all
close all

disp('  ');
disp('  ');
disp('* Images will be shown with re-scaled intensities, for display; ');
disp('   this does not affect the intensity measurements. ');
disp('* When the first image is shown, select the region to examine by defining a polygon.');


% image should be a grayscale .tif; 

% ----------------------------
% Load image file

% Options
prompt = {'Enter 1 to choose the filename from a dialog box, 0 to type it manually ', ...
    'If manual, .tif image filename (present directory): ', ...
    'Threshold to rescale display image (1.0 = none)'};
dlg_title = 'Series option'; num_lines= 1;
def     = {'1','image1.tif', '1.0'};  % default values
answer  = inputdlg(prompt,dlg_title,num_lines,def);
loadopt = logical(str2double(answer(1)));
thresh = str2double(answer(3));
pFileName = char(answer(2));
if (loadopt==1)
    % dialog box for filenames
    [pFileName,pPathName] = uigetfile('*.*', '.tif image...'); 
else
    pPathName = pwd;
end
fs = sprintf('Image file %s.', pFileName); disp(fs);
A1 = imread(strcat(pPathName, pFileName),'tif');
scrsz = get(0,'ScreenSize');  % screen size
% Set figure 1 to occupy the upper right quadrant of the screen
h1 = figure('Position', [scrsz(3)/2 scrsz(4)/2 scrsz(3)/2 scrsz(4)/2]);
dA1 = double(A1);
dA1thresh = dA1.*(dA1 < thresh*max(dA1(:)));
imagesc(dA1thresh); title('Image 1');
colormap('gray');

% ----------------------------------------------------------------------
% Load calibration file

% Options
prompt = {'Enter 1 to normalize intensities using a second, "calibration" image.', ...
    'Enter 1 to choose the filename from a dialog box, 0 to type it manually ', ...
    'If manual, .tif image filename (present directory): '};
dlg_title = 'Series option'; num_lines= 1;
def     = {'0','1','image2.tif'};  % default values
answer  = inputdlg(prompt,dlg_title,num_lines,def);
calibopt = logical(str2double(answer(1)));
loadopt2 = logical(str2double(answer(2)));
pFileName2 = char(answer(3));
if calibopt
    if (loadopt2==1)
        % dialog box for filenames
        [pFileName2,pPathName2] = uigetfile('*.*', '.tif image...');
    else
        pPathName2 = pwd;
    end
    fs = sprintf('Calibration image file %s.', pFileName); disp(fs);
    A2 = imread(strcat(pPathName2, pFileName2),'tif');
    h2 = figure;
    dA2 = double(A2);
    imagesc(dA2); title('Calibration image');
    colormap('gray');
end


% ----------------------------
% Image flipping  options
flipopt=1;

while (flipopt > 0)
    prompt = {'Flip image (will repeat):  0=none, 1=fliplr, 2=flipud, 3=rotate90CW, 4=rotate90CCW, 5=rotate180.'};
    dlg_title = 'Flipping option'; num_lines= 1;
    def     = {'0'};  % default values
    answer  = inputdlg(prompt,dlg_title,num_lines,def);
    flipopt = int8(str2double(answer(1)));
    % Flip all -- image and calibration
    switch flipopt
        case 1
            A1 = fliplr(A1);
            dA1 = fliplr(dA1);
            if calibopt
                A2 = fliplr(A2);
                dA2 = fliplr(dA2);
            end
        case 2
            A1 = flipud(A1);
            dA1 = flipud(dA1);
            if calibopt
                A2 = flipud(A2);
                dA2 = flipud(dA2);
            end
        case 3
            A1 = rot90(A1,3);
            dA1 = rot90(dA1,3);
            if calibopt
                A2 = rot90(A2,3);
                dA2 = rot90(dA2,3);
            end
        case 4
            A1 = rot90(A1);
            dA1 = rot90(dA1);
            if calibopt
                A2 = rot90(A2);
                dA2 = rot90(dA2);
            end
        case 5
            A1 = rot90(A1,2);
            dA1 = rot90(dA1,2);
            if calibopt
                A2 = rot90(A2,2);
                dA2 = rot90(dA2,2);
            end
    end

    if (flipopt>0)
        figure(h1)
        dA1thresh = dA1.*(dA1 < thresh*max(dA1(:)));
        imagesc(dA1thresh); title('Image 1');
        colormap('gray');
    end
end

% array of row and column numbers
sA = size(dA1);
row = repmat((1:sA(1))',1,sA(2));
col = repmat(1:sA(2), sA(1), 1);


% ----------------------------
% Regions of interest

disp('Select regions of interest: ');
disp('   -- For each region, define a polygon (return or dbl-click to end)');
disp('   -- Double-click (no polygon) to finish selecting regions');
figure(h1);
Nreg = 1;
j=0;
if calibopt==1
    fs = sprintf('Region \t Intensity \t StdDev \t CalibIntensity'); disp(fs);
else
    fs = sprintf('Region \t Intensity \t StdDev'); disp(fs);
end

while (Nreg > 0)
    % select region of interest -- user-specified polygon (return to end)
    mistake = true;
    while mistake
        reg = roipoly;  % a binary mask of the region of interest
        Nreg = sum(reg(:));  % number of pixels in region
        if Nreg==0
            % zero region -- make sure not accidental
            button = questdlg('Stop selecting regions?','Stop selection');
            mistake = ~strcmp(button,'Yes');  % false if 'Yes' entered
        else
            mistake=false;  % finite region; continue
        end
    end
    if (Nreg > 0)
        j = j+1;
        % mean position of the region selected
        meanrow(j) = sum(sum(row.*reg)) / Nreg;
        meancol(j) = sum(sum(col.*reg)) / Nreg;
        % calc mean, std of intensity in region selected in image 1
        meanint(j) = sum(sum(dA1.*reg)) / Nreg;
        stdint(j) = sqrt(sum(sum((dA1-meanint(j)).*...
            (dA1-meanint(j)).*reg))/(Nreg-1));
        if calibopt==1,
            meanintcalib(j) = sum(sum(dA2.*reg)) / Nreg;
            fs = sprintf('%d \t %.1f  \t %.1f  \t %.1f',...
                j, meanint(j), stdint(j), meanintcalib(j));
            disp(fs);
        else
            fs = sprintf('%d \t %.1f  \t %.1f', j, meanint(j), stdint(j));
            disp(fs);
        end
    end
end
Nsamples = j;
if calibopt==1,
    meanallcalib = mean(meanintcalib(:));
    disp(' ');
    disp('Normalizing by calibration image intensities: ');
    for j=1:Nsamples,
        fs = sprintf('%d \t %.1f  \t %.1f', ...
            j, meanint(j)*meanallcalib/meanintcalib(j), stdint(j));
        disp(fs);
    end
end

% ------------------------------------------------
% Option: Subtract a plane from the region -- not weighted by region size
% Load calibration file

% Options
prompt = {'Enter 1 to subtract a non-weighted plane from the selected regions, and output image.', ...
    '.tif image filename (present directory): '};
dlg_title = 'Plane subtraction option'; num_lines= 1;
def     = {'0','image_sub.tif'};  % default values
answer  = inputdlg(prompt,dlg_title,num_lines,def);
planesubopt = logical(str2double(answer(1)));
pFileNamesub = char(answer(2));
if planesubopt
    [a1, a2, c] = fitplane(meanint, meanrow, meancol);
    [Asub planecoeffout] = subtractplane(A1, true, [a1 a2 c]);
    imwrite(Asub, pFileNamesub, 'tiff');
end


% ------------------------------------------------

% close figures
close(h1);  
if calibopt==1
    close(h2); 
end



