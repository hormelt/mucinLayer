% ImageIntensities.m
% 
% Program to calculate the mean intensity over the same region in a series
% of images.
% Can also be used for a single image.
%
% Also asks for the 'time' (e.g. exposure time) values for the images.
%
% Prints the results, puts them into arrays, and makes a simple intensity vs. time plot
% Arrays:  exp = time values
%          meanintarray = mean intensity in region selected
%          stdintarray  = standard deviation of intensity in region
%          selected

% Raghuveer Parthasarathy
% 7 March, 2004 (begun)
% March 2006 -- some modifications
% April 5, 2011 (very minor changes)
% Last modified April 5, 2011 

 
clear all
close all

disp('  ');
disp('  ');
disp('  ');
disp('* Images will be shown with re-scaled intensities, for display; ');
disp('  this does not affect the intensity measurements. ');
disp('* When the first image is shown, select the region to examine by defining a polygon.');

nimages = 0;
while (nimages < 1) 
    nimages = input('Enter the number of images to examine (>= 1):  '); 
end

% images should be a grayscale .tif; saved as array A of type uint8; BMP also works.
loadopt = input('Enter 1 to choose filenames from a dialog box, 0 to type them manually:  ');

meanintarray = zeros(1,nimages);
stdintarray = zeros(1,nimages);

for j=1:nimages,
    if (loadopt==1)
        % dialog box for filename
        [pFileName,pPathName] = uigetfile('*.*', 'Image to load...'); 
        fs = sprintf('Path Name: %s', pPathName); disp(fs);
        fs = sprintf('File Name: %s', pFileName); disp(fs);
        A1 = imread(strcat(pPathName, pFileName),'tif');
        if (j==1)
            % switch working directory to (first) image directory
            presentDir = pwd;
            cd(pPathName);
        end
    else
        pFileName = input('Enter first image filename (assumes current directory) -- Don"t forget the extension!:  ','s');
        A1 = imread(pFileName, 'tif');
    end
    if (j>1)
        close(j-1);  % close previous figure
    end
    fs = sprintf('Image file %s.', pFileName); disp(fs)
    imgfig = j;  % the number of the figure with the image in it
    figure(imgfig);
    imshow(A1, []);

    if (j==1)
        disp('Select region of interest -- user-specified polygon (return to end)');
        % select region of interest -- user-specified polygon (return to end)
        reg = roipoly;  % a binary mask of the region of interest
        Nreg = sum(reg(:));  % number of pixels in region
    end
    
    % calc mean, std of intensity in region selected in image 1
    dA1 = double(A1);
    meanintarray(j) = mean(dA1(reg));
    stdintarray(j) = std(dA1(reg));
    fs = sprintf('Image %d:  Intensity in region = %.1f +/- %.1f', ...
        j, meanintarray(j), stdintarray(j));
    disp(fs);
end

exp = input('Enter the exposure time for each image, in the form "[1 2 3 4]" or "1:4": '); 

if nimages > 1
    % fit intensity v. exposure time
    fs = strcat('Exposure Times:  ', sprintf(' %.2f  ', exp)); disp(fs)
    fs = strcat('Intensities   :  ', sprintf(' %.2f  ', meanintarray)); disp(fs)
    [A, sigA, B, sigB] = fitline(exp, meanintarray, false);
    fs = sprintf('Intensity = A + B*exposure'); disp(fs);
    fs = sprintf('\tA = %.2f +/- %.2f, B = %.2f +/- %.2f', A, sigA, B, sigB); disp(fs);
    fs = sprintf('  A  \t sigma_A   \t   B  \t sigma_B'); disp(fs);
    fs = sprintf('  %.2f \t %.2f    \t  %.2f \t %.2f', A, sigA, B, sigB);
    disp(fs);
    
    fs = strcat('Mean Intensity over all images:  ', ...
        sprintf(' %.2f +/- %.2f ', mean(meanintarray), std(meanintarray))); disp(fs)
    
    figure; plot(exp, meanintarray, 'ko', 'markerfacecolor', [0.9 0.6 0.2]);
    hold on; plot(exp, A + B*exp, 'b-');
    
end
