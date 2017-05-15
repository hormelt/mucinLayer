function filtered_image = gaussianbpf(I,d0,d1)
% Guassian Bandpass Filter
%
% modified by TH 5/13/16 from
% Leonardo O. Iheme (leonardo.iheme@cc.emu.edu.tr)
% 24th of March 2011
%
%   I = The input grey scale image
%   d0 = Lower cut off frequency
%   d1 = Higher cut off frequency
%


f = I;
[nx ny] = size(f);
% f = uint8(f);
fftI = fft2(f,2*nx-1,2*ny-1);
fftI = fftshift(fftI);

% subplot(2,1,1)
% imshow(f,[]);
% title('Original Image')

% subplot(2,1,2)
% fftshow(fftI,'log')
% title('Fourier Spectrum of Image')
% Initialize filter.
filter1 = ones(2*nx-1,2*ny-1);
filter2 = ones(2*nx-1,2*ny-1);
filter3 = ones(2*nx-1,2*ny-1);

for i = 1:2*nx-1
    for j =1:2*ny-1
        dist = ((i-(nx+1))^2 + (j-(ny+1))^2)^.5;
        % Use Gaussian filter.
        filter1(i,j) = exp(-dist^2/(2*d1^2));
        filter2(i,j) = exp(-dist^2/(2*d0^2));
        filter3(i,j) = 1.0 - filter2(i,j);
        filter3(i,j) = filter1(i,j).*filter3(i,j);
    end
end
% Update image with passed frequencies
filtered_image = fftI + filter3.*fftI;

% subplot(2,2,3)
% fftshow(filter3,'log')
% title('Frequency Domain Filter Function Image')
filtered_image = ifftshift(filtered_image);
filtered_image = ifft2(filtered_image,2*nx-1,2*ny-1);
filtered_image = real(filtered_image(1:nx,1:ny));
% filtered_image = uint8(filtered_image);

% subplot(2,1,2)
% imshow(filtered_image,[])
% title('Bandpass Filtered Image')

