% imagebinRP.m
% Averages an NxN image into "bins" of mxm pixels, i.e. rescaling by 1/m
%   If N is not an integer multiple of m, ignore "remainder" pixels at
%   right, bottom edges
% Similar to imresize, but more reliable since avoids "unknown" filtering
% Input: A = input image (can be double, uint8, or uint16)
% Input: m = bin size (forced to be integer)
%
% Raghuveer Parthasarathy
% August 17, 2007

function Anew = imagebinRP(A, m)

m = round(m);  % force to be integer

imclass = class(A);  % Determine the image class
A = double(A);

N = size(A);
Nnew = [floor(N(1)/m) floor(N(2)/m)];
Anew = zeros(Nnew);
for j=1:Nnew(1)
    for k=1:Nnew(2)
        Abox = A((j-1)*m+1:(j*m),(k-1)*m+1:(k*m));  % pixels of A
        Anew(j,k) = mean(Abox(:));
    end
end

if strcmp(imclass, 'uint8')       % 8-bit image (0-255)
    Anew = uint8(Anew);
elseif strcmp(imclass, 'uint16')  % 16-bit image (0-65535)
    Anew = uint16(Anew);
end
