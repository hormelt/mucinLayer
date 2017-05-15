function [filtim ] = gaussbp( im, Fst1, Fst2 )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

n1 = 2^nextpow2(size(im,1)+1);
n2 = 2^nextpow2(size(im,2)+1);

fftim = zeros(n1,n2,size(im,3),size(im,4));

for j = 1:size(im,3)
    for k = 1:size(im,4)
        fftim(:,:,j,k) = fft2(im(:,:,j,k),n1,n2);
    end
end

shiftedfft = zeros(size(fftim,1),size(fftim,2),size(fftim,3),size(fftim,4));
filteredfft = zeros(size(fftim,1),size(fftim,2),size(fftim,3),size(fftim,4));
filteredShiftedfft = zeros(size(fftim,1),size(fftim,2),size(fftim,3),size(fftim,4));

Fst1 = 2;
Fst2 = 60;
d = fdesign.bandpass('Fst1,Fp1,Fp2,Fst2,Ast1,Ap,Ast2',1/512,5/512,100/512,500/512,60,1,60);
%d = fdesign.bandpass('Fst1,Fp1,Fp2,Fst2,Ast1,Ap,Ast2',Fst1,5,30,Fst2,60,1,60);
d = fdesign.bandpass('Fst1,Fp1,Fp2,Fst2,Ast1,Ap,Ast2',1/512,5/512,50/512,1,256,200,256);
H = design(d);

%For loop here prevents memory issues
for j = 1:size(im,3)
    for k = 1:size(im,4)
        shiftedfft(:,:,j,k) = fftshift(fftim(:,:,j,k));
        filteredfft(:,:,j,k) = filter(H,shiftedfft(:,:,j,k));
        filteredShiftedfft(:,:,j,k) = ifftshift(filteredfft(:,:,j,k));
    end
end


filtimtemp = ifft2(filteredShiftedfft,n1,n2);
filtim = real(filtimtemp(1:size(im,1),1:size(im,2),:,:));

end

