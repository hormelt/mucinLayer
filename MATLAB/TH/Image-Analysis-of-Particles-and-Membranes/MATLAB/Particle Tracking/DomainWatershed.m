function [ L, fgm, bgm ] = DomainWatershed( img, objs )

%%method 1: full technique from mathworks: http://www.mathworks.com/help/images/examples/marker-controlled-watershed-segmentation.html
hy = fspecial('sobel');
hx = hy';
Iy = imfilter(double(img), hy, 'replicate');
Ix = imfilter(double(img), hx, 'replicate');
gradmag = sqrt(Ix.^2 + Iy.^2);
se = strel('disk', 20); %3 works well
Ie = imerode(img, se);
Iobr = imreconstruct(Ie, img);
Iobrd = imdilate(Iobr, se);
Iobrcbr = imreconstruct(imcomplement(Iobrd), imcomplement(Iobr));
Iobrcbr = imcomplement(Iobrcbr);
fgm = imregionalmax(Iobrcbr);
se2 = strel(ones(5, 5)); %was 5
fgm2 = imclose(fgm, se2);
fgm3 = imerode(fgm2, se2);
fgm4 = bwareaopen(fgm3, 20);
Iobrcbr = Iobrcbr/max(max(Iobrcbr));
%bw = im2bw(Iobrcbr, graythresh(Iobrcbr));
bw = im2bw(Iobrcbr, .9);
D = bwdist(bw);
DL = watershed(D);
bgm = DL == 0;
gradmag2 = imimposemin(gradmag, bgm | fgm4);
L = watershed(gradmag2);

%% plot for debug
I4 = img;
I4(imdilate(L==0, ones(1,1)) | bgm | fgm4) = 1;
figure; imagesc(I4); colormap('gray');

%% method 2: use tracked positions for fgms 
hy = fspecial('sobel');
hx = hy';
Iy = imfilter(double(img), hy, 'replicate');
Ix = imfilter(double(img), hx, 'replicate');
gradmag = sqrt(Ix.^2 + Iy.^2);
gradmag = calc_image_gradient(img,[],2);
objs_pixel = round(objs(1:2,:));
fgm = zeros(size(img,1),size(img,2));
for j = 1:size(objs_pixel,2)
    fgm(objs_pixel(2,j),objs_pixel(1,j)) = 1;
end
%se = strel('disk',3);
%fgm(imdilate(fgm==1,se)) = 1;
D = bwdist(fgm);
DL = watershed(D);
bgm = DL == 0;
gradmag2 = imimposemin(gradmag, bgm | fgm);
L = watershed(gradmag2);

end

