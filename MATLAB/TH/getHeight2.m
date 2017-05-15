function [ height,gof ] = getHeight2( im, noisesz, kernelsz, rescale,xreg,plotopt )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

h = fspecial('gaussian', 3*kernelsz, 3*noisesz); %First pass filter is large
g = fspecial('gaussian',kernelsz,noisesz);
fim = double(imfilter(im,h));
cfim = double(imfilter(im,g));

for j = 1:size(im,1)
    %running average- better to lose resolution?
    sfim(j,:) = smooth(fim(j,:),3*rescale); csfim(j,:) = smooth(cfim(j,:),rescale);
end

thresh = multithresh(sfim,2); %non-random, lower threshold than otsu
bw = imbinarize(sfim,min(thresh));

ft = fittype('a*erf(b*(x-c))+d'); %note that top of the image is at 0

goodcols = (1+xreg):(size(sfim,2)-xreg);

h = waitbar(0,'Fitting Heights');

for j = goodcols
    data = sfim(:,j);
    startheight = find(bw(:,j)>0,1);
    cdata = csfim(max([(startheight-50),1]):min([numel(data),(startheight+50)]),j);
    
    x = -min([50, startheight-1]):min([50, numel(data)-50]); %this indexing is gnar
    
    
    %Doing this herebecause startpoint change
    fo = fitoptions('Method','NonlinearLeastSquares','Lower',[0, 0,-50, -50],...
        'Upper',[1.2*max(cdata),100,50,50],'StartPoint',[std(cdata), 1, 0,min(cdata)]);
    
    [coeff, gof] = fit(x',cdata,ft,fo);

    %c is the only parameter we care about- the others don't have meaning
    %particularly since blurring/averaging
    ci = confint(coeff);
    height(1,j-xreg) = startheight+coeff.c;
    height(2,j-xreg) = ci(1,3);
    height(3,j-xreg) = ci(2,3);
    
    if j == plotopt
    
    figure; plot(cdata)
    hold on
    plot(coeff.a*erf(coeff.b*(x-coeff.c))+coeff.d)
    plot(startheight+coeff.c,coeff.d,'ro')
    
    end
    
    waitbar(j/max(goodcols),h)
    
end

close(h)

if plotopt>0
    
    figure; imagesc(im); colormap('gray');
    hold on
    plot(goodcols,height(1,:),'Color',[.5 0 0],'LineWidth',2);
    
    figure; imagesc(sfim); colormap('gray');
    hold on
    plot(goodcols,height(1,:),'Color',[.5 0 0],'LineWidth',2);
    
    data = csfim(:,1000);
    x = 1:(numel(data)-xreg);
    [coeff, gof] = fit(x',data(x),ft,fo);
    
end


