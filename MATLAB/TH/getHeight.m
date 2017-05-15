function [ height ] = getHeight( im, noisesz, kernelsz, rescale,xreg,plotopt )
%Function to obtain height of a bright layer by fitting to error function
%
%   INPUTS
%im- image array
%noisesz- scale of noise (pixel or bead);
%kernelsz- scale of gaussian filter kernel;
%rescale- rescaling factor for running average in x-direction;
%xreg- exluded region; since large filters will mess up edges;
%plotopt- array containing frames to display fit to illumination profile
%   over. To suppress plotting output, make 0;
%
%   OUTPUTS
%height- height at each point in x, stored as [estimate; upperbound;
%lowerbound];

h = fspecial('gaussian', kernelsz, noisesz);

fim = double(imfilter(im,h));

for j = 1:size(im,1)
    %running average- better to lose resolution?
    sfim(j,:) = smooth(fim(j,:),rescale);
end

thresh = multithresh(sfim,2); %non-random, lower threshold than otsu
bw = imbinarize(sfim,min(thresh));

ft = fittype('a*erf(b*(x-c))+d'); %note that top of the image is at 0

goodcols = (1+xreg):(size(sfim,2)-xreg);

h = waitbar(0,'Fitting Heights');

for j = goodcols
    data = sfim(:,j);
    startheight = find(bw(:,j)>0,1);
    
    %Doing this herebecause startpoint change
    fo = fitoptions('Method','NonlinearLeastSquares','Lower',[0, 0,0, -50],...
        'Upper',[100,100,600,50],'StartPoint',[std(data), 1, startheight,min(data)]);
    
    %with large filters the edges will get fucked up
    x = 1:(numel(data)-xreg);
    [coeff, gof] = fit(x',data(x),ft,fo);
    
    %c is the only parameter we care about- the others don't have meaning
    %particularly since blurring/averaging
    ci = confint(coeff);
    height(1,j-xreg) = coeff.c; %could also use gof to obtain uncertainty
    height(2,j-xreg) = ci(1,3);
    height(3,j-xreg) = ci(2,3);
    
    if j == plotopt
        figure; plot(data)
        hold on
        plot(x,coeff.a*erf(coeff.b*(x-coeff.c))+coeff.d)
        plot(coeff.c,coeff.d,'ro')
    end
    
    waitbar(j/max(goodcols),h)
    
end

close(h)

%plot overal fit

if sum(plotopt)>0
    
    figure; imagesc(im); colormap('gray');
    hold on
    plot(goodcols,height(1,:),'Color',[.5 0 0],'LineWidth',2);
    
    figure; imagesc(sfim); colormap('gray');
    hold on
    plot(goodcols,height(1,:),'Color',[.5 0 0],'LineWidth',2);
    
end
