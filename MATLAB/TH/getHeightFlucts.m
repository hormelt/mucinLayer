function [ a0, a1, d2, z2, I, rmse, rsquare, adjrsquare, xind, yind ] = getHeightFlucts( im, samplesize, xind, yind, Rmat, channel )
%getHeightFlucts Still not exactly certain what this is doing!
%INPUTS: im - image stack from loadZSeriesWithChannels

%% Input Catches

if ~exist('samplesize','var') || isempty(samplesize)
    samplesize = 1000;
end

if ~exist('channel','var') || isempty(channel)
    channel = 1;
end

%% Parameters and Fit Options
%
% options = fitoptions('Method','NonlinearLeastSquares',...  %why these numbers?
%     'Lower',[0 0 0 1],...
%     'Upper',[250 250 5 30],...
% %     'StartPoint',[250 0 .1 5 ],...
%     'MaxFunEvals',1E8,...
%     'MaxIter',1E8);

% f = fittype('a1*(erfc(d2*(x-z2*ones(size(x,1),size(x,2)))))+a0*ones(size(x,1),size(x,2))','independent','x','options',options);

% f = fittype('a1*exp(-(d2*(x-z2))^2)+a0','independent','x','options',options);

% initialize outputs

h = zeros(samplesize,size(im,3)); herror = zeros(samplesize,size(im,3));
a0 = zeros(samplesize,size(im,3)); a1 = zeros(samplesize,size(im,3)); % amplitude, offset
d2 = zeros(samplesize,size(im,3)); z2 = zeros(samplesize,size(im,3)); % witdth, height
A0error = zeros(samplesize,size(im,3)); A1error = zeros(samplesize,size(im,3)); %uncertainties in these
d2error = zeros(samplesize,size(im,3)); z2error = zeros(samplesize,size(im,3));
I = zeros(size(im,4),size(im,3),samplesize);
z = 1:size(im,4);
rmse = zeros(samplesize,size(im,3)); rsquare = zeros(samplesize,size(im,3)); adjrsquare = zeros(samplesize,size(im,3));
rotz = zeros(1,size(im,4));

%% Downsampling Data

rng('shuffle');

if ~exist('xind','var') || isempty(xind)
    xind = round(rand(1,samplesize)*(size(im,1)-1))+1; %slightly complex counting avoids getting an index of 0 and ensures we are correctly bounded in the image
end

if ~exist('yind','var') || isempty(yind)
    yind = round(rand(1,samplesize)*(size(im,2)-1))+1;
end

%% Fits and Outputs

for j = 1:samplesize
    tic
    for k = 1:size(im,3)
        zstack = squeeze(im(yind(j),xind(j),k,:,channel));
        [starta0, startz2] = max(zstack);
        %
        options = fitoptions('Method','NonlinearLeastSquares',...  %why these numbers?
            'Lower',[0 0 0 1],...
            'Upper',[250 250 5 30],...
            'StartPoint',[starta0 0 1/std(zstack) startz2],...
            'MaxFunEvals',1E5,...
            'MaxIter',1E5);
        
        f = fittype('a1*exp(-(d2*(x-z2))^2)+a0','independent','x','options',options);
        
        if exist('Rmat','var')
            for jj = 1:length(z)
                rotz(jj) = Rmat(3,:)*[xind(j) yind(j) z(jj)]';
            end
            [c2,gof]=fit(rotz',zstack,f);
        else
            [c2,gof]=fit(z',zstack,f);
        end
        
        a0(j,k)=c2.a0;
        a1(j,k)=c2.a1;
        d2(j,k)=c2.d2;
        z2(j,k)=c2.z2;
        rmse(j,k)= gof.rmse;
        rsquare(j,k) = gof.rsquare;
        adjrsquare(j,k) = gof.adjrsquare;
        
        %         h(j,k)=z2(j,k)-z(2); %why z(2)?
        
        errors=confint(c2);
        
        herror(j,k)=(errors(2,4)-errors(1,4))/2+0.75/2; %what is the hardcoded number?
        A0error(j,k)=(errors(2,1)-errors(1,1))/2;
        A1error(j,k)=(errors(2,2)-errors(1,2))/2;
        d2error(j,k)=(errors(2,3)-errors(1,3))/2;
        z2error(j,k)=(errors(2,4)-errors(1,4))/2;
        
        I(:,k,j) = a1(j,k)*exp(-(d2(j,k)*(z-z2(j,k))).^2)+a0(j,k);
        
    end
    %     I(:,k,j) = a1(j,k)*(1-erf(d2(j,k)*(z-z2(j,k))))+a0(j,k);
    
    toc
end

a0(:,:,1) = a0; a0(:,:,2) = A0error;
a1(:,:,1) = a1; a1(:,:,2) = A1error;
d2(:,:,1) = d2; d2(:,:,2) = d2error;
z2(:,:,1) = z2; z2(:,:,2) = z2error;

