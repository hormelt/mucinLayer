function [ GoodTracks] = getGoodTracks( arr, lastTime, tol, criterion, plotopt, scale)
%Function to remove shitty tracks from the good ones
%   Detailed explanation goes here

count = 0;

switch criterion
    
    case 'slope'
        
        fline = fittype('m*x');
        checks = [];
        for j = 1:numel(arr)
            thisMSDs = arr(j).MSDs;
            thisMSDsLow = thisMSDs(:,thisMSDs(3,:)<lastTime);
            goodMSDs = thisMSDsLow;
            unqtracks = unique(thisMSDsLow(4,:));
            figure; hold on
            tracks = arr(j).Tracks;
            hold on
            for k = unqtracks
                firstPoint = thisMSDsLow(1,thisMSDsLow(4,:)==k);
                firstPoint = firstPoint(1);
                fitCheck = fit(log(1+thisMSDsLow(3,thisMSDsLow(4,:)==k))',log(thisMSDsLow(1,thisMSDsLow(4,:)==k)/firstPoint)',fline,'StartPoint',1/2);
                checks = [checks fitCheck.m];
                loglog(thisMSDsLow(3,thisMSDsLow(4,:)==k),thisMSDsLow(1,thisMSDsLow(4,:)==k),'r.')
                if fitCheck.m>(1/2+tol)
                    goodMSDs(:,goodMSDs(4,:)==k) = [];
                    tracks(:,tracks(6,:)==k) = [];
                else
                    loglog(goodMSDs(3,goodMSDs(4,:)==k),goodMSDs(1,goodMSDs(4,:)==k),'b.')
                end
            end
            GoodTracks(j).MSDs = goodMSDs;
            GoodTracks(j).tracks = tracks;
            set(gca,'XScale','log','YScale','log')
            xlabel 'Time (s)'
            ylabel 'MSD (\mum^2)'
        end
        
    case 'simplestepdist'
        
        bins = 100:100:6000;
        
        for j = 1:numel(arr)
            tempTracks = arr(j).Tracks;
            figure; hold on
            for k = unique(tempTracks(6,:));
                thisTrack = tempTracks(:,tempTracks(6,:)==k);
                
                for k = 1:round(0.5*size(tempTracks,2))
                    
                    %                     dx_temp = thisTrack((
                    
                    [N,X] = hist(dr);
                    [~,maxind] = max(N);
                    loglog(arr(j).MSDs(1,arr(j).MSDs(4,:)==k),'.','Color',[2*X(maxind)/1500 0 (1500-2*X(maxind))/1500])
                    
                end
                
                set(gca,'XScale','log','YScale','log')
                axis([1e0 1e2 0 8000])
                
            end
            
            GoodTracks = 1;
            
        end
        
    case 'closedformdist'
        
        tic
        nu = repmat([-4.99:.01:5],1000,1);
        c = repmat([0.01:.01:10]',1,1000);
        for j = 1:numel(arr)
            tempTracks = arr(j).Tracks;
            %             figure; hold on
            bestnu = zeros(numel(unique(tempTracks(6,:))),100);
            b = zeros(numel(unique(tempTracks(6,:))),100);
            a = zeros(numel(unique(tempTracks(6,:))),100);
            unqtracks = unique(tempTracks(6,:));
            for k = 1:numel(unique(tempTracks(6,:)));
                
                thisTrack = tempTracks(:,tempTracks(6,:)==unqtracks(k));
                
                for t = 1:100
                    
                    dx = thisTrack(1,(1+t):end)-thisTrack(1,1:(end-t));
                    dy = thisTrack(2,(1+t):end)-thisTrack(2,1:(end-t));
                    
                    dr = dx.^2+dy.^2;
                    mu1 = mean(dr); mu2 = mean(dr.^2); mu3 = mean(dr.^3);
                    
                    mse1 = (mu2/mu1^2-bessely(nu+2,c).*bessely(nu,c)./bessely(nu+1,c).^2).^2;
                    mse2 = (mu3/mu1^3-bessely(nu+3,c).*bessely(nu,c).^2./bessely(nu+3,c).^3).^2;
                    totmse = mse1+mse2;
                    [~, BestFitInd] = min(totmse(:));
                    [I,J] = ind2sub(size(totmse),BestFitInd);
                    bestnu(k,t) = nu(I);
                    bestc = c(J);
                    
                    b(k,t) = 2*mu1*bessely(bestnu(k,t)+1,bestc)/bestc/bessely(bestnu(k,t),bestc);
                    a(k,t) = b(k,t)*bestc^2/4;
                    
                    toc
                    
                end
                
            end
            GoodTracks(j).a = a;
            GoodTracks(j).b = b;
            GoodTracks(j).nu = bestnu;
        end
        
        
        
end


