function [ avmsd ] = calc_msd_av( msd )
%Calculates an average msd from several particle tracks

unique_frames = unique(msd(3,:));

for j = length(unique_frames)
   avmsd(1,j) = mean(msd(1,msd(3,:)==unique_frames(j)));
   avmsd(2,j) = mean(msd(2,msd(3,:)==unique_frames(j)));
   avmsd(3,j) = j;
end
    
avmsd(4,:) = ones(1,length(avmsd(1,:)));

end

