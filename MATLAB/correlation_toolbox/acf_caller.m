clear

for frame=1:50
    
    data_array(:,:,frame)=double(imread('sample_data.tif',frame));
    
end

dim_list=[1 2 3];
pad_flag = 0;

auto_corr=acf(data_array,dim_list,pad_flag);

figure
imagesc(squeeze(auto_corr(79,:,:)))
figure
imagesc(auto_corr(:,:,25))
