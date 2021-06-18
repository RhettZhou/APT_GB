%% Collect image stacks from experimental result. by Rhett, 09/17/19, modified 05/15/21
%  e-mail: x.zhou@mpie.de;
%% This is a function to collect images stacks for machine learning prediction.
%% Input parameter: 
% BW is the 2D filtered image at the optimized orientation.
% imageSize is the size of each sub-image
% The program requires to select a folder to restore all generated images.
%% Output parameter (stored in a subfolder “APTdata”); 
% All the generated images are stored in a folder named as “APTdata”
%% |Try|: collectImg_R(BW,imageSize);
%% 
function fileOut = collectImg_R(BW,subSize,gap,fileOut)

if( exist(fileOut, 'file') )
   dos_cmd = sprintf( 'del /F /Q /S "%s"',fileOut);
   [~, ~] = system(dos_cmd);
end

rng('shuffle')

imgSize = size(BW);
num_height_split = floor((imgSize(1)-subSize)/gap);
num_width_split = floor((imgSize(2)-subSize)/gap);
sH = floor((imgSize(1) - num_height_split*gap -subSize)/2);
sW = floor((imgSize(2) - num_width_split*gap -subSize)/2);

digit = max(ceil(log10(num_height_split*num_width_split)),1);
digitN = ['%0' num2str(digit) 'u'];

for i=1:num_height_split
   for j=1:num_width_split
       subImg=zeros(subSize,subSize);
       subImgName = (i-1)*num_width_split+j;
       for ii=1:subSize
         for jj=1:subSize
            subImg(ii,jj) = BW(sH+(i-1)*gap+ii,sW+(j-1)*gap+jj);
         end
       end
       dataName = ['/' num2str(subImgName,digitN)];
       h5create(fileOut,dataName,[subSize subSize]);
       h5write(fileOut, dataName, subImg);
    end
end
end
