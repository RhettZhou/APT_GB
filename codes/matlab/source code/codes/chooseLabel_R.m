%% Choose a label. by Rhett, 09/17/19   modified 05/16/21
%  e-mail: x.zhou@mpie.de;
%% 
function [BWB,BWTr] = chooseLabel_R(BW,subSize,gap,labelName)

fileID = fopen(labelName,'r');
label = fscanf(fileID,'%d\n');
fclose(fileID);

imgSize = size(BW);
num_height_split = floor((imgSize(1)-subSize)/gap);
num_width_split = floor((imgSize(2)-subSize)/gap);
sH = floor((imgSize(1) - num_height_split*gap - subSize)/2);
sW = floor((imgSize(2) - num_width_split*gap - subSize)/2);
labelTr = label;

BWB = zeros(size(BW));
BWTr = zeros(size(BW));

for i=1:num_height_split
   for j=1:num_width_split
      subImgName = (i-1)*num_width_split+j;
      if label(subImgName) > 0
         label(subImgName) = 1;
      end
      if labelTr(subImgName) > 1
         labelTr(subImgName) = 1;
      else
         labelTr(subImgName) = 0;
      end
      BWB(sH+(i-1)*gap+1:sH+(i-1)*gap+subSize,sW+(j-1)*gap+1:sW+(j-1)*gap+subSize)...
          = BWB(sH+(i-1)*gap+1:sH+(i-1)*gap+subSize,sW+(j-1)*gap+1:sW+(j-1)*gap+subSize)...
          +label(subImgName)*ones(subSize,subSize);
      BWTr(sH+(i-1)*gap+1:sH+(i-1)*gap+subSize,sW+(j-1)*gap+1:sW+(j-1)*gap+subSize)...
          = BWTr(sH+(i-1)*gap+1:sH+(i-1)*gap+subSize,sW+(j-1)*gap+1:sW+(j-1)*gap+subSize)...
          +labelTr(subImgName)*ones(subSize,subSize);
    end
end
end