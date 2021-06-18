function numOfImgPerUnit = subImage(conProj,subSize,ifGB,grainInd,pathName,group,imGapSize)
mkdir(sprintf('%s%s%s',pathName,'\processImg\',num2str(group)))
imgSize=size(conProj);                              % Raw image size
num_width_split = (imgSize(1)-subSize)/imGapSize + 1;        % Num of sub-image in a row
num_height_split = (imgSize(2)-subSize)/imGapSize + 1;       % Num of sub-image in a column
numOfImgPerUnit = num_width_split*num_height_split;                   % Total image number
label = [];
for i=1:num_width_split
   for j=1:num_height_split
       subImg=zeros(subSize,subSize);
       GBIndList = [];
       grainIndList = [];
       for ii=1:subSize
         for jj=1:subSize
            subImg(ii,jj) = conProj((i-1)*imGapSize+ii,(j-1)*imGapSize+jj);
            GBIndList = [GBIndList,ifGB((i-1)*imGapSize+ii,(j-1)*imGapSize+jj)];
            grainIndList = [grainIndList,grainInd((i-1)*imGapSize+ii,(j-1)*imGapSize+jj)];
         end
       end
       subImgName = (i-1)*num_height_split+j;
       GBIndCheck = unique(GBIndList);
       grainIndCheck = length(unique(grainIndList(grainIndList>0)));
       
       if ismember(3,GBIndCheck) && grainIndCheck >= 3
           label(subImgName) = 2; %% Triple
       elseif ismember(3,GBIndCheck) && grainIndCheck == 2
           label(subImgName) = 1; %% GB
       elseif ismember(3,GBIndCheck) && grainIndCheck == 1
           label(subImgName) = 0; %% In Grain
       else
         if ismember(1,GBIndCheck) && grainIndCheck >= 2 
            label(subImgName) = 1;
         elseif ismember(1,GBIndCheck) && grainIndCheck == 1
            label(subImgName) = 0;
         else
             if ismember(2,GBIndCheck) && grainIndCheck == 1
               label(subImgName) = 0;
             else
                a = 1
             end
         end
       end
       
       fileOut = [pathName '\processImg\' num2str(group) '\' num2str(subImgName) '.tif'];
       imwrite(subImg, fileOut);
    end
end

fileOut = [pathName '\processImg\' num2str(group) '\label.txt'];
fid=fopen(fileOut,'wt');
for ii = 1:length(label)
    fprintf(fid,'%d\n',label(ii));
end
fclose(fid);

end

