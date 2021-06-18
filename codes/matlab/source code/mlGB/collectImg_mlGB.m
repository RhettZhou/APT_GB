%% separate training and testing image group
function collectImg_mlGB(segment,numOfGBnet,numOfNoise,numOfImgPerUnit,numOfGBnetTest,pathName)
if ~exist('pathName','var')
    pathName = uigetdir({},'Select a folder');
end
rng('shuffle')

% segment = 10;
% numOfGBnet = 150;
% numOfNoise = 4;
% numOfImgPerUnit = 25;
% numOfGBnetTest = 10;

%% ------------------------------------------------------------------------
folderName = [pathName '\trainImg\all'];
if( exist(folderName, 'dir') )
   dos_cmd = sprintf( 'rmdir /S /Q "%s"',folderName);
   [~, ~] = system(dos_cmd);
end
mkdir(sprintf('%s%s',pathName,'\trainImg\all'))
numOfImg = numOfGBnet*numOfImgPerUnit*numOfNoise;
imgShuffle = randperm(numOfImg);
for i = 1:numOfImg
   unitID = ceil(i/numOfImgPerUnit);
   subImgID = mod(i,numOfImgPerUnit);
   if subImgID == 1
     fileID = fopen( [pathName '\processImg\' num2str(unitID) '\label.txt'],'r');
     labeli = fscanf(fileID,'%d\n');
     fclose(fileID);
     label(i:i+numOfImgPerUnit-1) = labeli;
   elseif subImgID == 0
     subImgID = 25;  
   end
   copyfile(sprintf('%s%s%s%s%s%s',pathName,'\processImg\', num2str(unitID),'\',num2str(subImgID),'.tif'),sprintf('%s%s%s%s',pathName,'\trainImg\all\', num2str(imgShuffle(i)),'.tif'));
   labelShuffle(imgShuffle(i)) = label(i);
end
fileOut = [pathName '\trainImg\all.txt'];
fid=fopen(fileOut,'wt');
for i = 1:numOfImg
    fprintf(fid,'%d\n',labelShuffle(i));
end
fclose(fid);

%% ------------------------------------------------------------------------
folderName = [pathName '\testImg\all'];
if( exist(folderName, 'dir') )
   dos_cmd = sprintf( 'rmdir /S /Q "%s"',folderName);
   [~, ~] = system(dos_cmd);
end
mkdir(sprintf('%s%s',pathName,'\testImg\all'))
numOfImgTest = numOfGBnetTest*numOfImgPerUnit*numOfNoise;
imgShuffleTest = randperm(numOfImgTest);
for i = numOfImg+1:numOfImg+numOfImgTest
   unitID = ceil(i/numOfImgPerUnit);
   subImgID = mod(i,numOfImgPerUnit);
   if subImgID == 1
     fileID = fopen( [pathName '\processImg\' num2str(unitID) '\label.txt'],'r');
     labeli = fscanf(fileID,'%d\n');
     fclose(fileID);
     labelT(i-numOfImg:i-numOfImg+numOfImgPerUnit-1) = labeli;
   elseif subImgID == 0
     subImgID = 25;  
   end
   copyfile(sprintf('%s%s%s%s%s%s',pathName,'\processImg\', num2str(unitID),'\',num2str(subImgID),'.tif'),sprintf('%s%s%s%s',pathName,'\testImg\all\', num2str(imgShuffleTest(i-numOfImg)),'.tif'));
   labelShuffleTest(imgShuffleTest(i-numOfImg)) = labelT(i-numOfImg);
end
fileOut = [pathName '\testImg\all.txt'];
fid=fopen(fileOut,'wt');
for i = 1:numOfImgTest
    fprintf(fid,'%d\n',labelShuffleTest(i));
end
fclose(fid);


%% ------------------------------------------------------------------------
for i = segment:segment:numOfGBnet
   writeGroupImg(pathName,i,numOfNoise,numOfImgPerUnit,labelShuffle,labelShuffleTest,numOfGBnetTest)
end


%% ------------------------------------------------------------------------
function writeGroupImg(pathName,i,numOfNoise,numOfImgPerUnit,labelShuffle,labelShuffleTest,numOfGBnetTest)
maxTrainImg = i*numOfNoise*numOfImgPerUnit;

folderNameTemp = [pathName '\S' num2str(maxTrainImg)];                      % Remove previous existing training data
if( exist(folderNameTemp, 'dir') )
   dos_cmd = sprintf( 'rmdir /S /Q "%s"',folderNameTemp);
   [~, ~] = system(dos_cmd);
end

mkdir(sprintf('%s%s%s%s',pathName,'\S',num2str(maxTrainImg),'\train'))
mkdir(sprintf('%s%s%s%s',pathName,'\S',num2str(maxTrainImg),'\test'))

fileOutTrain = [pathName '\S' num2str(maxTrainImg) '\train.txt'];
fileOutTest = [pathName '\S' num2str(maxTrainImg) '\test.txt'];
fidTrain=fopen(fileOutTrain,'wt');
fidTest=fopen(fileOutTest,'wt');

for i = 1:maxTrainImg
      copyfile(sprintf('%s%s%s%s',pathName,'\trainImg\all\', num2str(i),'.tif'),sprintf('%s%s%s%s%s%s',pathName,'\S',num2str(maxTrainImg),'\train\', num2str(i),'.tif'));
      fprintf(fidTrain,'%d\n',labelShuffle(i));
end

for i = 1:numOfGBnetTest*numOfNoise*numOfImgPerUnit
      copyfile(sprintf('%s%s%s%s',pathName,'\testImg\all\', num2str(i),'.tif'),sprintf('%s%s%s%s%s%s',pathName,'\S',num2str(maxTrainImg),'\test\', num2str(i),'.tif'));
      fprintf(fidTest,'%d\n',labelShuffleTest(i)); 
end

fclose(fidTrain);
fclose(fidTest);
end

end