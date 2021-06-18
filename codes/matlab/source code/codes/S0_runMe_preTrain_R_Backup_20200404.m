%% Generate training dataset. by Rhett, 09/17/19
%  e-mail: x.zhou@mpie.de;
%% This is a function to generate training and testing data for machine learning application. It takes about 3 days. However, this step is only required one time. Once the data is generated, they could be used repeatedly. 
%% Input parameter: 
% the program requires to select one folder to restore the generated images.
%% Output parameter (In the subfolder named as “1500”); 
% Training and testing images for machine learning application.
% The corresponding labels for machine learning application.
%% |Try|: S0_runMe_preTrain_R;
%%
function S0_runMe_preTrain_R(numOfGBnetRestar,numOfGBnet,segment,pName)

if ~exist('pName','var')
    pathName = uigetdir({},'Select a folder');
end

mkdir(sprintf('%s%s',pathName,'\processImg'))

addpath(genpath('mlGB'));
% required subfiles: addNoiseVolumeFraction; mlGB; proj22D; subImage; 
% collectImg_mlGB; 
% occasionally used: saveLammps;

tic
%% Constant variables
rng('shuffle');

%numOfGBnetRestar = 1;                                                      % The first GBnet for generating training data 
%numOfGBnet = 15;                                                           % Structure units
%segment = 15;                                                              % How many GBnets to output a "Traing and testing" dataset.

numOfNoise = 4;                                                            % For the same GB strucutre, how many different solute distribution can have.
filterFactor = 1.8;                                                        % determine the threshold to file out atoms, used in function proj22D

GBnetSize = 400;  %%% 400                                                  % Atomic structure setup, unit Angstrom
GBnetHeight = 50; %%% 50                                                   % Height of unit box
GBnNetGrainNum = 3;                                                        % How many grains inside the box, periodic condiction

resolu2D = 4;                                                              % Resolution in 2D images, unit Angstrom
subSizeNm = 80;  %%% 80                                                          % Subimage size in unit Angstrom

a = 3.924;                                                                 % Do not need to change, the result is not dependent on this value
stru = 'FCC';                                                              % Also do not matter

subSize = subSizeNm/resolu2D;                                              % 20 pixel
numOfImgPerUnit = (GBnetSize/subSizeNm)^2;                                 % 25 images

%%
for i = numOfGBnetRestar:numOfGBnet+ceil(numOfGBnet/10)                       % There will be one extra for the training dataset.

GBNetWidth = 4 + round(35*rand);                                           % Generate A random GB width
M = mlGB(GBnetSize,GBnetSize,GBnetHeight,GBnNetGrainNum,GBNetWidth,a,stru,i+1,0);   % Generate Polycrystalline
%saveLammps(M);  % if any atom structure needs to be exported.
toc
  for j = 1:numOfNoise
    BbaseCon = 5 + round(4*rand);                                          % base concentration before adding noise
    GBbaseCon = 12 + round(8*rand);                                        % base GB concentration before adding noise
    cluVol = 1 + round(5*rand);                                            % volume fraction of clusters in the bulk
    cluRad = 6 + round(GBNetWidth*0.2*rand);                               % cluster radius in the bulk
    cluCon = 10 + round(15*rand);                                          % concentration in selected cluster in the bulk
    GBcluVol = 70 + round(30*rand);                                        % volume fraction of clusters in the GBs
    GBcluRad = 8 + round(GBNetWidth*0.4*rand);                             % cluster radius in the GBs
    GBcluCon = 70 + round(30*rand);                                        % concentration in selected cluster in the GBs
    MNoise = addNoiseVolumeFraction(M,BbaseCon,cluVol,cluRad,cluCon,GBbaseCon,GBcluVol,GBcluRad,GBcluCon,a,stru);    % Add noise to the system
    group = (i-1)*numOfNoise + j;
    [conProjMS,ifGB,grainInd] = proj22D(MNoise,resolu2D,pathName,group,j,filterFactor); % Generate 2D images
    numOfImgPerSub = subImage(conProjMS,subSize,ifGB,grainInd,pathName,group);     % Cut and save subImages, generate labels
    toc
  end
end

collectImg_mlGB(segment,numOfGBnet,numOfNoise,numOfImgPerUnit,pathName);

end