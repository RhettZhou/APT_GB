% This function is used to add solute with noise into the system
% Based on clusters information, i.e. number, radius, concentration (Bulk&GB)
% MNew = addNoiseClusterInform(M,8,50,4,20,400,4,50);

function M = addNoiseVolumeFraction(M,BbaseCon,cluVol,cluRad,cluCon,GBbaseCon,GBcluVol,GBcluRad,GBcluCon,lattice,stru)

rng('shuffle')
aNum = length(M);                                 % total atom number
GBM = M(M(:,6) ~= 2,:);                           % atoms in GBs or triple junctions
BM = M(M(:,6) == 2,:);                           % atoms in Bulk
GBaNum = length(GBM);                             % num of atoms in clusters
GBfraction = 100*GBaNum/aNum                      % Calculate the fraction of GB atoms


%% Calculate cluster numbers based on volume fraction
if strcmp(stru, 'FCC')
  atomVolume = lattice^3/4;
elseif strcmp(stru, 'BCC')
  atomVolume = lattice^3/2;
end

numAtomPClu = (4*pi/3)*cluRad^3/atomVolume;
cluNum = floor((cluVol/100)*aNum/numAtomPClu)
GBnumAtomPClu = (4*pi/3)*GBcluRad^3/atomVolume;
GBcluNum = floor((GBcluVol/100)*GBaNum/GBnumAtomPClu)


%% Add solute into the bulk
M(assignSolute(M(BM(:,1),1),BbaseCon),2) = 2;
% Add solute into the GBs
M(assignSolute(M(GBM(:,1),1),GBbaseCon),2) = 2;

%% Use cluster to add noice to the bulk
cluRadSig = 0.4*cluRad;                           % Define the cluster sigma
cluNumList = randperm(aNum,cluNum);                 % find the atom line of all clusters
cluRadList = normrnd(cluRad,cluRadSig,[1,cluNum]);             % define cluster size
for i = 1:cluNum                                   % Define the center of clusters
     cluster{i}.center.x = M(cluNumList(i),3);
     cluster{i}.center.y = M(cluNumList(i),4);
     cluster{i}.center.z = M(cluNumList(i),5);
end

cluSelectCount = 0;                                % select atoms in the clusters
for i = 1:aNum
   for j = 1:cluNum
     if ((cluster{j}.center.x - M(i,3))^2 + (cluster{j}.center.y - M(i,4))^2 + ...
         (cluster{j}.center.z - M(i,5))^2) <= cluRadList(j)^2
        cluSelectCount = cluSelectCount + 1;
        cluSelect(cluSelectCount) = M(i,1);
        break
     end 
   end
end

M(assignSolute(cluSelect,cluCon),2) = 2;          % Add solute atoms into clusters
perInClu = 100*length(cluSelect)/aNum
solCon = 100*length(M(M(:,2) == 2,1))/aNum

%% Use clusters to add noice to the GB
GBcluRadSig = 0.4*GBcluRad;                           % Define the cluster sigma

GBcluNumList = GBM(randperm(GBaNum,GBcluNum),1);                 % find the atom line of all GB clusters
GBcluRadList = normrnd(GBcluRad,GBcluRadSig,[1,GBcluNum]);       % define GB cluster size

for i = 1:GBcluNum                                   % Define the center of GB clusters
     GBclu{i}.center.x = M(GBcluNumList(i),3);
     GBclu{i}.center.y = M(GBcluNumList(i),4);
     GBclu{i}.center.z = M(GBcluNumList(i),5);
end

GBcluSelectCount = 0;                                % select atoms in the GB clusters

for i = 1:GBaNum
   for j = 1:GBcluNum
     if ((GBclu{j}.center.x - GBM(i,3))^2 + (GBclu{j}.center.y - GBM(i,4))^2 + ...
         (GBclu{j}.center.z - GBM(i,5))^2) <= GBcluRadList(j)^2
        GBcluSelectCount = GBcluSelectCount + 1;
        GBcluSelect(GBcluSelectCount) = GBM(i,1);
        break
     end 
   end
end

GBcluSelectInbulk = assignSolute(GBcluSelect,GBcluCon);
M(GBcluSelectInbulk,2) = 2;              % Add solute atoms into clusters / Overall database
GBMAfterClu = M(M(:,6) == 1,:);          % Add solute atoms into clusters / GB database

GBperInClu = 100*length(GBcluSelect)/GBaNum
GBsolCon = 100*length(GBMAfterClu(GBMAfterClu(:,2)==2,1))/GBaNum

solCon = 100*length(M(M(:,2) == 2,1))/aNum

end

function B = assignSolute(A,percentage)

Anum = length(A);
Bnum = round((percentage*Anum)/100);

B = A(randperm(Anum,Bnum));

end

