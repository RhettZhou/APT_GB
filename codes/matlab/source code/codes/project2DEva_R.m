%% Find the idea theta value. by Rhett, 09/21/19
%  e-mail: x.zhou@mpie.de;
%% This is a function to find the best orientation for 2D projection and generate 2D concentration image.
%% 
function [BW,conProjR,conProj,X,Y,xNum,yNum,dValue,dCValue,dNValue] = ...
    project2DEva_R(theta,resolu,factor,pAll,pSol,meanAll,stdAll)

rho = deg2rad(0);
theta = deg2rad(theta);

[conProj,xNum,yNum,X,Y] = proj(pAll,pSol,theta,rho,resolu);

conProjR = conProj;        %For raw C map image
conProjR(isnan(conProjR)) = 0;  % Avoid some bad point in the image

BW = conProjR;            % For filtered C map image, also used for IE mapping

conMax = max(conProjR(:));
conMin = min(conProjR(:));
if meanAll + factor*stdAll >= conMax || meanAll + factor*stdAll <= conMin
  BW(BW < meanAll) = 0;     
  BW(BW >= meanAll) = 1;
else
  BW(BW < (meanAll + factor*stdAll)) = 0;     
  BW(BW >= (meanAll + factor*stdAll)) = 1;
end

conProjF = conProj;                                   % For accounting
conProjF(conProjF < (meanAll + factor*stdAll)) = NaN;

conProj(conProj < (meanAll + factor*stdAll)) = 0;     % For filtered C map image, but empty region is NaN
conProj(conProj >= (meanAll + factor*stdAll)) = 1;
                                    
meanAllP = mean(conProjF(:),'omitnan');
dValue = meanAllP*(sum(conProj(:) == 1)/sum(~(isnan(conProj(:)))));
dCValue = meanAllP;
dNValue = (sum(conProj(:) == 1)/sum(~(isnan(conProj(:)))));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function p = normProj(p,theta,rho)

[p(:,1),p(:,2),p(:,3)] = rotZ(p(:,1),p(:,2),p(:,3),theta);
[p(:,1),p(:,2),p(:,3)] = rotY(p(:,1),p(:,2),p(:,3),rho);

end

function [xN,yN,zN] = rotZ(x,y,z,theta)

xN = x*cos(theta) - y*sin(theta);
yN = x*sin(theta) + y*cos(theta);
zN = z;

end

function [xN,yN,zN] = rotY(x,y,z,rho)

xN = x*cos(rho) + z*sin(rho);
yN = y;
zN = z*cos(rho) - x*sin(rho);

end

function [conProj, xNum, yNum, X, Y] = proj(pAll,pSol,theta,rho,resolu)

pAllProj = normProj(pAll,theta,rho);
pSolProj = normProj(pSol,theta,rho);
xMax = ceil(max(pAllProj(:,1))/resolu)*resolu + resolu;
xMin = floor(min(pAllProj(:,1))/resolu)*resolu - resolu;
yMax = ceil(max(pAllProj(:,3))/resolu)*resolu + resolu;
yMin = floor(min(pAllProj(:,3))/resolu)*resolu - resolu;

xNum = (xMax-xMin)/resolu;
yNum = (yMax-yMin)/resolu;

[X,Y] = meshgrid(xMin:resolu:xMax,yMin:resolu:yMax);

countsAll(size(X,1),size(X,2)) = zeros;
countsSol(size(X,1),size(X,2)) = zeros;

pAllNum = length(pAllProj);

for j = 1:pAllNum
  xCou = ceil((pAllProj(j,1) - xMin)/resolu);
  yCou = ceil((pAllProj(j,3) - yMin)/resolu);
  countsAll(yCou,xCou) = countsAll(yCou,xCou)+1;
end
 
pSolNum = length(pSolProj);

for j = 1:pSolNum
  xCou = ceil((pSolProj(j,1) - xMin)/resolu);
  yCou = ceil((pSolProj(j,3) - yMin)/resolu);
  countsSol(yCou,xCou) = countsSol(yCou,xCou)+1;
end

conProj = countsSol./countsAll;
       
end

end