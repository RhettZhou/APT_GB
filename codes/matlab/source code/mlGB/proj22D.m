function [conProjMS,ifGB, grainInd] = proj22D(M,resolu,pathName,group,imgKey,filterFactor)

% project simulated 3D atom map into 2D concentration map

aTheta = 0;
rho = 0;

pAll = M(:,3:5);
pSol = M(M(:,2) == 2,3:5);
GBinfor(:) = M(:,6);                     % Whether it is a GB, 1 is GB, 2 is in bulk, 3 is triple
grainInfor(:) = M(:,7);                     % Grain index

[conProj,conProjMS,ifGB, grainInd,xNum,yNum,X,Y] = proj(pAll,pSol,GBinfor,grainInfor,aTheta,rho,resolu,filterFactor);


showImg(conProj,conProjMS,ifGB, grainInd,xNum,yNum,X,Y, pathName,group,imgKey);


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

function [conProj,conProjMS,ifGB,grainInd,xNum,yNum,X,Y] = proj(pAll,pSol,GBinfor,grainInfor,aTheta,rho,resolu,filterFactor)

theta = deg2rad(aTheta);
pAllProj = normProj(pAll,theta,rho);
pSolProj = normProj(pSol,theta,rho);
xMax = ceil(max(pAllProj(:,1))/resolu)*resolu;
xMin = floor(min(pAllProj(:,1))/resolu)*resolu;
yMax = ceil(max(pAllProj(:,2))/resolu)*resolu;
yMin = floor(min(pAllProj(:,2))/resolu)*resolu;

xNum = (xMax-xMin)/resolu;
yNum = (yMax-yMin)/resolu;

[X,Y] = meshgrid(xMin+resolu:resolu:xMax,yMin+resolu:resolu:yMax);

countsAll = NaN(size(X,1),size(X,2));
ifGB = NaN(size(X,1),size(X,2));
grainInd = NaN(size(X,1),size(X,2));
countsSol = NaN(size(X,1),size(X,2));

pAllNum = length(pAllProj);

for j = 1:pAllNum
  xCou = ceil((pAllProj(j,1) - xMin)/resolu);
  yCou = ceil((pAllProj(j,2) - yMin)/resolu);
  if isnan(countsAll(yCou,xCou))
     countsAll(yCou,xCou) = 0;   
  end
  countsAll(yCou,xCou) = countsAll(yCou,xCou)+1;
  if isnan(ifGB(yCou,xCou)) || ifGB(yCou,xCou) == 2
     ifGB(yCou,xCou) = GBinfor(j);
  end
  if isnan(grainInd(yCou,xCou))
     grainInd(yCou,xCou) = grainInfor(j);
  end 
end
 
pSolNum = length(pSolProj);

for j = 1:pSolNum
  xCou = ceil((pSolProj(j,1) - xMin)/resolu);
  yCou = ceil((pSolProj(j,2) - yMin)/resolu);
  if isnan(countsSol(yCou,xCou))
     countsSol(yCou,xCou) = 0;   
  end
  countsSol(yCou,xCou) = countsSol(yCou,xCou)+1;
end

conProj = countsSol./countsAll;


conProjMS = conProj;
stdAll = std(conProjMS(ifGB==2), 'omitnan');
meanAll = mean(conProjMS(ifGB==2), 'omitnan');
conProjMS(conProjMS < (meanAll + filterFactor*stdAll)) = 0;
conProjMS(conProjMS >= (meanAll + filterFactor*stdAll)) = 1;
conProjMS(isnan(conProjMS)) = 0;
conProj(isnan(conProj)) = 0;
       
end

function showImg(conProj,conProjMS,ifGB, grainInd,xNum,yNum,X,Y,pathName,group,imgKey)

f1 = figure('Name','Concentration','Visible','Off');
surf(X,Y,conProj)
shading interp
pbaspect ([xNum yNum 1])
axis([min(X(:)) max(X(:)) min(Y(:)) max(Y(:))])
view(2)
hold off
fileNameOut = [pathName '\processImg\' num2str(group) '_Concentration.jpg'];
set(f1, 'CreateFcn', 'set(gcbo,''Visible'',''On'')');
saveas (f1, fileNameOut);

f2 = figure('Name','ConcentrationMS','Visible','Off');
surf(X,Y,conProjMS)
shading interp
pbaspect ([xNum yNum 1])
axis([min(X(:)) max(X(:)) min(Y(:)) max(Y(:))])
view(2)
hold off
fileNameOut = [pathName '\processImg\' num2str(group) '_ConcentrationMS.jpg'];
set(f2, 'CreateFcn', 'set(gcbo,''Visible'',''On'')');
saveas (gca, fileNameOut);

if imgKey == 1
f3 = figure('Name','GB','Visible','Off');
surf(X,Y,ifGB)
shading interp
pbaspect ([xNum yNum 1])
axis([min(X(:)) max(X(:)) min(Y(:)) max(Y(:))])
view(2)
hold off
fileNameOut = [pathName '\processImg\' num2str(group) '_GB.jpg'];
set(f3, 'CreateFcn', 'set(gcbo,''Visible'',''On'')');
saveas (f3, fileNameOut);

f4 = figure('Name','GBInd','Visible','Off');
surf(X,Y,grainInd)
shading interp
pbaspect ([xNum yNum 1])
axis([min(X(:)) max(X(:)) min(Y(:)) max(Y(:))])
view(2)
hold off
fileNameOut = [pathName '\processImg\' num2str(group) '_GBInd.jpg'];
set(f4, 'CreateFcn', 'set(gcbo,''Visible'',''On'')');
saveas (f4, fileNameOut);
end

end

end