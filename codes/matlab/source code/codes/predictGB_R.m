%% Generate GB plane mesh. by Rhett, 11/10/19
%  e-mail: x.zhou@mpie.de;
%% This is a function to automatically generate *.obj file for GB concentration mapping.
%% Input parameter: 
% factor,0.5, is the image filter factor
% disMin, 12, is the minimum distance for a line segment (nm).
% BW is the 2D filtered image at the optimized orientation.
% imageSize is the size of each sub-image
% X and Y are the dimensions of the 2D image
% PaTheta is the ideal theta value.
% Default rho value, 0. 
% Then the program requires to select the label file *.txt.
%  The program requires to name an output *.obj file.
%% Output parameter; 
%   A *.obj file to restore the generated GB plane mesh
%% |Try|: predictGB_R(0.5,12,BW,imageSize,X,Y,PaTheta,0);
%% 
function [BWBN,BWBNS,BWBNSS,BWTrF,BWTrRegmax,vert2DRaw,vert2D,lineSeg,CCBTrCentroidRaw,CCBTrCentroid,LabMin,LabMax,resetKey]...
    = predictGB_R(GBIF,GBSkel,minL,maxL,Tr1,Tr2,BWB,BWTr,X,Y,BW,GBGap)

resetKey = 0;
imgSize = size(BW);

%% Find tripple junction

BWTrF = smoothdata(BWTr,1,'gaussian',Tr1);
BWTrF = smoothdata(BWTrF,2,'gaussian',Tr1);
BWTrF = smoothdata(BWTrF,1,'gaussian',Tr2);
BWTrF = smoothdata(BWTrF,2,'gaussian',Tr2);
BWTrF = round(BWTrF);
BWTrF(BWTrF<(mean(nonzeros(BWTrF(:)))/2)) = 0;
BWTrRegmax = imregionalmax(BWTrF);
BWTrRegmax = bwmorph(BWTrRegmax,'thicken',2);

%% Find the cooridate of Tri points from Tri detection
[CCTrCentroidT,CCTrNumPixels] = connComp(BWTrRegmax);
CCTrCentroid(:,1) = CCTrCentroidT(:,2);
CCTrCentroid(:,2) = CCTrCentroidT(:,1);
[~,CCTrOrder] = sort(CCTrNumPixels,'descend'); 
CCTrCentroid = CCTrCentroid(CCTrOrder,:);


%% Find the cooridate of ending or branching points from GB detection

BWBF = smoothdata(BWB,1,'gaussian',2*GBGap);
BWBF = smoothdata(BWBF,2,'gaussian',2*GBGap);

BWBN = BWBF/max(BWBF(:));
BWBNMean = mean(BWBN(:)); 
BWBNStd = std(BWBN(:));
BWBNMax = max(BWBN(:));
BWBNMin = min(BWBN(:));
if BWBNMean+GBIF*BWBNStd >= BWBNMax || BWBNMean+GBIF*BWBNStd <= BWBNMin
  BWBN(BWBN> BWBNMean)= 1;
  BWBN(BWBN<= BWBNMean)= 0;
  resetKey = 1;
else
  BWBN(BWBN>(BWBNMean+GBIF*BWBNStd))= 1;
  BWBN(BWBN<=(BWBNMean+GBIF*BWBNStd))= 0;
end
BWBNS = bwmorph(BWBN,'skel',Inf);

BWBNSS = bwmorph(BWBNS,'spur',GBSkel);
BWBNSS = bwmorph(BWBNSS,'clean');  % newly added, remove isolated pixels.

BWBNBr = bwmorph(BWBNSS,'branchpoints');

BWBNBrD = bwmorph(BWBNBr,'diag');   % add diagnal value;
BWBNBRC = BWBNBrD - BWBNBr;   % matrix of turning points;
BWBNBRC = BWBNBRC.*BWBNSS;     % from two possible point to one

BWBNBRCT = bwmorph(BWBNBRC,'thicken',1); % find 8 neighbors  
BWBNBRCT = bwmorph(BWBNBRCT,'thin',1);  % find 8 neighbors

BWBNBr = BWBNBrD.*(~BWBNBRCT) + BWBNBRC; % All turning points, 1 pixel

BWBNBrCT = bwmorph(BWBNBr,'thicken',1); % find 8 neighbors of all turning points 
BWBNBrCT = bwmorph(BWBNBrCT,'thin',1);  % find 8 neighbors of all turning points

BWBNE = bwmorph(BWBNSS,'endpoints'); % find all the end point, 1 pixel
BWBNN = BWBNE + BWBNBr;  % collection of ending and branch points

[Br(:,1),Br(:,2)] = find (BWBNBr ==1);
[E(:,1),E(:,2)] = find (BWBNE ==1);

% collect points from the end points
k = 1;
while k <= size(E,1)
  [nei3rd,nei3rdOut,Br,E] = find3rdN(E(k,:),Br,E,BWBNSS);
  if ~isnan(nei3rd)
    CCBTrCentroid(k,:) = E(k,:);
    nei3rdColl{k} = nei3rd;
    nei3rdOutColl{k} = nei3rdOut;
    k = k + 1;
  end 
end
k = k-1;
% collect points from the branch points
m = 1;
while m <= size(Br,1)
  [nei3rd,nei3rdOut,Br,E] = find3rdN(Br(m,:),Br,E,BWBNSS);
  if ~isnan(nei3rd)
    CCBTrCentroid(k+m,:) = Br(m,:);
    nei3rdColl{k+m} = nei3rd;
    nei3rdOutColl{k+m} = nei3rdOut;
    m = m + 1;
  end 
end
pNum = k+m-1;   % total points collected.

%% Draw GB lines
l = 0;
for i = 1:pNum-1
  P2 = findLine(i,nei3rdColl,nei3rdOutColl,BWBNSS,CCBTrCentroid);
  if ~isnan(P2)
    for j = 1:size(P2,1)
       l = l+1;
       lineSegIn(l,:) = P2(j,:);
    end
  end
end

%% maxium and minimum GB length check
if exist('lineSegIn','var')
    lineSegM = lineSegIn;
    lineSegM(isnan(lineSegM(:,1)),:) = [];

    LabMin = min(lineSegM(:,7));
    LabMax = max(lineSegM(:,7));

    if ~isnan(maxL)
      if maxL < round(LabMax)
        for i = 1:size(lineSegIn,1)
         if lineSegIn(i,7) > maxL
           lineSegNew = addNewLine(lineSegIn(i,:),maxL,BWBNSS,nei3rdColl);
           if ~isnan(lineSegNew)
             lineSegM(i,:) = [nan,nan,nan,nan,nan,nan,nan,nan];
             lineSegM = [lineSegM;lineSegNew];
           end
         end
        end
      end
      lineSegM(isnan(lineSegM(:,1)),:) = [];
    end

    lineSegC = [lineSegM(:,1:2);lineSegM(:,3:4)];
    [~,ia,~] = unique(lineSegC,'rows');
    CCBTrCentroidTemp = lineSegC(ia,:);

    CCBTrCentroidTemp = [CCBTrCentroidTemp,zeros(size(CCBTrCentroidTemp,1))'];

    for i = 1:size(CCBTrCentroidTemp,1)
      for j = 1:size(lineSegM,1)
         if lineSegM(j,1) == CCBTrCentroidTemp(i,1) && lineSegM(j,2) == CCBTrCentroidTemp(i,2)
            CCBTrCentroidTemp(i,3) = CCBTrCentroidTemp(i,3) + 1;
         end
         if lineSegM(j,3) == CCBTrCentroidTemp(i,1) && lineSegM(j,4) == CCBTrCentroidTemp(i,2)
            CCBTrCentroidTemp(i,3) = CCBTrCentroidTemp(i,3) + 1;
         end
      end
    end

    if ~isnan(minL)
      if minL > round(LabMin)
        for i = 1:size(lineSegM,1)
          if lineSegM(i,7) < minL        
            lineSegM = changeLine(i,lineSegM,CCBTrCentroidTemp);
            lineSegM(i,:) = [nan,nan,nan,nan,nan,nan,nan,nan];
            i = 1;
          end
        end
      end
      lineSegM(isnan(lineSegM(:,1)),:) = [];
    end

    lineSegC = [lineSegM(:,1:2);lineSegM(:,3:4)];
    [~,ia,~] = unique(lineSegC,'rows');
    CCBTrCentroid = lineSegC(ia,:);

    for i = 1:size(lineSegM,1)
        for j = 1:size(CCBTrCentroid,1)
            if lineSegM(i,1) == CCBTrCentroid(j,1) && lineSegM(i,2) == CCBTrCentroid(j,2)
                lineSeg(i,1) = j;
            end
            if lineSegM(i,3) == CCBTrCentroid(j,1) && lineSegM(i,4) == CCBTrCentroid(j,2)
                lineSeg(i,2) = j;
            end
        end
    end


    %% Find a better tripple junction. 
    CCBTrCentroidRaw = CCBTrCentroid;
    CCBTrKey = zeros(size(CCBTrCentroid,1),1);
    for i = 1:size(CCTrCentroid,1)
       for j = 1:size(CCBTrCentroid,1)
           dTr(j) = distance2D(CCTrCentroid(i,:),CCBTrCentroidRaw(j,:));
       end
       if exist('dTr','var')                         % Rhett changed  20200523
          [dTrMin,dTrIndex] = min(dTr);
          if CCBTrKey(dTrIndex) == 0 && dTrMin < 20
              CCBTrKey(dTrIndex) =  dTrMin;
              CCBTrCentroid(dTrIndex,:) = CCTrCentroid(i,:);
          end
       end
    end

    vert2D = [X(1,CCBTrCentroid(:,2))',Y(CCBTrCentroid(:,1),1)];
    vert2DRaw = [X(1,CCBTrCentroidRaw(:,2))',Y(CCBTrCentroidRaw(:,1),1)];

    BWTrRegmax = 255*BWTrRegmax;
else
    vert2D = [];
    vert2DRaw = [];
    lineSeg = [];
    CCBTrCentroidRaw = [];
    CCBTrCentroid = [];
    LabMin = [];
    LabMax = [];
end
if ~exist('lineSeg','var')
    lineSeg = [];
    
end
end









%%
function c = distance2D(a,b)
c = sqrt((b(1) - a(1))^2 + (b(2) - a(2))^2); 
end

%%
function [Centroid,CCnumPixels] = connComp(matrix)
CC = bwconncomp(matrix);
Co = regionprops(CC, 'Centroid' );
Centroid = cat(1,Co.Centroid);
Centroid = round(Centroid);
CCnumPixels = cellfun(@numel,CC.PixelIdxList);
end
%%
function [nei3rd,nei3rdOut,Br,E] = find3rdN(P,Br,E,BWBNSS)
matrixS = size(BWBNSS);
aa = ismember(P,E,'rows');
k = 0;
m = 0;
delete = 0;
for i = -2:2
  for j = -2:2
     if 0 < P(1)+i && P(1)+i <= matrixS(1) && 0 < P(2)+j && P(2)+j <= matrixS(2)
       if ~(i == 0 && j == 0)
         a = [P(1)+i,P(2)+j];
         if aa
           if ismember(a,Br,'rows')
              delete = 1;
              idx = any(E(:,1) == P(1) & E(:,2) == P(2),2);
              E(idx,:) = [];
              break;
           elseif ismember(a,E,'rows')
              delete = 1;
              idx = any(E(:,1) == P(1) & E(:,2) == P(2),2);
              E(idx,:) = [];
              break;
           else
              if BWBNSS(a(1),a(2)) ~=0
                k = k + 1;
                nei3(k,:) = a;
                if abs(i) == 2 || abs(j) == 2 
                  m = m + 1;
                  nei3O(m,:) = a;
                end 
              end
           end
         else
           if ismember(a,Br,'rows')
              idx = any(Br(:,1) == P(1) & Br(:,2) == P(2),2);
              Br(idx,:) = [];
              delete = 1;
              break;
           else
              if BWBNSS(a(1),a(2)) ~=0
                k = k + 1;
                nei3(k,:) = a;
                if abs(i) == 2 || abs(j) == 2 
                  m = m + 1;
                  nei3O(m,:) = a;
                end 
              end
           end
         end
       end
     end
  end
  if delete == 1
      nei3rd = NaN;
      nei3rdOut = NaN;
      break;
  end
end

if delete ~= 1
  if m ~= 0
    nei3rd = [nei3;P];
    nei3rdOut = nei3O;
  else
    idx = any(E(:,1) == P(1) & E(:,2) == P(2),2);
    E(idx,:) = [];
    idx = any(Br(:,1) == P(1) & Br(:,2) == P(2),2);
    Br(idx,:) = [];
    nei3rd = NaN;
    nei3rdOut = NaN; 
  end  
end
end
%%
function pE = findLine(pB,nei3rdColl,nei3rdOutColl,BWBNSS,CCBTrCentroid)
z = 0;
pE = nan;
for i = 1:size(nei3rdOutColl{pB},1) 
    dSeg = 0;
    x = nei3rdOutColl{pB}(i,1);
    y = nei3rdOutColl{pB}(i,2);
    P1 = [x,y];
    P2 = find2nd(P1,nei3rdColl{pB},BWBNSS);
    [a,pE1] = checkEnd(P2,nei3rdColl);
    if ~isnan(P2)
     if a == 1
       P2R = P2;
       P1R = CCBTrCentroid(pB,:);
       dSeg = dSeg + distance2D(P1R,P2);
       P3 = findNext(P2,P1,P1R,BWBNSS);
       [a,pE1] = checkEnd(P3,nei3rdColl);
       if ~isnan(P3)
        infiniteLoop = 0;
        while a
           infiniteLoop = infiniteLoop + 1;
           dSeg = dSeg + distance2D(P2,P3);
           P1 = P2;
           P2 = P3;
           P3 = findNext(P2,P1,P1R,BWBNSS);
           if isnan(P3) | (infiniteLoop > 1000) 
              a = 0;
              pE1 = nan;
           else
              [a,pE1] = checkEnd(P3,nei3rdColl); 
           end
        end
        if ~isnan(P3) & ~isnan(pE1)
           dSeg = dSeg + distance2D(P2,CCBTrCentroid(pE1,:));
        end
        if pE1 > pB
           z = z + 1;
           pE(z,1:2) = CCBTrCentroid(pB,:);
           pE(z,3:4) = CCBTrCentroid(pE1,:);
           pE(z,5:6) = P2R;
           pE(z,7) = dSeg;
           pE(z,8) = distance2D(pE(z,1:2),pE(z,3:4));
        end
       end
     else
       dSeg = dSeg + distance2D(CCBTrCentroid(pB,:),P2);  
       dSeg = dSeg + distance2D(P2,CCBTrCentroid(pE1,:));
       if pE1 > pB
         z = z + 1;
         pE(z,1:2) = CCBTrCentroid(pB,:);
         pE(z,3:4) = CCBTrCentroid(pE1,:);
         pE(z,5:6) = P2;
         pE(z,7) = dSeg;
         pE(z,8) = distance2D(pE(z,1:2),pE(z,3:4));
       end
     end
    end
end
end
%%
function P2 = find2nd(P,nei3rdColl,BWBNSS)
matrixS = size(BWBNSS);
P2 = NaN;
for i = -1:1
  for j = -1:1
    if 0 < P(1)+i && P(1)+i <= matrixS(1) && 0 < P(2)+j && P(2)+j <= matrixS(2)
      if ~(i == 0 && j == 0)
        a = [P(1)+i,P(2)+j];
        if BWBNSS(a(1),a(2)) ~=0
          if ~ismember(a,nei3rdColl,'rows')
             P2 = a;
          end
        end
      end
    end  
  end
end

end
%%
function P3 = findNext(P2,P1,P0,BWBNSS)
matrixS = size(BWBNSS);
P3 = nan;
for i = -1:1
  for j = -1:1
    if 0 < P2(1)+i && P2(1)+i <= matrixS(1) && 0 < P2(2)+j && P2(2)+j <= matrixS(2)
      if ~(i == 0 && j == 0)
        a = [P2(1)+i,P2(2)+j];
        if BWBNSS(a(1),a(2)) ~=0
          if ~(a(1) == P1(1) && a(2) == P1(2)) 
             key = 0;
             for l = -2:2
                 for m = -2:2
                   if a(1) == P0(1)+l && a(2) == P0(2)+m
                       key = 1; 
                   end
                 end
             end
             if key ~=1
               P3 = a;
             end
          end
        end
      end
    end  
  end
end

end
%%
function P3 = findNextNew(P2,P1,BWBNSS)
matrixS = size(BWBNSS);
P3 = nan;
for i = -1:1
  for j = -1:1
    if 0 < P2(1)+i && P2(1)+i <= matrixS(1) && 0 < P2(2)+j && P2(2)+j <= matrixS(2)
      if ~(i == 0 && j == 0)
        a = [P2(1)+i,P2(2)+j];
        if BWBNSS(a(1),a(2)) ~=0
          if ~(a(1) == P1(1) && a(2) == P1(2)) 
            P3 = a;
          end
        end
      end
    end  
  end
end

end
%%
function [a,pE] = checkEnd(P3,nei3rdColl)
a = 1; pE = NaN;
for i = 1:length(nei3rdColl)
   for j = 1:size(nei3rdColl{i},1)
      if nei3rdColl{i}(j,1) == P3(1) && nei3rdColl{i}(j,2) == P3(2)
        a = 0;
        pE = i;
      end   
   end
end
if isnan(P3)
   a = 0; 
end
end
%%
function lineSegN = addNewLine(lineSegIn,maxL,BWBNSS,nei3rdColl)
key = 0;
P1 = lineSegIn(1:2);
P2 = lineSegIn(5:6);
dSeg = distance2D(P1,P2);
P1R = P1;
P3 = findNext(P2,P1,P1R,BWBNSS);
[a,~] = checkEnd(P3,nei3rdColl);
while a
    if dSeg > maxL
      key = key + 1;
      lineSegN(key,1:2) = P1R;
      lineSegN(key,3:4) = P2;
      lineSegN(key,5:6) = [0,0];
      lineSegN(key,7) = dSeg;
      lineSegN(key,8) = distance2D(P1R,P2);
      P1 = P2;
      P2 = P3;
      dSeg = distance2D(P1,P2);
      P1R = P1;
      P3 = findNextNew(P2,P1,BWBNSS);
      [a,~] = checkEnd(P3,nei3rdColl);
    else
      dSeg = dSeg + distance2D(P2,P3);
      P1 = P2;
      P2 = P3;
      P3 = findNextNew(P2,P1,BWBNSS);
      [a,~] = checkEnd(P3,nei3rdColl);
    end
end

dSeg = dSeg + distance2D(P2,lineSegIn(3:4));
if dSeg > maxL
   key = key + 1;
   if key == 1
     lineSegN = nan;
   else
     lineSegN(key,1:2) = P1R;
     lineSegN(key,3:4) = lineSegIn(3:4);
     lineSegN(key,5:6) = [0,0];
     lineSegN(key,7) = dSeg;
     lineSegN(key,8) = distance2D(P1R,lineSegIn(3:4));
   end
else
   if key == 0 || key == 1
     lineSegN = nan;
   else
     lineSegN(key,3:4)= lineSegIn(3:4);
     lineSegN(key,7) = lineSegN(key,7) + dSeg;
     lineSegN(key,8) = distance2D(lineSegN(key-1,3:4),lineSegIn(3:4));
   end  
end



end
%%
function lineSegMax = changeLine(k,lineSegMax,CCBTrCentroidTemp)

for j = 1:size(CCBTrCentroidTemp,1)
   if lineSegMax(k,1) == CCBTrCentroidTemp(j,1) && lineSegMax(k,2) == CCBTrCentroidTemp(j,2)
      pointC1 = CCBTrCentroidTemp(j,3);
   end
   if lineSegMax(k,3) == CCBTrCentroidTemp(j,1) && lineSegMax(k,4) == CCBTrCentroidTemp(j,2)
      pointC2 = CCBTrCentroidTemp(j,3);
   end
end

for i = 1:size(lineSegMax,1)
  if i ~= k
    if (pointC1 == 1 && pointC2 ~= 2) || pointC1 == 2 || (pointC1 > 2 && pointC1 <= pointC2)  
      if lineSegMax(i,1) == lineSegMax(k,1) &&  lineSegMax(i,2) == lineSegMax(k,2) 
         lineSegMax(i,1:2) = lineSegMax(k,3:4);
         lineSegMax(i,5:6) = [nan,nan];
         lineSegMax(i,7) = lineSegMax(i,7) + lineSegMax(k,7);
         lineSegMax(i,8) = distance2D(lineSegMax(i,1:2),lineSegMax(i,3:4));
      end
      if lineSegMax(i,3) == lineSegMax(k,1) &&  lineSegMax(i,4) == lineSegMax(k,2) 
         lineSegMax(i,3:4) = lineSegMax(k,3:4);
         lineSegMax(i,5:6) = [nan,nan];
         lineSegMax(i,7) = lineSegMax(i,7) + lineSegMax(k,7);
         lineSegMax(i,8) = distance2D(lineSegMax(i,1:2),lineSegMax(i,3:4));
      end
    else 
      if lineSegMax(i,1) == lineSegMax(k,3) &&  lineSegMax(i,2) == lineSegMax(k,4) 
         lineSegMax(i,1:2) = lineSegMax(k,1:2);
         lineSegMax(i,5:6) = [nan,nan];
         lineSegMax(i,7) = lineSegMax(i,7) + lineSegMax(k,7);
         lineSegMax(i,8) = distance2D(lineSegMax(i,1:2),lineSegMax(i,3:4));
      end
      if lineSegMax(i,3) == lineSegMax(k,3) &&  lineSegMax(i,4) == lineSegMax(k,4) 
         lineSegMax(i,3:4) = lineSegMax(k,1:2);
         lineSegMax(i,5:6) = [nan,nan];
         lineSegMax(i,7) = lineSegMax(i,7) + lineSegMax(k,7);
         lineSegMax(i,8) = distance2D(lineSegMax(i,1:2),lineSegMax(i,3:4));
      end  
    end
  end
end
end

