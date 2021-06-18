% V2_June 19, 2017 Rhett @UA Developed after Traiviratana, S., University
% of California,San Diego, 2008. This script is used to create
% polycrystalline grain structure Input parameters:
% boxSize(Angstron):(pxl,pyl,pzl,px,py,pz),
% seedType(seed) grainNum(nog),
% lattice(latticei),structureOfCysralline(strui) Seed: 0, Random; 1,
% Default; >1, User defined rType: 0, No Texture; 1, 001 Texture; 2, 011
% Texture; 3, 111 Texture tol: 1/2 of GB width 
% change seed number (seedNum is useless but to change seed number)

function M = mlGB(pxi,pyi,pzi,nog,GBwidthi,latticei,strui,seed,rType)
clearvars -global
global pxl pyl pzl px py pz
global r0 majorlat minorlat lattice
global ano
global ph p f
global nopoly nop nol li
global stru
global gb
global tol GBwidth
pxl = 0.0; pyl = 0.0; pzl = 0.0; px = pxi; py = pyi; pz = pzi;
stru = strui; lattice = latticei; GBwidth = GBwidthi;                                                   

if seed == 0
  rng('shuffle')   % Time seed random
elseif seed  ==1
  rng(0)           % default seed 
else
  rng(seed)        % Any positive number, user defined
end

offset = [-1.0;0.0;1.0];
for i =1:nog
   site{i}.x = randomWithIn(pxl,px);
   site{i}.y = randomWithIn(pyl,py);
   site{i}.z = randomWithIn(pzl,pz);   
end
if rType ~= 0
   seedNum =  randomWithIn(pxl,px); %seedNum =  randomWithIn(pxl,px);       %not use, for generate good texture structure 
end
rectArea = px*py;

if strcmp(stru, 'FCC')
r0 =  0.25*sqrt(2.0*lattice*lattice);
majorlat = lattice;
minorlat = 0.5*lattice;
maxatom = 4.0*(px*py*pz)/(lattice*lattice*lattice) 
elseif strcmp(stru, 'BCC')
r0 = 0.25*sqrt(3.0*lattice*lattice);
majorlat = lattice;
minorlat = 0.5*lattice;
maxatom = 2.0*(px*py*pz)/(lattice*lattice*lattice) 
end

grainSize = 2.0*sqrt(rectArea/(pi*nog))
tol = r0;

k = 0;
for i = 1:3
  for j = 1:3
    for n = 1:nog
      k = k+1 ;
      siteExt{k}.x = site{n}.x+offset(i)*px;
      siteExt{k}.y = site{n}.y+offset(j)*py;
    end
  end
end

nopoly = length(siteExt);
for i = 1:nopoly
  siteExtArray(i,1) = siteExt{i}.x;
  siteExtArray(i,2) = siteExt{i}.y;
end

[pArray,fOld] = voronoin(siteExtArray);
for i = 1:size(pArray,1)
    pOld{i}.x = pArray(i,1);
    pOld{i}.y = pArray(i,2);
end

nopOld = length(pOld);

j = 0;  
for i = 1:nopoly  %Which is the number of grain in exted condition
   if inTheRect(siteExt{i})
     j = j + 1;
     siteCenter(j) = i;
   end
end

noj = 0;
for i = 1:nopOld % Check
   if inTheRect(pOld{i})
     noj = noj + 1;
     junction(noj) = i;          %  not use---------------------------------------
   end
end

nop = 0;

for j = 1:nog % Update the polygons (points) in the center box
   ph{siteCenter(j)}.nop = length(fOld{siteCenter(j)});
   ph{siteCenter(j)}.com = siteExt{siteCenter(j)};
   
   for i = 1:ph{siteCenter(j)}.nop  % assign each point in a polygon
      a = whetherExistInP(pOld{fOld{siteCenter(j)}(i)}); % need to check whether we have the point already
      if  a == 0
        nop = nop + 1;
        p{nop}.x =  pOld{fOld{siteCenter(j)}(i)}.x;
        p{nop}.y =  pOld{fOld{siteCenter(j)}(i)}.y;
        f{siteCenter(j)}(i) = nop;
      else
        f{siteCenter(j)}(i) = a; 
      end
      ph{siteCenter(j)}.p(i) = f{siteCenter(j)}(i);
   end 
   
end

k = 0;
for i = 1:3  % Update the polygons (points) outside the center box
  for j = 1:3
    for n = 1:nog 
      k = k+1;
      if ~(i == 2 && j == 2) 
        ph{k}.nop = ph{siteCenter(n)}.nop;
        ph{k}.com = siteExt{k}; 
        for ii = 1:ph{k}.nop
          pTemp.x = p{ph{siteCenter(n)}.p(ii)}.x + offset(i)*px;
          pTemp.y = p{ph{siteCenter(n)}.p(ii)}.y + offset(j)*py;
          a = whetherExistInP(pTemp);
          if  a == 0
            nop = nop + 1;
            p{nop}.x =  pTemp.x;
            p{nop}.y =  pTemp.y;
            f{k}(ii) = nop;
          else
            f{k}(ii) = a; 
          end
          ph{k}.p(ii) = f{k}(ii);
        end
      end
    end
  end
end

%---------------------lines------------------------------------------------
countLines;  % some line doesn't have a second polygon attached, that value is 0;
nol = length(li);

for j = 1:nog  % Update the polygons (lines) in the center box
  ph{siteCenter(j)}.blank = 1;
  ph{siteCenter(j)}.nol = 0;
  ph{siteCenter(j)}.min.x = NaN; ph{siteCenter(j)}.max.x = NaN;
  ph{siteCenter(j)}.min.y = NaN; ph{siteCenter(j)}.max.y = NaN;
  ph{siteCenter(j)}.min.z = NaN; ph{siteCenter(j)}.max.z = NaN;
  ph{siteCenter(j)}.grain = j;
  gb{siteCenter(j)} = [];
  for i = 1:nol % assgin each line in a polygon
    if li{i}.p(1) == siteCenter(j) || li{i}.p(2) == siteCenter(j)
      ph{siteCenter(j)}.nol = ph{siteCenter(j)}.nol + 1;
      ph{siteCenter(j)}.c{ph{siteCenter(j)}.nol} = center2(siteExt{li{i}.p(2)},siteExt{li{i}.p(1)});
      if li{i}.p(1) == siteCenter(j)
        r = va2b(siteExt{li{i}.p(1)},siteExt{li{i}.p(2)});
        ph{siteCenter(j)}.n{ph{siteCenter(j)}.nol} = unit(r);
        ph{siteCenter(j)}.ngh(ph{siteCenter(j)}.nol) = li{i}.p(2);
      else
        r = va2b(siteExt{li{i}.p(2)},siteExt{li{i}.p(1)});
        ph{siteCenter(j)}.n{ph{siteCenter(j)}.nol} = unit(r);
        ph{siteCenter(j)}.ngh(ph{siteCenter(j)}.nol) = li{i}.p(1); 
      end
    end
  end   
end

k = 0;
for i = 1:3  % Update the polygons (lines) outside the center box
  for j = 1:3
    for n = 1:nog 
      k = k+1;
      if ~(i == 2 && j == 2) 
        ph{k}.blank = 1;
        ph{k}.nol = ph{siteCenter(n)}.nol;
        ph{k}.min.x = NaN; ph{k}.max.x = NaN;
        ph{k}.min.y = NaN; ph{k}.max.y = NaN;
        ph{k}.min.z = NaN; ph{k}.max.z = NaN;
        ph{k}.grain = n;
        gb{k} = [];
        a = [];      
        for ii = 1:ph{k}.nol
          if ii == ph{k}.nol
            iii = 1;  
          else
            iii = ii + 1;
          end
          [a(ii,1),a(ii,2)] = nghPoly(k,ph{k}.p(ii),ph{k}.p(iii));
        end
        a = sortrows(a,1);
        for ii = 1:ph{k}.nol
          ph{k}.c{ii}.x = ph{siteCenter(n)}.c{ii}.x + offset(i)*px;
          ph{k}.c{ii}.y = ph{siteCenter(n)}.c{ii}.y + offset(j)*py;
          ph{k}.n{ii} = ph{siteCenter(n)}.n{ii};
          ph{k}.ngh(ii) = a(ii,2);
        end  
      end
    end
  end
end


%%-Fill the atoms----------------------------------------------------------
ano = 0;

for i = 1:nopoly
  if inTheRect(siteExt{i}) && ph{i}.blank
    or = makeOrientation(rType);
    orR = makeOrientationR(rType);
    fillAt(i,i,or,orR,ph{i}.com,rType);
    for j = 1:ph{i}.prCount
      fillAt(0,i,or,orR,ph{i}.periodic{j},rType);
    end
  end
end

%ano = ano                                     % atom number before deleted
%% Much time, not useful for this purse. There might be overlapped atoms---
%[fill,filledList] = updateFilledList;
% n = 0;
% for i = 1:fill  %  Delete the overlapped atoms
%   for j = 1:ph{filledList(i)}.nol
%     if ~(ph{filledList(i)}.ngh(j) == 0)
%       a = search2Boxes(filledList(i),ph{filledList(i)}.ngh(j));
%       n = n + a;
%     end    
%   end
% end

%%--Output------------------------------------------------------------------


M = recordAtoms(junction,pOld);

end

%% functions
function r = randomWithIn(a, b)
r = rand*(max(a,b)-min(a,b)) + min(a,b);
end

function a = inTheRect(r)
global pxl pyl px py 
if r.x > pxl && r.x < pxl + px
  if r.y > pyl && r.y < pyl + py
    a = 1;
  else
    a = 0;
  end
else
  a = 0;
end
end

function or = makeOrientation(k)
or.rotM = [];
if k == 1  % <100>//Z-axis
  or.rotM = EulerAngles(or.rotM,0,0,0);
  or.or = 1;
elseif k == 2 % <110> // Z-axis
  or.rotM = EulerAngles(or.rotM,0,45,0);
  or.or = 2;
elseif k == 3 % <111>//Z-axis
  theta = 54.735610317245;
  or.rotM = EulerAngles(or.rotM,0,theta,45);
  or.or = 3;
elseif k == 0  % Euler parameter
  a = randomWithIn(-1.0,1.0);
  phi =acos(a);
  b = randomWithIn(-1.0,1.0); u.x = b;
  c = randomWithIn(-1.0,1.0); u.y = c;
  d = randomWithIn(-1.0,1.0); u.z = d;
  u = unit(u);
  b = u.x*sin(phi);
  c = u.y*sin(phi);
  d = u.z*sin(phi);
  or.rotM = EulerParameters(or.rotM,a,b,c,d);
  Z.x = 0; Z.y = 0; Z.z = 1.0;
  Z = rotationT(Z,or.rotM);
  p = unit(Z);
  p = vPositive(p);
  i = 0;
  while ~inTriangle(p) && i < 5
    [p,i] = vReorder(p,i);
  end
  or.or = triangle2Color(p);
  fi = fopen('euler_parameters','at');
  fprintf(fi,'%5d%14.7f%14.7f%14.7f%14.7f%14.7f%14.7f%14.7f%14.7f%14.7f%14.7f\n'...
    ,or.or,a,b,c,d,Z.x,Z.y,Z.z,p.x,p.y,p.z); 
  fclose(fi);
end



end

function or = makeOrientationR(k)
or.rotM = [];
if k == 1  % <100>//Z-axis
  rZRandom = randomWithIn(0,90);
  or.rotM = EulerAngles(or.rotM,rZRandom,0,0);
  or.or = 4;
elseif k == 2 % <110> // Z-axis
  rZRandom = randomWithIn(0,180);
  or.rotM = EulerAngles(or.rotM,rZRandom,0,0);
  or.or = 4;
elseif k == 3 % <111>//Z-axis
  rZRandom = randomWithIn(0,60);
  or.rotM = EulerAngles(or.rotM,rZRandom,0,0);
  or.or = 4;
elseif k == 0 % random, not use
  rZRandom = randomWithIn(0,180);
  or.rotM = EulerAngles(or.rotM,rZRandom,0,0);
  or.or = 4;
end

end


function rotM = EulerAngles(rotM,phi,theta,psi)
phi = phi*pi/180; theta = theta*pi/180; psi = psi*pi/180;

rotM(1,1) = cos(psi)*cos(phi) - cos(theta)*sin(phi)*sin(psi);
rotM(1,2) = cos(psi)*sin(phi) + cos(theta)*cos(phi)*sin(psi);
rotM(1,3) = sin(psi)*sin(theta);
rotM(2,1) = -sin(psi)*cos(phi) - cos(theta)*sin(phi)*cos(psi);
rotM(2,2) = -sin(psi)*sin(phi) + cos(theta)*cos(phi)*cos(psi);
rotM(2,3) = cos(psi)*sin(theta);
rotM(3,1) = sin(theta)*sin(phi);
rotM(3,2) = -sin(theta)*cos(phi);
rotM(3,3) = cos(theta);

end

function rotM = EulerParameters(rotM,a,b,c,d)

rotM(1,1) = a*a + b*b -c*c -d*d;
rotM(1,2) = 2*(b*c + a*d);
rotM(1,3) = 2*(b*d - a*c);
rotM(2,1) = 2*(b*c - a*d);
rotM(2,2) = a*a - b*b + c*c - d*d;
rotM(2,3) = 2*(c*d + a*b);
rotM(3,1) = 2*(b*d + a*c);
rotM(3,2) = 2*(c*d - a*b);
rotM(3,3) = a*a - b*b - c*c + d*d;

end

function u = unit(a) % 2D and 3D 
if length(fieldnames(a)) == 3
l = len(a);
u.x = a.x/l;
u.y = a.y/l;
u.z = a.z/l;
else
l = len(a);
u.x = a.x/l;
u.y = a.y/l;   
end
end  

function l = len(a) % 2D and 3D
if length(fieldnames(a)) == 3
l = a.x*a.x + a.y*a.y + a.z*a.z;
l = sqrt(l);
else
l = a.x*a.x + a.y*a.y;
l = sqrt(l);    
end
end

function u = rotationT(p,r)
u.x = r(1,1)*p.x + r(2,1)*p.y + r(3,1)*p.z;
u.y = r(1,2)*p.x + r(2,2)*p.y + r(3,2)*p.z;
u.z = r(1,3)*p.x + r(2,3)*p.y + r(3,3)*p.z;

end

function a = vPositive(a)
if a.x < 0
    a.x = -1.0*a.x;
end
if a.y < 0
    a.y = -1.0*a.y;
end
if a.z < 0
    a.z = -1.0*a.z;
end
end

function a = inTriangle(v)
% that 100, 110, 111 triangle
r = spherical(v);
o.x = 0.0; o.y = 0.0; o.z = 0.0;    % origin
n.x = 0.0; n.y = -1.0; n.z = 1.0;   % normal = cross(100,111)
a2x = va2b(o,v);
if r.z < pi/4 && dot(a2x,n) <= 0
    a = 1;
else
    a = 0;
end

end

function r = spherical(v)
r.x = len(v);  % radial distance
r.y = acos(v.z/r.x); % zenith angle from  pos-z (north pole); 
r.z = atan(v.y/v.x); % azimuth angle from pos-x
end

function u = va2b(r,s) % 2D and 3D
if length(fieldnames(r)) == 3
u.x = s.x - r.x;
u.y = s.y - r.y;
u.z = s.z - r.z;
else
u.x = s.x - r.x;
u.y = s.y - r.y;   
end
end

function [v,i] = vReorder (v,i)
a = v.x; b =v.y; c = v.z;
if i == 0
  v.x = b; v.y = c; v.z =a; i = 1;
elseif i == 1
  v.x = b; v.y = c; v.z =a; i = 2;
elseif i == 2
  v.x = b; v.y = a; v.z =c; i = 3;
elseif i == 3
  v.x = b; v.y = c; v.z =a; i = 4;
elseif i == 4
  v.x = b; v.y = c; v.z =a; i = 5;
else
  fprintf('Something must be wrong!!\n');
end
end

function d = dot(a,b)% 2D and 3D 
if length(fieldnames(a)) == 3
d = a.x*b.x + a.y*b.y + a.z*b.z;
else
d = a.x*b.x + a.y*b.y;    
end
end

function color = triangle2Color(r)
d100 = unit3(1.0,0.0,0.0);
d110 = unit3(1.0,1.0,0.0);
d111 = unit3(1.0,1.0,1.0);
l1 = int16(15*(1.0 - lenA2B(r,d100)));
l2 = int16(15*(1.0 - lenA2B(r,d110)));
l3 = int16(15*(1.0 - lenA2B(r,d111)));
color = l1*16*16 + l2*16 + l3;

end

function u = unit3(x,y,z)
l = sqrt(x*x + y*y + z*z);
u.x = x/l;
u.y = y/l;
u.z = z/l;
end

function l = lenA2B(r,s) % 2D and 3D
if length(fieldnames(r)) == 3
u.x = s.x - r.x;
u.y = s.y - r.y;
u.z = s.z - r.z;
l = u.x*u.x + u.y*u.y +  u.z*u.z;
l = sqrt(l);
else
u.x = s.x - r.x;
u.y = s.y - r.y;
l = u.x*u.x + u.y*u.y;
l = sqrt(l);    
end

end

function fillAt(j,domain,or,orR,r,rType)
global ph p
global nopoly
global ano
global majorlat minorlat
global stru

noa = 0; ext = 0;   %  need check  Not sure 0/1
extend = [];
rotM = or.rotM;
rotMR = orR.rotM;

if j == 0
  for i = 1:nopoly
    [sum,~] = inPolygon(r,ph{i});
    if sum < 1
      j = i;
      r3D = r;
    end
  end
else
   r3D = exclude3D(r); 
end

if ~(rType == 0)
  rr = rotation(r3D,rotMR);
  rr = rotation(rr,rotM); % move it down, because need to change from 2D to 3D  
else
  rr = rotation(r3D,rotM); % move it down, because need to change from 2D to 3D  
end

if ph{j}.blank
  ph{j}.start = ano;   
  for i = 1:ph{j}.nop  % The boundary of a polygon
    [zLo, zHi] = excludeZ(p{ph{j}.p(i)}); 
    if ~(rType == 0)
      p1 = rotation(zLo,rotMR);
      p2 = rotation(zHi,rotMR);
      p1 = rotation(p1,rotM);
      p2 = rotation(p2,rotM);
    else
      p1 = rotation(zLo,rotM);
      p2 = rotation(zHi,rotM);  
    end
    [ph{j}.min, ph{j}.max] = bound(ph{j}.min,ph{j}.max,p1);
    [ph{j}.min, ph{j}.max] = bound(ph{j}.min,ph{j}.max,p2);
  end
  box = makeBox;
  ph{j}.min  = latSize(ph{j}.min,rr,0); % Make bound multiple of lattice size
  ph{j}.max  = latSize(ph{j}.max,rr,1);
  xyz.z = ph{j}.min.z;
  while xyz.z <= ph{j}.max.z
    xyz.y = ph{j}.min.y;
    while xyz.y <= ph{j}.max.y
      xyz.x = ph{j}.min.x;
      while xyz.x <= ph{j}.max.x
        if ~(rType == 0)
          rt = rotationT(xyz,rotM);
          rt = rotationT(rt,rotMR);          
        else
          rt = rotationT(xyz,rotM);
        end
        [inph,inphGB] = inPolygon(rt,ph{j});
        if inph < 1
          if inPolyhedron(rt,box) < 1
            noa = saveAtom(rt,j,noa,inph,inphGB,domain,ph{j}.grain);
          else
            pr = periodic(rt);
            k = polygonHasThis(pr);
            if notInThisList(extend,ext,k) && ph{k}.blank      % extend
              ext = ext + 1; extend(ext) = k; ph{j}.periodic{ext} = pr; 
            end
          end
        end
        
        if strcmp(stru, 'FCC')
          tmp.x = xyz.x + minorlat;
          tmp.y = xyz.y + minorlat;
          tmp.z = xyz.z;
          if ~(rType == 0)
            rt = rotationT(tmp,rotM);
            rt = rotationT(rt,rotMR);
          else
            rt = rotationT(tmp,rotM);
          end
          [inph,inphGB] = inPolygon(rt, ph{j});
          if inph < 1 && inPolyhedron(rt,box) < 1
            noa = saveAtom(rt,j,noa,inph,inphGB,domain,ph{j}.grain);
          end
      
          tmp.x = xyz.x + minorlat;
          tmp.y = xyz.y;
          tmp.z = xyz.z + minorlat;
          if ~(rType == 0)
            rt = rotationT(tmp,rotM);
            rt = rotationT(rt,rotMR);
          else
            rt = rotationT(tmp,rotM);
          end
          [inph,inphGB] = inPolygon(rt, ph{j});
          if inph < 1 && inPolyhedron(rt,box) < 1
            noa = saveAtom(rt,j,noa,inph,inphGB,domain,ph{j}.grain);
          end
       
          tmp.x = xyz.x;
          tmp.y = xyz.y + minorlat;
          tmp.z = xyz.z + minorlat;
          if ~(rType == 0)
            rt = rotationT(tmp,rotM);
            rt = rotationT(rt,rotMR);
          else
            rt = rotationT(tmp,rotM);
          end
          [inph,inphGB] = inPolygon(rt, ph{j});
          if inph < 1 && inPolyhedron(rt,box) < 1
            noa = saveAtom(rt,j,noa,inph,inphGB,domain,ph{j}.grain);
          end
      
        elseif strcmp(stru, 'BCC')
          tmp.x = xyz.x + minorlat;
          tmp.y = xyz.y + minorlat;
          tmp.z = xyz.z + minorlat;
          if ~(rType == 0)
            rt = rotationT(tmp,rotM);
            rt = rotationT(rt,rotMR);
          else
            rt = rotationT(tmp,rotM);
          end
          [inph,inphGB] = inPolygon(rt, ph{j});
          if inph < 1 && inPolyhedron(rt,box) < 1
            noa = saveAtom(rt,j,noa,inph,inphGB,domain,ph{j}.grain);
          end
        end
        xyz.x = xyz.x + majorlat;
      end
      xyz.y = xyz.y + majorlat;
    end
    xyz.z = xyz.z + majorlat;
  end

  if noa == 0
    ph{j}.blank = 1;
  else
    ph{j}.noa = noa;
    ph{j}.prCount = ext;
    ph{j}.blank = 0;
  end

end

end

function u = rotation(p,r)
u.x = r(1,1)*p.x + r(1,2)*p.y + r(1,3)*p.z;
u.y = r(2,1)*p.x + r(2,2)*p.y + r(2,3)*p.z;
u.z = r(3,1)*p.x + r(3,2)*p.y + r(3,3)*p.z;

end

function countLines
global nopoly
global f
global li

k = 0;
for j = 1:nopoly
  noppf = length(f{j});
  for i = 1:noppf
     k = k + 1; 
     key = 1;
     if i == noppf
        kMin = min(f{j}(i),f{j}(1));
        kMax = max(f{j}(i),f{j}(1));  
     else
        kMin = min(f{j}(i),f{j}(i+1));
        kMax = max(f{j}(i),f{j}(i+1));
     end
     
     if k > 1
     for kCheck = 1:k-1
       if li{kCheck}.p(3) == kMin && li{kCheck}.p(4) == kMax
          k = k - 1;
          li{kCheck}.p(2) = j;
          key = -key;
          break
       end         
     end
     end
     if key == 1
       li{k}.p(1) = j;        
       li{k}.p(3) = kMin; li{k}.p(4) = kMax;
     end
   end   
end
end

function c = center2(r,s) % 2D and 3D
if length(fieldnames(r)) == 3
  c.x = (r.x + s.x)/2.0;
  c.y = (r.y + s.y)/2.0;
  c.z = (r.z + s.z)/2.0;    
else
  c.x = (r.x + s.x)/2.0;
  c.y = (r.y + s.y)/2.0;  
end

end

function [sum,sumGB] = inPolygon(r,ph) % r can be 2D or 3D
global pzl pz
global tol GBwidth
global p

x.x = r.x; x.y = r.y;
sum = 0;
v = NaN;

for i = 1:ph.nol
  a2x = va2b(ph.c{i},x);
  val = dot(a2x,ph.n{i});
  if i == ph.nol
     a = [p{ph.p(i)}.x p{ph.p(i)}.y];
     b = [p{ph.p(1)}.x p{ph.p(1)}.y];
  else
     a = [p{ph.p(i)}.x p{ph.p(i)}.y];
     b = [p{ph.p(i+1)}.x p{ph.p(i+1)}.y];
  end
  c = [x.x x.y]; 
  valH = abs(det([b-a; c-a]))/norm(b-a);
  v = min(v,valH);
  if val >= 0.0
    sum = i + 1;   
  end
end
if v <= tol && sum == 0
  sum = -1;
end

sumGB = sum;                      % Another filter to make GB larger
if v <= GBwidth/2 && sum == 0
  sumGB = -1;
end

if r.z < pzl || r.z > pzl+ pz  % should always in the box region
  sum = 1;
  sumGB = 1;
end

end


function [minV,maxV] = bound(minV,maxV,p) % 2D and 3D
if length(fieldnames(p)) == 3
minV.x = min(minV.x,p.x);
minV.y = min(minV.y,p.y);
minV.z = min(minV.z,p.z);
maxV.x = max(maxV.x,p.x);
maxV.y = max(maxV.y,p.y);
maxV.z = max(maxV.z,p.z);
else
minV.x = min(minV.x,p.x);
minV.y = min(minV.y,p.y);
maxV.x = max(maxV.x,p.x);
maxV.y = max(maxV.y,p.y);   
end

end

function [zLo,zHi] = excludeZ(r)
global pzl pz

zLo.x = r.x; zLo.y = r.y;zLo.z = pzl;
zHi.x = r.x; zHi.y = r.y;zHi.z = pzl + pz;

end

function box = makeBox  % something is questionable
global pxl pyl pzl px py pz
minX = pxl; minY = pyl; minZ = pzl;
maxX = pxl + px; maxY = pyl + py; maxZ = pzl + pz;

pb{1}.x = minX; pb{1}.y = minY; pb{1}.z = minZ;    % the orignal author -0.01 for each term?
pb{7}.x = maxX; pb{7}.y = maxY; pb{7}.z = maxZ;
box.com.x = (minX + maxX)/2.0;
box.com.y = (minY + maxY)/2.0;
box.com.z = (minZ + maxZ)/2.0;
pb{2}.x = pb{7}.x; pb{2}.y = pb{1}.y; pb{2}.z = pb{1}.z;
pb{3}.x = pb{7}.x; pb{3}.y = pb{7}.y; pb{3}.z = pb{1}.z;
pb{4}.x = pb{1}.x; pb{4}.y = pb{7}.y; pb{4}.z = pb{1}.z;
pb{5}.x = pb{1}.x; pb{5}.y = pb{1}.y; pb{5}.z = pb{7}.z;
pb{6}.x = pb{7}.x; pb{6}.y = pb{1}.y; pb{6}.z = pb{7}.z;
pb{8}.x = pb{1}.x; pb{8}.y = pb{7}.y; pb{8}.z = pb{7}.z;
minX = 0; minY = 0; minZ = 0; maxX = 0; maxY = 0; maxZ = 0;

fb{1}.center = center4(pb{2},pb{3},pb{6},pb{7});
fb{2}.center = center4(pb{1},pb{4},pb{5},pb{8});
fb{3}.center = center4(pb{3},pb{4},pb{7},pb{8});
fb{4}.center = center4(pb{1},pb{2},pb{5},pb{6});
fb{5}.center = center4(pb{5},pb{6},pb{7},pb{8});
fb{6}.center = center4(pb{1},pb{2},pb{3},pb{4});
fb{1}.normal = vNormal(pb{2},pb{3},pb{6},box,fb{1});
fb{2}.normal = vNormal(pb{1},pb{4},pb{5},box,fb{2});
fb{3}.normal = vNormal(pb{3},pb{4},pb{7},box,fb{3});
fb{4}.normal = vNormal(pb{1},pb{2},pb{5},box,fb{4});
fb{5}.normal = vNormal(pb{5},pb{6},pb{7},box,fb{5});
fb{6}.normal = vNormal(pb{1},pb{2},pb{3},box,fb{6});

for i = 1:6
  box.c{i} = fb{i}.center;
  box.n{i} = fb{i}.normal;
end

box.nof = 6;

end

function c = center4(r,s,t,u)
c.x = (r.x + s.x +t.x + u.x)/4.0;
c.y = (r.y + s.y +t.y + u.y)/4.0;
c.z = (r.z + s.z +t.z + u.z)/4.0;

end

function n = vNormal(r,s,t,ph,fc)
u = va2b(r,s);
v = va2b(s,t);
n = cross(u,v);
o = va2b(ph.com,fc.center);
d = dot(o,n);
if d < 0.0
n = cross(v,u);
end

end

function u =cross(r,s)
u.x = r.y*s.z - r.z*s.y;
u.y = r.z*s.x - r.x*s.z;
u.z = r.x*s.y - r.y*s.x;

end

function v = latSize(a,r,mm)
global majorlat

v = r;

if mm == 0
  while v.x > a.x 
     v.x = v.x - majorlat;
  end
  while v.y > a.y 
     v.y = v.y - majorlat;
  end
  while v.z > a.z 
     v.z = v.z - majorlat;
  end
elseif mm == 1
  while v.x < a.x 
     v.x = v.x + majorlat;
  end
  while v.y < a.y 
     v.y = v.y + majorlat;
  end
  while v.z < a.z 
     v.z = v.z + majorlat;
  end  
end

end

function u = exclude3D(r)
global pzl pz

u.x = r.x; u.y = r.y; u.z = pzl + pz/2.0;

end

function sum = inPolyhedron(x,p) % 3D only,
sum = 0;
for i = 1:p.nof
  a2x = va2b(p.c{i},x);
  val = dot(a2x,p.n{i});
  if val > 0.0
    sum = i + 1;
  end
end
end

function noa = saveAtom(rt,j,noa,tno,inphGB,domain,grain) % I don't understand how GB is defined
global ano atom
global gb

ano = ano + 1;      
noa = noa + 1;

atom{ano}.r = rt;
atom{ano}.type = tno + 2;                     %  GB =1, internal = 2 For delete overlapp atoms
atom{ano}.GB = inphGB + 2;                     %  GB =1, internal = 2  For indicate GB postions
atom{ano}.ph = domain;
atom{ano}.grain = grain;
atom{ano}.grainNP = j;                        % GB Number without periodic

if tno == -1
  gb{j} = [gb{j} ano];
end     

end

function u = periodic(r)
global pxl pyl px py 
u = r;

if r.x < pxl
  u.x = r.x + px;
elseif r.x > pxl + px
  u.x = r.x - px;
end
if r.y < pyl
  u.y = r.y + py;
elseif r.y > pyl + py
  u.y = r.y - py;
end


end

function j = polygonHasThis(r)
global ph nopoly
for i =1:nopoly
  [sum,~] = inPolygon(r,ph{i});
  if  sum < 1
    j = i;
    break  
  end
end

end

function a = notInThisList(list, size, check)
go = 0;
for i = 1:size
  if check == list(i)
    go = go + 1;
  end
end

if go == 0
    a = 1;
else
    a = 0;
end

end

function [fill,filledList]  = updateFilledList
global nopoly
global ph

fill = 0;
for i =1:nopoly
  if ~ph{i}.blank
     fill = fill + 1;
     filledList(fill) = i;     
  end
end

end

function count = search2Boxes(boxA,boxB)
global gb
global atom
global tol

count = 0;
if ~isempty(gb{boxA})
 for i = 1:length(gb{boxA})
   if ~isempty(gb{boxB})
     for j = 1:length(gb{boxB})
       dist = lenA2B(atom{gb{boxA}(i)}.r,atom{gb{boxB}(j)}.r);
       if dist <= tol
         count = count + 1;
         atom{gb{boxA}(i)}.ph = 0;
         break
       end
     end
   end
 end
end
end

function M = recordAtoms(junction,pOld)
global atom
global ano
global GBwidth

TriR2 = (GBwidth/3)^2;

if TriR2 <= 16
    TriR2 = 16;
end

number = 0;
for i = 1:ano
  if ~(atom{i}.ph == 0)
    number = number + 1;
    for j = 1:length(junction)
       if (pOld{junction(j)}.x - atom{i}.r.x)^2 + (pOld{junction(j)}.y - atom{i}.r.y)^2 <= TriR2
          atom{i}.GB = 3;             %tri = 3; GB =1, internal = 2
       end
    end
    
    M(number,1) = number;
    M(number,2) = 1;                    %M(number,2) = atom{i}.ph;
    M(number,3) = atom{i}.r.x;
    M(number,4) = atom{i}.r.y;
    M(number,5) = atom{i}.r.z;
    M(number,6) = atom{i}.GB;        
    if atom{i}.GB == 2
       M(number,7) = atom{i}.grainNP;
    else
       M(number,7) = 0;
    end
  end
end

end

function a = whetherExistInP(r)
global p
global nop

a = 0;
if ~(nop == 0)
  for i = 1:nop
     if abs(r.x - p{i}.x) < 0.0001 && abs(r.y - p{i}.y) < 0.0001
        a = i;
        break
     else
        a = 0;
     end
  end
end

end

function [i, b] = nghPoly (k,r1i,r2i) 
global nol
global li

r1 = min(r1i,r2i); r2 = max(r1i,r2i);

i = 0;
while  i <= nol
  i = i + 1;
  if r1 == li{i}.p(3) && r2 == li{i}.p(4)
     if k == li{i}.p(1)
       b = li{i}.p(2);
     else
       b = li{i}.p(1);
     end
     break
  end
end
end





















