function obj = exclude3D_R(vert2D,lines,theta)

rho = deg2rad(0);
theta = deg2rad(theta);

depth = 50.0;
numOfV = size(vert2D,1);
numOfL = size(lines,1);
xUp = depth*ones(numOfV,1);
xDw = -depth*ones(numOfV,1);
vert3DUp = [xUp vert2D];
vert3DDw = [xDw vert2D];
vert3D = [vert3DUp;vert3DDw];
faces = [lines(:,1) lines(:,2) lines(:,2)+ numOfV*ones(numOfL,1) lines(:,1)+ numOfV*ones(numOfL,1)];
vert3D = normProj(vert3D,theta,rho);
obj.vertices = vert3D;
obj.faces = faces;

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

end