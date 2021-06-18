function saveLammps(M)

[m,~]=size(M);

pxl = min(M(:,3));
pxh = max(M(:,3));
pyl = min(M(:,4));
pyh = max(M(:,4));
pzl = min(M(:,5));
pzh = max(M(:,5));

fileOut = ['test.s'];
fid=fopen(fileOut,'wt');
fprintf(fid,'%s\n', 'ITEM: TIMESTEP');
fprintf(fid,'%s\n', '0');
fprintf(fid,'%s\n', 'ITEM: NUMBER OF ATOMS');
fprintf(fid,'%d\n', m);
fprintf(fid,'%s\n', 'ITEM: BOX BOUNDS pp pp pp');
fprintf(fid,'%f\t', pxl); fprintf(fid,'%f\n', pxh);
fprintf(fid,'%f\t', pyl); fprintf(fid,'%f\n', pyh);
fprintf(fid,'%f\t', pzl); fprintf(fid,'%f\n', pzh);
fprintf(fid,'%s\n', 'ITEM: ATOMS id type x y z GBtype grain');
for i=1:1:m
     fprintf(fid,'%d\t',M(i,1));
     fprintf(fid,'%d\t',M(i,2));
     fprintf(fid,'%8.6f\t',M(i,3));
     fprintf(fid,'%8.6f\t',M(i,4));
     fprintf(fid,'%8.6f\t',M(i,5));
     fprintf(fid,'%d\t',M(i,6));
     fprintf(fid,'%d\n',M(i,7));
end
end