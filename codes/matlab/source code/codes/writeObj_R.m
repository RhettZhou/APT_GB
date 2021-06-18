function writeObj_R(patch,OutFileName)
% exports obj file into Wavefront obj.

fid = fopen(OutFileName,'wt');

if( fid == -1 )
    error('Cant open file.');
    return;
end

fprintf(fid, '# Rhett Zhou MPIE 2019 after Peter Felfer \n');

% writing each patch
fprintf(fid, 'v %f %f %f\n', patch.vertices');
% writing faces for individual objects
oName = 's 1\n';
fprintf(fid, oName);
fprintf(fid, 'f %u %u %u %u\n', patch.faces');
fclose(fid);
clear fid;
end