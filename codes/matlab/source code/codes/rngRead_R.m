% This function is used to change structure file from web to ovita;
function rng = rngRead_R(fileName)

if ~exist('fileName','var')
    [rngName, pathRngName] = uigetfile('*.rng','Select a file');
    fileName = [pathRngName rngName];
    
%     if ~file
%         disp('no file selected');
%         return
%     end
end

fid = fopen(fileName);
tLine = fgets(fid);

while ~strcmp(tLine(1),'-')
   tLine = fgets(fid);
end

tLineN = tLine;

Count = 0;
rng = {};
[~,tLineN] = strtok(tLineN);
while ~isempty(tLineN)                                                     % split line into tokens, incase, we have multiple interested elements.
  [tmp,tLineN] = strtok(tLineN);
  if ~isempty(tmp)
     Count = Count + 1;                                                    % count tokens
     rng{1,Count+2} = tmp;
  end
end

rngNum = 1;

tLine = fgets(fid);
tLineN = tLine;
[tmpK,tLineN] = strtok(tLineN);

while ~isempty(tLineN)
  Count = 0;
  rngNum = rngNum + 1;
  while ~isempty(tLineN)
    [tmp,tLineN] = strtok(tLineN);
    if ~isempty(tmp)
      Count = Count + 1;                                                   % count tokens
      rng{rngNum,Count} = str2double(tmp);
    end
  end
  tLine = fgets(fid);
  tLineN = tLine;
  [tmpK,tLineN] = strtok(tLineN);
end

fclose(fid);

end
 