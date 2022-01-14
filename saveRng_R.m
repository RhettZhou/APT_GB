%%% Written by Rhett Xuyang Zhou 20210524
function saveRng_R(rng, fileName)
if ~exist('fileName','var')
    [file, path] = uiputfile('*.rng','Save rng file');
    fileName = [path file];
end
%% Get the file length
rngLen = length(rng);
atoms = {};
atomLen = 0;
atomsCol = {};

plAtoms = {};
plAtomLen = 0;
plAtomsCol = {};
plAtomsName = {};

mcBegin = [];
mcEnd = [];

plAtomsMcLen = 0;
plAtomMcBegin = [];
plAtomMcEnd = [];
plAtomMcSymb = {};

mcLine = {};
plAtomsMcLine = {};

for i = 1:rngLen
   mcBegin(i) = rng(i).mcbegin; 
   mcEnd(i) = rng(i).mcend;
   atomInIonLen = length(rng(i).atoms);
   atomInIonMax = 0;
   for j = 1:atomInIonLen
       if rng(i).atoms(j).count > 1
           atomInIonMax = 1;
       end
       
       if ~any(strcmp(atoms,rng(i).atoms(j).atom))
           atomLen = atomLen + 1;
           atoms(atomLen) = cellstr(rng(i).atoms(j).atom);
           atomsCol(atomLen) = cellstr([rng(i).atoms(j).atom...
                                ' ' sprintf('%0.2f',rng(i).color(1))...
                                ' ' sprintf('%0.2f',rng(i).color(2))...
                                ' ' sprintf('%0.2f',rng(i).color(3))]);
       end
   end
   if atomInIonLen > 1 || atomInIonMax == 1
       plAtomsMcLen = plAtomsMcLen + 1;
       plAtomMcBegin(plAtomsMcLen) = rng(i).mcbegin; 
       plAtomMcEnd(plAtomsMcLen) = rng(i).mcend;
       symb = split(rng(i).rangeName,' ');
       symb = join(symb,'');
       plAtomMcSymb{plAtomsMcLen} = cellstr(symb);
       if ~any(strcmp(plAtoms,cellstr(symb)))
           plAtomLen = plAtomLen + 1;
           plAtomsName(plAtomLen) = cellstr(rng(i).rangeName);
           plAtoms(plAtomLen) = cellstr(symb);
           plAtomsCol(plAtomLen) = cellstr([char(symb)...
                                ' ' sprintf('%0.2f',rng(i).color(1))...
                                ' ' sprintf('%0.2f',rng(i).color(2))...
                                ' ' sprintf('%0.2f',rng(i).color(3))]);
       end
   end
end

mcMax = max(mcEnd);
[mcHead, mcDigit,mcEleStr] = mcHeadGet(mcMax,atomLen,atoms);

plAtomsMcMax = max(plAtomMcEnd);
[plAtomsMcHead, plAtomMcDigit,plAtomMcEleStr] = mcHeadGet(plAtomsMcMax,plAtomsMcLen,plAtoms);

for i = 1:rngLen
   mc = ['. ' char(sprintf(mcDigit,mcBegin(i))) ' ' char(sprintf(mcDigit,mcEnd(i)))];
   mcEleCount = zeros(1,atomLen);
   atomInIonLen = length(rng(i).atoms);
   mcEleCountMerge = [];
   for j = 1:atomInIonLen
      index = find(mcEleStr==string(rng(i).atoms(j).atom));
      mcEleCount(index) = rng(i).atoms(j).count;
   end
   for j = 1:atomLen
      eleDigit = eleDigitGet(atoms(j));
      mcEleCountStr = string(sprintf(eleDigit,mcEleCount(j)));
      mcEleCountMerge = [char(mcEleCountMerge) ' ' char(mcEleCountStr)];
   end
   mcLine{i} = [char(mc) char(mcEleCountMerge)];
end

for i = 1:plAtomsMcLen
   plAtomsMc = ['. ' char(sprintf(plAtomMcDigit,plAtomMcBegin(i))) ' '...
                     char(sprintf(plAtomMcDigit,plAtomMcEnd(i)))];
   plAtomsMcEleCount = zeros(1,plAtomLen);
   index = find(plAtomMcEleStr==string(plAtomMcSymb{i}));
   plAtomsMcEleCount(index) = 1;
   plAtomsMcEleCountMerge = [];
   for j = 1:plAtomLen
      eleDigit = eleDigitGet(plAtoms(j));
      plAtomsMcEleCountStr = string(sprintf(eleDigit,plAtomsMcEleCount(j)));
      plAtomsMcEleCountMerge = [char(plAtomsMcEleCountMerge) ' ' char(plAtomsMcEleCountStr)];
   end
   plAtomsMcLine{i} = [char(plAtomsMc) char(plAtomsMcEleCountMerge)];
end

%% write the matrix as floats
%% opens the file
fid = fopen(fileName,'wt');
fprintf(fid,'%d%s%d\n', atomLen,' ',rngLen);
for i = 1:atomLen
    fprintf(fid,'%s\n', string(atoms(i)));
    fprintf(fid,'%s\n', string(atomsCol(i)));
end
fprintf(fid,'%s\n', mcHead);
for i = 1:rngLen
    fprintf(fid,'%s\n', string(mcLine(i)));
end
fprintf(fid,'\n');
if plAtomsMcLen ~= 0
    fprintf(fid,'\n');
    fprintf(fid,'%s\n', '--- polyatomic extension');
    fprintf(fid,'%d%s%d\n', plAtomLen,' ',plAtomsMcLen);
    for i = 1:plAtomLen
        fprintf(fid,'%s\n', string(plAtoms(i)));
        fprintf(fid,'%s\n', string(plAtomsCol(i)));
    end
    fprintf(fid,'%s\n', plAtomsMcHead);    
    for i = 1:plAtomsMcLen
        fprintf(fid,'%s\n', string(plAtomsMcLine(i)));
    end
    fprintf(fid,'\n');
end
fclose(fid);
disp('rngfile saved');

end

function [mcHead,mcDigit,mcEleStr] = mcHeadGet(mcMax,atomLen,atoms)
if mcMax < 10
    mcDigit = '%6.4f';
    mcHead = '---------------';
elseif any(mcMax >= 10) && any(mcMax <100)
    mcDigit = '%7.4f';
    mcHead = '-----------------';
elseif any(mcMax >= 100) && any(mcMax <1000)
    mcDigit = '%8.4f';
    mcHead = '-------------------';
else
    mcDigit = '%9.4f';
    mcHead = '---------------------';
end
mcEleStr = [];
for i = 1:length(atoms)
    mcHead = [char(mcHead) ' ' char(atoms(i))];
    mcEleStr = [mcEleStr atoms(i)];
end
end

function eleDigit = eleDigitGet(ele)
    eleDigit = ['%' sprintf('%d', strlength(ele)) 'd'];
end

 



