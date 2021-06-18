function rng = chooseRNG_R(solute,excelRng,atomIon)
% read ranges out of a variable pasted from Excel. 
% default: atomIon = 0: atoms; atomIon = 1: ions;
% Create by Rhett after Felfer, 11/15/2019

%% variable setup
elements = excelRng(1,3:end);
numEl = length(elements);

elementCount = cell2mat(excelRng(2:end,3:end));
ranges = cell2mat(excelRng(2:end,1:2));
numRng = length(ranges(:,1));

[ions, ~, idxc] = unique(elementCount,'rows');
numIon = length(ions(:,1));

rng = [];
%% select the ranges for ions
if strcmp(atomIon,'ion')
  for i = 1:numIon
   %% determining the name of the ion
   ionName = '';
   for el = 1:numEl
       if ions(i,el) > 0
           if ions(i,el) == 1
               ionName = [ionName, elements{el}];
           else
               ionName = [ionName, elements{el}, num2str(ions(i,el))];
           end
       end
   end    
   %% determining the ranges of the solute ion
   if strcmp(solute,ionName)
      rngIdx = idxc == i;
      rng(:,1:2) = ranges(rngIdx,:);
      rng(:,3) = i*ones(size(rng,1),1);
   end     
  end
  rng = [rng,ones(size(rng,1),1)];

%% select the ranges for atoms
else
   for el = 1:numEl
     if strcmp(solute,elements{el})
       count = elementCount(:,el); %multiplicity of element in each range
       rng = [];      
       for k = 1:numRng
         %% determining the ranges of the solute atoms
          if count(k) > 0        
            rng = [rng;ranges(k,:),el,count(k)];
          end
       end
     end
   end
end
   
end

