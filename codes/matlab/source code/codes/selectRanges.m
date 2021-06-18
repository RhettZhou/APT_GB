function posOut = selectRanges(pos,ranges)
%selects all atoms within a list of ranges of the form:
%[mcbegin mcend]
posOut = [];

for i = 1:size(ranges,1)
  posOut = [posOut; pos(pos(:,4) > ranges(i,1) & pos(:,4) <= ranges(i,2), :)];
end