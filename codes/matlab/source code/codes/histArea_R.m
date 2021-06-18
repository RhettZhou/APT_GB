function valuesW = histArea_R(values,bins,weight)
if weight == 1
    weight = ones(numel(values),1);
end

valuesW = zeros(numel(bins),1);

for id = 1:numel(values)
  [~,binId] = histc(values(id),bins);
  if binId ~= 0
     valuesW(binId) = valuesW(binId) + weight(id);
  end 
end

valuesW = valuesW./sum(weight);

end

