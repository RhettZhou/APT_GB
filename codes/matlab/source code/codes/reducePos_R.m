function posNew = reducePos_R(pos,percentage)
rng('shuffle')

for i = 1:length(pos)
  a = rand;
  if a <= percentage/100
    ra(i) = 1;
  else
    ra(i) = 0;
  end
end

posTemp = [ra',pos];
posNew = posTemp(posTemp(:,1)==1,2:5);

end