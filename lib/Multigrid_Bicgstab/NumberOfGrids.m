function numGRIDS = NumberOfGrids(N)

numGRIDS = 1;
Nnew = N;

while (Nnew > 6)
Nnew = (ceil(Nnew/2)-1);
numGRIDS = numGRIDS+1;
end

end