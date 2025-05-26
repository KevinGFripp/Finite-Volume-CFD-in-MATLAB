function map = WhiteBlueWhite()

map = zeros(256,3);

Blue = [0.1 0.3 0.6].*256;
White = [256 256 256];


Indexes = [1 78 128 178 256];

for n = 1:3
map(:,n) = round(interp1(Indexes, ...
                         [White(n) Blue(n) White(n) Blue(n) White(n)], ...
                         1:256,"linear"))./256;
end

end