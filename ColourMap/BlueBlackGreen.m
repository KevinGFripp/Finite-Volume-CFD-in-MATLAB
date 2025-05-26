function map = BlueBlackGreen()

map = zeros(256,3);

Blue = [0 120 256];
Black = [0 0 0];
Green = [0 256 120];


Indexes = [1 128 256];

for n = 1:3
map(:,n) = round(interp1(Indexes,[Blue(n) Black(n) Green(n)],1:256,"linear"))./256;
end

end