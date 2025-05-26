function map = BlueWhiteRed()

map = zeros(256,3);

Blue = [0 0 256];
White = [256 256 256];
Red = [256 0 0];


Indexes = [1 128 256];

for n = 1:3
map(:,n) = round(interp1(Indexes,[Blue(n) White(n) Red(n)],1:256,"linear"))./256;
end

end