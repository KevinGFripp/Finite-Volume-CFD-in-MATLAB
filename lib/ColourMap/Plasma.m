function map = Plasma()

map = zeros(256,3);

Blue = [0 10 90];
Purple = [120 5 100];
Orange = [255 130 0];
Yellow = [255 255 0];
White = [255 255 255];

Indexes = [1 64 128 192 256];

for n = 1:3
map(:,n) = round(interp1(Indexes,[Blue(n) Purple(n) Orange(n) ...
                            Yellow(n) White(n)],1:256,"linear"))./256;
end

end