function map = BlackGreenBlack()

map = zeros(256,3);
Black = [0 0 0];
Green = [16 230 162];


Indexes = [1 128 256];

for n = 1:3
map(:,n) = round(interp1(Indexes,[Black(n) Green(n) Black(n)],1:256,"linear"))./256;
end

end