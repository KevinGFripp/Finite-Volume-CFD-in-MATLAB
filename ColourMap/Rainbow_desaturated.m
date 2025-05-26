function map = Rainbow_desaturated()

map = zeros(256,3);

Black = [0 50 85];
DarkBlue = [0 75 180];
Cyan = [0 255 255];
DarkGreen = [0 140 0];
Yellow = [255 255 0];
Orange = [255 170 0];
Red = [150 0 0];
DarkRed = [255 0 0];

Indexes = [1 40 75 114 150 178 220 256];

RedInterp = interp1(Indexes,[Black(1) DarkBlue(1) Cyan(1) DarkGreen(1) ...
                             Yellow(1) Orange(1) DarkRed(1) Red(1)], ...
                             1:256,"linear");
GreenInterp = interp1(Indexes,[Black(2) DarkBlue(2) Cyan(2) DarkGreen(2) ...
                             Yellow(2) Orange(2) DarkRed(2) Red(2)], ...
                             1:256,"linear");
BlueInterp = interp1(Indexes,[Black(3) DarkBlue(3) Cyan(3) DarkGreen(3) ...
                              Yellow(3) Orange(3) DarkRed(3) Red(3)], ...
                              1:256,"linear");

map(:,1) = round(RedInterp)./256;
map(:,2) = round(GreenInterp)./256;
map(:,3) = round(BlueInterp)./256;

end