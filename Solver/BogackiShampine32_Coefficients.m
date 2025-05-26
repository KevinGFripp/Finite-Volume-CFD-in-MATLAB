function [ci,a_s,b,b_tilde] = BogackiShampine32_Coefficients()
%% Coefficients

a_s = zeros(4,4);

a_s(2,1) = 1/2;

a_s(3,1) = 0;
a_s(3,2) = 3/4;

a_s(4,1) = 2/9;
a_s(4,2) = 1/3;
a_s(4,3) = 4/9;

b = zeros(4,1);
b(1) = 2/9;
b(2) = 1/3;
b(3) = 4/9;
b(4) = 0;

b_tilde = zeros(4,1);
b_tilde(1) = 7/24;
b_tilde(2) = 1/4;
b_tilde(3) = 1/3;
b_tilde(4) = 1/8;

ci = zeros(4,1);
ci(2) = 1/2;
ci(3) = 3/4;
ci(4) = 1;


end