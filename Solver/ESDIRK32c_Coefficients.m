function [ci,a_s,b,b_tilde] = ESDIRK32c_Coefficients()

%Singly diagonally implicit Runge-Kutta methods 
% with an explicit first stage by Anne Kværnø, 2004

%% Coefficients
alpha = 1/2;

a_s = zeros(4,4);

a_s(2,1) = alpha;
a_s(2,2) = alpha;

a_s(3,1) = 5/8;
a_s(3,2) = 3/8;
a_s(3,3) = alpha;

a_s(4,1) = 7/18;
a_s(4,2) = 1/3;
a_s(4,3) = -2/9;
a_s(4,4) = alpha;

ci = zeros(4,1);
ci(2) = 2*alpha;
ci(3) = 3/2;
ci(4) = 1;

b_tilde = zeros(4,1);
b_tilde(1) = 1/2;
b_tilde(2) = 1/2;
b_tilde(3) = 0;

b = zeros(4,1);
b(1) = 7/18;
b(2) = 1/3;
b(3) = -2/9;
b(4) = alpha;



end