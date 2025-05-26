function [ci,a_s,b,b_tilde] = ESDIRK32_Coefficients()

%Singly diagonally implicit Runge-Kutta methods 
% with an explicit first stage by Anne Kværnø, 2004

%% Coefficients
alpha = 0.4358665215;

a_s = zeros(4,4);

a_s(2,1) = alpha;
a_s(2,2) = alpha;

a_s(3,1) = (-4*alpha*alpha +6*alpha -1)/(4*alpha);
a_s(3,2) = (1-2*alpha)/(4*alpha);
a_s(3,3) = alpha;

a_s(4,1) = (6*alpha - 1)/(12*alpha);
a_s(4,2) = -1/((24*alpha - 12)*alpha);
a_s(4,3) = (6*alpha -1 -6*alpha*alpha)/(6*alpha -3);
a_s(4,4) = alpha;

ci = zeros(4,1);
ci(2) = 2*alpha;
ci(3) = 1;
ci(4) = 1;

b_tilde = zeros(4,1);
b_tilde(1) = (-4*alpha*alpha +6*alpha -1)/(4*alpha);
b_tilde(2) = (1-2*alpha)/(4*alpha);
b_tilde(3) = alpha;

b = zeros(4,1);
b(1) = (6*alpha - 1)/(12*alpha);
b(2) = -1/((24*alpha - 12)*alpha);
b(3) = (6*alpha -1 -6*alpha*alpha)/(6*alpha -3);
b(4) = alpha;



end