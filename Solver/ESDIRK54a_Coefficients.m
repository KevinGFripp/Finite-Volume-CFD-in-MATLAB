function [ci,a_s,b,b_tilde] = ESDIRK54a_Coefficients()

%% Coefficients
alpha = 0.26;

a_s = zeros(7,7);

a_s(2,1) = alpha;
a_s(2,2) = alpha;

a_s(3,1) = 0.13;
a_s(3,2) = 0.84033320996790809;
a_s(3,3) = alpha;

a_s(4,1) = 0.22371961478320505;
a_s(4,2) = 0.47675532319799699;
a_s(4,3) = -0.06470895363112615;
a_s(4,4) = alpha;

a_s(5,1) = 0.16648564323248321;
a_s(5,2) = 0.10450018841591720;
a_s(5,3) = 0.03631482272098715;
a_s(5,4) = -0.13090704451073998;
a_s(5,5) = alpha;

a_s(6,1) = 0.13855640231268224;
a_s(6,2) = 0;
a_s(6,3) = -0.04245337201752043;
a_s(6,4) = 0.02446657898003141;
a_s(6,5) = 0.61943039072480676;
a_s(6,6) = alpha;

a_s(7,1) = 0.13659751177640291;
a_s(7,2) = 0;
a_s(7,3) = -0.05496908796538376;
a_s(7,4) = -0.04118626728321046;
a_s(7,5) = 0.62993304899016403;
a_s(7,6) = 0.06962479448202728;
a_s(7,7) = alpha;

b = zeros(7,1);
b(1) = a_s(7,1);
b(2) = a_s(7,2);
b(3) = a_s(7,3);
b(4) = a_s(7,4);
b(5) = a_s(7,5);
b(6) = a_s(7,6);
b(7) = a_s(7,7);

b_tilde = zeros(7,1);
b_tilde(1) = a_s(6,1);
b_tilde(2) = a_s(6,2);
b_tilde(3) = a_s(6,3);
b_tilde(4) = a_s(6,4);
b_tilde(5) = a_s(6,5);
b_tilde(6) = a_s(6,6);
b_tilde(7) = 0;

ci = zeros(7,1);
ci(2) = 2*alpha;
ci(3) = sum(squeeze(a_s(3,:)));
ci(4) = sum(squeeze(a_s(4,:)));
ci(5) = sum(squeeze(a_s(5,:)));
ci(6) = sum(squeeze(a_s(6,:)));
ci(7) = sum(squeeze(a_s(7,:)));

end