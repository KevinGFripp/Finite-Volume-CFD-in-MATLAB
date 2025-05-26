function [SVP_u,SVP_v] =StageValuePredictor(uvector2,uvector1, ...
                                            dudt_2,dudt_1,...
                                            vvector2,vvector1, ...
                                            dvdt_2,dvdt_1,...
                                            dt_old)
%% Hermite cubic spline u(t) = at^3 + bt^2 +ct +d 
%% unit time interval of [0,1]
%% Extrapolation is t = 1 + ci dt/dt_old, ci <=1

%% coefficients
SVP_u = struct('a',0,'b',0,'c',0','d',0,'dt_old',dt_old);

SVP_u.a = 2*(uvector1 - uvector2) + dudt_1 + dudt_2;
SVP_u.b = 3*(uvector2 - uvector1) - 2*dudt_1 - dudt_2;
SVP_u.c = dudt_1;
SVP_u.d = uvector1;

SVP_v = struct('a',0,'b',0,'c',0','d',0,'dt_old',dt_old);

SVP_v.a = 2*(vvector1 - vvector2) + dvdt_1 + dvdt_2;
SVP_v.b = 3*(vvector2 - vvector1) - 2*dvdt_1 - dvdt_2;
SVP_v.c = dvdt_1;
SVP_v.d = vvector1;

end