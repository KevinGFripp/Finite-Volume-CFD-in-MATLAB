function [estimate] = HermiteCubicSplinePredictor(h,SVP)

t = 1 + h/SVP.dt_old;
estimate = (t^3).*SVP.a + (t^2).*SVP.b + t.*SVP.c + SVP.d;

end