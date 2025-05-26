function dtnew = AdaptTimeStep(Error,dt,maxdt,TOL,Solver)

dtnew = min(dt * (min(0.92*TOL/Error,1.3))^(1/(Solver.Order)), maxdt);

 dtnew = max(1e-5,dtnew);

end