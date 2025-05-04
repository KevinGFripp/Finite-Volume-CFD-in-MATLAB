function dx_est = EstimateMaxCellSize(w,N,L)
dx_est = L/2 * tanh(2*w/(N))/tanh(w);
end