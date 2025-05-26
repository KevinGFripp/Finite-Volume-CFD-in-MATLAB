function dx_est = EstimateMinCellSize(w,N,L)
dx_est = L/2 * (1+tanh((2/N-1)*w)/tanh(w));
end
