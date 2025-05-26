function dt = Fixed_dt(maxdt,Re)

dt = min(0.95*maxdt/log(Re),0.25*maxdt);

end