function maxdt = CFL_Condition(Mesh,U)

maxdt = min(min(Mesh.u_dx),min(Mesh.u_dy))/U;

end