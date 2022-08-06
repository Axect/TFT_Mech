using DifferentialEquations, NCDataFrame, DataFrames, TensorCast

function damped_sho!(du,u,p,t)
    du[1] = u[2] * exp(-p[1] * t)
    du[2] = -p[2] * u[1] * exp(p[1] * t)
end

u0 = [0.1, 0]
p0 = [0.0, 200]
p = [0.1, 200]
tspan = (0.0, 2.0)

prob1 = ODEProblem(damped_sho!, u0, tspan, p0)
prob2 = ODEProblem(damped_sho!, u0, tspan, p)

sol1 = solve(prob1,alg_hints=[:stiff],reltol=1e-10,abstol=1e-10)
sol2 = solve(prob2,alg_hints=[:stiff],reltol=1e-10,abstol=1e-10)

#plot(sol)

@cast u1[i,j] := sol1.u[i][j]
@cast u2[i,j] := sol2.u[i][j]

df1 = DataFrame(
    :t => sol1.t,
    :x => u1[:,1],
    :v => u1[:,2],
)

df2 = DataFrame(
    :t => sol2.t,
    :x => u2[:,1],
    :v => u2[:,2],
)

writenc(df1, "data/01_sho.nc")
writenc(df2, "data/02_damped_sho.nc")