using ModelingToolkit, MethodOfLines, OrdinaryDiffEq, DomainSets, NPZ



# Defining the variables and their derivatives

@parameters r t
@variables phi(..) Pphi(..)
Dt = Differential(t)
Dr = Differential(r)
Drr = Differential(r)^2



# Defining all parameters

r_min =  0
t_min = -20
r_max = 60
t_max = 60
t0=0
r0=5
l0=1
M0=30
alpha=1
a=2
b=3
G=1


# Defining the parameters of the initial Gaussian

gaussian_peak = r_max-20
width = 1
height = 1





# Defining the initial conditions

phi0(r, t) = height*exp(-width*(r-gaussian_peak)^2)
Pphi0(r, t) = 0




# Defining the velocity and mass functions along with their derivativeas

vv(r,t) = (1/(1+r))*tanh((t - t0)/l0)

M(r,t) = M0*(1+tanh((r - r0 - vv(r,t)*t)/(alpha*l0)))^a * (tanh(r/l0))^b

Mderiv(r,t) = 32*(1-tanh(((t*tanh(t))/(r+1)) - r + 5))*(2*tanh(r^3)*((t*tanh(t)/(r+1)^2)+1))*(sech(((t*tanh(t))/(r+1)) - r + 5))^2 - 3*r^2*(sech(r^3))^3 * (tanh(((t*tanh(t))/(r+1)) - r + 5)-1)

g(r,t) = sqrt(2*G*Mderiv(r,t)/r^2)
    
gderiv(r,t) = (1/2)*(2*G*M(r,t)/r^2)^(-1/2) * (2*G*Mderiv(r,t)*r^2 - 4*G*M(r,t)*r)/(r^4)
    




# Defining the differential equation

eq = [Dt(phi(r,t)) ~ Pphi(r, t)/(r^2) + g(r,t)*Dr(phi(r,t)), 
      Dt(Pphi(r,t)) ~ gderiv(r,t)*Pphi(r,t) + g(r,t)*Dr(Pphi(r,t)) + 2*r*Dr(phi(r,t)) + r^2 * Drr(phi(r,t))]




# Setting the interval in space and time

domains = [r ∈ Interval(r_min, r_max),
           t ∈ Interval(t_min, t_max)]




# Defining the boundary and intial conditions

bcs = [phi(r,t_min) ~ phi0(r,t_min),
       phi(r_max,t) ~ 0,

       Pphi(r,t_min) ~ Pphi0(r,t_min),
       Pphi(r_min,t) ~ 0,
       Pphi(r_max,t) ~ 0] 





# Initializing the PDE system

@named pdesys = PDESystem(eq,bcs,domains,[r,t],[phi(r,t),Pphi(r,t)])






# Defining the number of steps in the discretization

N = 50

order = 2 # This may be increased to improve accuracy of some schemes

# Integers for x and y are interpreted as number of points. Use a Float to directtly specify stepsizes dx and dy.
discretization = MOLFiniteDifference([r=>N], t)






# Convert the PDE problem into an ODE problem
println("Discretization:")
@time prob = discretize(pdesys,discretization)



# Solving the system

println("Solve:")
@time sol = solve(prob, TRBDF2(), saveat=0.1)













discrete_r = sol[r]
discrete_t = sol[t]

solu = sol[phi(r, t)]
solv = sol[Pphi(r , t)]






npzwrite("disc_r.npy", discrete_r)

npzwrite("disc_t.npy", discrete_t)

npzwrite("solu.npy", solu)

npzwrite("solv.npy", solv)



