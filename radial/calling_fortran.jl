using DelimitedFiles

cd("/home/alok/.julia/packages/JAC/yvKzj/src/radial/")

struct Dfree
    R::Vector{Float64}      # Radial grid
    P::Vector{Float64}      # Large component
    Q::Vector{Float64}      # Small Component
    Phase::Float64          # Inner Phase shift
    Delta::Float64          # Coulomb Phase shift
end

x2=readdlm("./dhfs079.tab")
R0=convert.(Float64,x2[4:size(x2)[1],1])
RV0=convert.(Float64,x2[4:size(x2)[1],2])
#R0=vec(readdlm("radial.txt"))
#RV0=vec(readdlm("potential.txt"))

"""
Takes the user input radial grid (R), R*V, Energy, kappa
"""
function dfree(R0::Vector{Float64}, RV0::Vector{Float64}, E::Float64, k::Int64)

        npts=size(R0)   # No of radila points
        phase=Ref{Float64}(0.0)
        delta=Ref{Float64}(0.0)
        r=zeros(Float64,25000)
        p=zeros(Float64,25000)
        q=zeros(Float64,25000)

        ccall((:mydfree, "./mod_dfree.so"), Cvoid, 
                (Ref{Int64},Ptr{Float64}, Ptr{Float64}, Ref{Float64}, Ref{Int64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ref{Float64}, Ref{Float64}), 
                npts, R0, RV0, E, k, r, p, q, phase, delta)
        
        #wvfree=readdlm("dirac.dat")
        #nrow=size(wvfree)[1]
        print("phase = ", phase, "delta = ", delta)
        print(typeof(phase),typeof(delta))

        return Dfree(r, p, q, phase.x, delta.x)
end

@time el1 = dfree(R0,RV0,200.,50)
