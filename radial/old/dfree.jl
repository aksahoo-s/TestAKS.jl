using DelimitedFiles

cd("/home/alok/julia-codes/fortran/radial")

struct Dfree
    R::Vector{Float64}      # Radial grid
    P::Vector{Float64}      # Large component
    Q::Vector{Float64}      # Small Component
    Phase::Float64          # Inner Phase shift
    Delta::Float64          # Coulomb Phase shift
end

#x2=readdlm("dhfs079.tab")
#R0=convert.(Float64,x2[4:size(x2)[1],1])
#RV0=convert.(Float64,x2[4:size(x2)[1],2])
R0=vec(readdlm("radial.txt"))
RV0=vec(readdlm("potential.txt"))

"""
Takes the user input radial grid (R), R*V, Energy, kappa, ch 0/1 for incoming or out going
"""
function dfree(R0::Vector{Float64}, RV0::Vector{Float64}, E::Float64, k::Int64, ch::Int64)

        npts=size(R0)   # No of radila points
        
        ccall((:mydfree_, "./mydfree.so"), Cvoid, 
                (Ref{Int64},Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Int64}, Ref{Int64}), npts, R0, RV0, E, k, ch)
        
        if k < 10
            wvfree=readdlm("dirac_00$(k)_$(ch).dat")
        elseif k <100
            wvfree=readdlm("dirac_0$(k)_$(ch).dat")
        else
            wvfree=readdlm("dirac_$(k)_$(ch).dat")
        end

        nrow=size(wvfree)[1]

        return Dfree(wvfree[6:nrow,1], wvfree[6:nrow,2], wvfree[6:nrow,3],wvfree[2,1], wvfree[4,1])
end

@time el1 = dfree(R0,RV0,200.,8,0)