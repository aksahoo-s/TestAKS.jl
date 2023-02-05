This is for the calculation of Phase and Radial P(r) and Q(r) 

The inputs are given in "radial.txt" and "potential.txt"
The subroutine "mydfree.f90" is called from the julia file "dfree.jl"
generate a shared binary file with the following command using gfortran

gfortran -fPIC -shared mydfree.f90 radial.f -o mydfree.so

radial.txt  :- containes the radial grid
potential.txt :- containes the R*V(r) data
contants.mod  :- module with all the contants used
dhfs079.tab :- (not required)table for data file

mydfree.f90 :- My modified subroutine which needs to be called from julia
myprogram.f90 :- The basic check for bug free run of mydfree.f90 file
myradialsub.f :- (not required any more)The modified file from the original radial package

myjuliafort.f90 :- the best method to use
calling_fortran.jl :- the best method to use
gfortran -shared -fPIC radial.f dfree.f90 -o mod_dfree.so