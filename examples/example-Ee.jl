
#
println("Ef) DR rate coefficients for neon-like Ar: AverageSCA model.")
using JLD
#
if  true
    # Calculation of the Ne^+  2s^2 2p^6; 2s, 2p --> 3s, 3p, 3d, ... 5d electron capture and subsequent autoionization and stabilization
    setDefaults("print summary: open", "Ne-like-DR-rate.sum")
    setDefaults("method: continuum, Galerkin");            setDefaults("method: normalization, pure sine")
    grid = Radial.Grid(Radial.Grid(false), rnt = 1.0e-5, h = 5.0e-2, hp = 2.0e-2, rbox = 15.0)
    name = "DR rate coefficient for neon-like argon and for temperatures < 500 K"
    asfSettings = AsfSettings(AsfSettings(), maxIterationsScf = 36)
    wa   = Cascade.Computation(Cascade.Computation(); name=name, nuclearModel=Nuclear.Model(10.), grid=grid, asfSettings=asfSettings,
                               approach=Cascade.AverageSCA(),
                               scheme=Cascade.DrRateCoefficientScheme([ElecCapture(), Auger(), Radiative()], [E1], 2.5, 0.0, 1, 3, 2, 
                                                                      [Shell("2s"), Shell("2p")]),
                               initialConfigs=[Configuration("1s^2 2s^2 2p^6")] )
    println(wa)
    wb = perform(wa; output=true)
    setDefaults("print summary: close", "")
    #
elseif false
    # Calculation of the Ne^+  2s^2 2p^5; photo-ionization cross sections
    setDefaults("print summary: open", "Ne-plus-photoionization.sum")

    name = "Photoionization of Si- "
    grid = Radial.Grid(Radial.Grid(false); rnt = 3.0e-6, h = 2.0e-2, hp = 3.0e-2, NoPoints=1110)
    wa   = Cascade.Computation(Cascade.Computation(); name=name, nuclearModel=Nuclear.Model(10.), grid=grid, approach=Cascade.AverageSCA(),
                               scheme=Cascade.PhotonIonizationScheme([Photo()], 1, [1.0, 3.0, 5.0]),
                               initialConfigs=[Configuration("1s^2 2s^2 2p^5")] )
    println(wa)
    wb = perform(wa; output=true)
    setDefaults("print summary: close", "")
    
else
    # Simulation of the Ne^+  2s^2 2p^5; photo-ionization cross sections
    setDefaults("print summary: open", "Ne-plus-photoabsorption-simulation.sum")

    data = [JLD.load("zzz-cascade-ionizing-computations-2020-07-10T19.jld"), JLD.load("zzz-cascade-excitation-computations-2020-07-10T19.jld")]
    ##x JLD2.@load "zzz-cascade-ionizing-computations-2020-06-23T09.jld"
    ##x resIon  = results
    ##x JLD2.@load "zzz-cascade-excitation-computations-2020-06-23T09.jld"
    ##x resExc = results
    ##x data = [resExc]

    name = "Simulation of photoabsorption for Ne^+"

    wc   = Cascade.Simulation(Cascade.Simulation(), name=name, property=Cascade.PhotoAbsorption(), 
                              settings=Cascade.SimulationSettings(0., 0., 1., 4., 3., [(1, 1.0)]), computationData=data )
    println(wc)
    wd = perform(wc; output=true)
    setDefaults("print summary: close", "")

end
