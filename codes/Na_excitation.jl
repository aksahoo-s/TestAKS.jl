using Main.JAC

#setDefaults("print summary: open", "zzz-ImpactExcitation.sum")
setDefaults("method: continuum, Galerkin")           ## setDefaults("method: continuum, Galerkin")  "method: continuum, asymptotic Coulomb"
setDefaults("method: normalization, pure sine")      ## setDefaults("method: normalization, pure Coulomb")    setDefaults("method: normalization, pure sine")

grid = Radial.Grid(Radial.Grid(false), rnt = 4.0e-6, h = 5.0e-2, hp = 1.0e-2, rbox = 30.0)

# Calculation of collision strengths and impact excitation cross sections for sodium-like Ar
wa = Atomic.Computation(Atomic.Computation(), name="impact excitation", grid=grid, nuclearModel=Nuclear.Model(11.), 
            initialConfigs  = [Configuration("1s^2 2s^2 2p^6 3s")],
            finalConfigs    = [Configuration("1s^2 2s^2 2p^6 3p")], 
            processSettings = ImpactExcitation.Settings([20.], true, true, true,
                              LineSelection(true, indexPairs=[(1,1), (1,2)]), 0., 18, CoulombInteraction()) )
@show wa
wb = perform(wa)


#
#setDefaults("print summary: close", "")