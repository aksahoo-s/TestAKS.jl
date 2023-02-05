#
println("Dl) Test of the PhotoDoubleIonization module with ASF from an internally generated initial- and final-state multiplet.")
#
setDefaults("print summary: open", "zzz-PhotoDoubleIonization.sum")
setDefaults("method: continuum, Galerkin")           ## setDefaults("method: continuum, Galerkin")  "method: continuum, asymptotic Coulomb"
                                                     ## setDefaults("method: normalization, Ong-Russek") 
setDefaults("method: normalization, pure sine")      ## setDefaults("method: normalization, pure Coulomb")    setDefaults("method: normalization, pure sine")
setDefaults("unit: rate", "a.u.")   

grid = Radial.Grid(Radial.Grid(false), rnt = 4.0e-6, h = 5.0e-2, hp = 1.0e-2, rbox = 10.0)

if  false
    # Green function for single-photon 1s^2 photoionization of helium
    name             = "Photo-double ionization of helium"
    refConfigs       = [Configuration("1s^2")]
    levelSymmetries  = [LevelSymmetry(0, Basics.plus), LevelSymmetry(1, Basics.minus)]
    greenSettings    = GreenSettings(10, [0, 1, 2], 0.1, true, LevelSelection())
    greenexpansion   = GreenExpansion( AtomicState.DampedSpaceCI(), Basics.DeExciteSingleElectron(), levelSymmetries, 2, greenSettings)
    greenRep         = Representation(name, Nuclear.Model(2.), grid, refConfigs, greenexpansion) 
    greenOut         = generate(greenRep, output=true)
    doubleGreen      = greenOut["Green channels"]
    
elseif true
    # Single-photon 1s^2 photoionization of helium
    doubleSettings   = PhotoDoubleIonization.Settings(EmMultipole[E1], UseGauge[UseCoulomb, UseBabushkin], [3., 3.5, 4., 4.5, 5., 6., 7., 8., 9.], 
                                                      doubleGreen, 4, false, true, 3, LineSelection(true, indexPairs=[(1,0)]))
    
    wa = Atomic.Computation(Atomic.Computation(), name="xx", grid=grid, nuclearModel=Nuclear.Model(2.), 
                            initialConfigs  =[Configuration("1s^2")],
                            finalConfigs    =[Configuration("1s^0")], 
                            processSettings = doubleSettings )

    wb = perform(wa)
elseif  false
    # Green function for single-photon 2p^2 photoionization of neon
    name             = "Photo-double ionization of Be-like neon"
    refConfigs       = [Configuration("1s^2 2s^2")]
    levelSymmetries  = [LevelSymmetry(0, Basics.plus), LevelSymmetry(1, Basics.minus)]
    greenSettings    = GreenSettings(10, [0, 1, 2], 0.01, true, LevelSelection())
    greenexpansion   = GreenExpansion( AtomicState.DampedSpaceCI(), Basics.DeExciteSingleElectron(), levelSymmetries, 4, greenSettings)
    greenRep         = Representation(name, Nuclear.Model(5.), grid, refConfigs, greenexpansion) 
    greenOut         = generate(greenRep, output=true)
    doubleGreen      = greenOut["Green channels"]
    
elseif true
    # Single-photon 2p^2 photoionization of neon
    doubleSettings   = PhotoDoubleIonization.Settings(EmMultipole[E1], UseGauge[UseCoulomb, UseBabushkin], [10., 20., 30.], doubleGreen, 2, 
                                                      false, true, 2, LineSelection(true, indexPairs=[(1,2)]))
    
    wa = Atomic.Computation(Atomic.Computation(), name="xx", grid=grid, nuclearModel=Nuclear.Model(5.), 
                            initialConfigs  =[Configuration("1s^2 2s^2")],
                            finalConfigs    =[Configuration("1s 2s")], 
                            processSettings = doubleSettings )

    wb = perform(wa)
    
    
end
setDefaults("print summary: close", "")


