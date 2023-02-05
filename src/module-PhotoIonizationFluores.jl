
"""
`module  JAC.PhotoIonizationFluores`  
    ... a submodel of JAC that contains all methods for computing photo-excitation-autoionization cross 
        sections and rates.
"""
module PhotoIonizationFluores 

    using ..AutoIonization, ..Basics, ..ManyElectron, ..Radial, ..Nuclear, ..PhotoEmission

    """
    `struct  PhotoIonizationFluores.Settings`  
        ... defines a type for the details and parameters of computing photon-impact excitation-autoionization 
            pathways |i(N)>  --> |m(N)>  --> |f(N-1)>.

        + multipoles              ::Array{EmMultipole,1}               ... Specifies the multipoles of the radiation field that are to be included.
        + gauges                  ::Array{UseGauge,1}                  ... Specifies the gauges to be included into the computations.
        + printBefore  ::Bool                               ... True, if all energies and lines are printed before their evaluation.
        + selectPathways          ::Bool                               ... True if particular pathways are selected for the computations.
        + selectedPathways        ::Array{Tuple{Int64,Int64,Int64},1}  ... List of list of pathways, given by tupels (inital, inmediate, final).
        + maxKappa                ::Int64                              ... Maximum kappa value of partial waves to be included.
    """
    struct Settings
        multipoles                ::Array{EmMultipole,1}
        gauges                    ::Array{UseGauge,1} 
        printBefore    ::Bool
        selectPathways            ::Bool
        selectedPathways          ::Array{Tuple{Int64,Int64,Int64},1}
        maxKappa                  ::Int64 
    end 


    """
    `PhotoIonizationFluores.Settings()`  ... constructor for the default values of photon-impact excitation-autoionizaton settings.
    """
    function Settings()
        Settings( Basics.EmMultipole[], UseGauge[], false,  false, Tuple{Int64,Int64,Int64}[], 0)
    end


    # `Base.show(io::IO, settings::PhotoIonizationFluores.Settings)`  
    #		... prepares a proper printout of the variable settings::PhotoIonizationFluores.Settings.  
    function Base.show(io::IO, settings::PhotoIonizationFluores.Settings) 
        println(io, "multipoles:              $(settings.multipoles)  ")
        println(io, "gauges:                  $(settings.gauges)  ")
        println(io, "printBefore:  $(settings.printBefore)  ")
        println(io, "selectPathways:          $(settings.selectPathways)  ")
        println(io, "selectedPathways:        $(settings.selectedPathways)  ")
        println(io, "maxKappa:                $(settings.maxKappa)  ")
    end



    """
    `struct  PhotoIonizationFluores.Channel`  
        ... defines a type for a photon-impact excitaton & autoionization channel that specifies all quantum numbers, phases 
            and amplitudes.

        + excitationChannel  ::PhotoEmission.Channel       ... Channel that describes the photon-impact excitation process.
        + augerChannel       ::AutoIonization.Channel      ... Channel that describes the subsequent Auger/autoionization process.
    """
    struct  Channel
        excitationChannel    ::PhotoEmission.Channel
        augerChannel         ::AutoIonization.Channel
    end 


    """
    `struct  PhotoIonizationFluores.Pathway`  
        ... defines a type for a photon-impact excitation pathway that may include the definition of different excitation 
            and autoionization channels and their corresponding amplitudes.

        + initialLevel        ::Level                  ... initial-(state) level
        + intermediateLevel   ::Level                  ... intermediate-(state) level
        + finalLevel          ::Level                  ... final-(state) level
        + photonEnergy        ::Float64                 ... energy of the (incoming) electron
        + electronEnergy      ::Float64                 ... energy of the (finally outgoing, scattered) electron
        + crossSection        ::EmProperty              ... total cross section of this pathway
        + hasChannels         ::Bool                    ... Determines whether the individual excitation and autoionization 
                                                            channels are defined in terms of their multipole, gauge, 
                                                            free-electron kappa, phases and the total angular momentum/parity 
                                                            as well as the amplitude, or not.
        + channels            ::Array{PhotoIonizationFluores.Channel,1}  ... List of channels of this pathway.
    """
    struct  Pathway
        initialLevel          ::Level
        intermediateLevel     ::Level
        finalLevel            ::Level
        photonEnergy          ::Float64
        electronEnergy        ::Float64
        crossSection          ::EmProperty
        hasChannels           ::Bool
        channels              ::Array{PhotoIonizationFluores.Channel,1}
    end 


    """
    `PhotoIonizationFluores.Pathway()`  
        ... 'empty' constructor for an photon-impact excitation-autoionization pathway between a specified 
            initial, intermediate and final level.
    """
    function Pathway()
        Pathway(Level(), Level(), Level(), 0., 0., 0., false, PhotoIonizationFluores.Channel[] )
    end


    # `Base.show(io::IO, pathway::PhotoIonizationFluores.Pathway)`  
    #		... prepares a proper printout of the variable pathway::PhotoIonizationFluores.Pathway.
    function Base.show(io::IO, pathway::PhotoIonizationFluores.Pathway) 
        println(io, "initialLevel:               $(pathway.initialLevel)  ")
        println(io, "intermediateLevel:          $(pathway.intermediateLevel)  ")
        println(io, "finalLevel:                 $(pathway.finalLevel)  ")
        println(io, "photonEnergy                $(pathway.photonEnergy)  ") 
        println(io, "electronEnergy              $(pathway.electronEnergy)  ")
        println(io, "crossSection:               $(pathway.crossSection)  ")
        println(io, "hasChannels:                $(pathway.hasChannels)  ")
        println(io, "channels:                   $(pathway.channels)  ")
    end



    """
    `PhotoIonizationFluores.computePathways(finalMultiplet::Multiplet, intermediateMultiplet::Multiplet, initialMultiplet::Multiplet, 
                                            grid::Radial.Grid, settings::PhotoIonizationFluores.Settings; output=true)`  
        ... to compute the photo-excitation-fluorescence amplitudes and all properties as requested by the given settings.
            A list of lines::Array{PhotoIonizationFluores.Lines} is returned.
    """
    function  computePathways(finalMultiplet::Multiplet, intermediateMultiplet::Multiplet, initialMultiplet::Multiplet, grid::Radial.Grid, 
                              settings::PhotoIonizationFluores.Settings; output=true)
        println("")
        printstyled("PhotoIonizationFluores.computePathways(): The computation of photo-excitation-fluorescence amplitudes starts now ... \n", color=:light_green)
        printstyled("-------------------------------------------------------------------------------------------------------------------- \n", color=:light_green)
        println("")
        #
        println("Not yet implemented: Data structures and properties still need to be worked out.")
        #
        pathways = "Not yet implemented !"
        #
        if    output    return( pathways )
        else            return( nothing )
        end
    end

end # module
