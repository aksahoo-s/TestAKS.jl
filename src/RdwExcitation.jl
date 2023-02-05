
"""
`module  JAC.RdwExcitation`  
    ... alok modified submodel of JAC that contains all methods for computing electron impact excitation cross sections and collision strengths.
"""
module RdwExcitation 

    using Printf, ..AngularMomentum, ..Basics, ..Continuum, ..Defaults, ..InteractionStrength, ..ManyElectron, 
                  ..Nuclear, ..Radial, ..SpinAngular, ..TableStrings, ..RadialIntegrals, GSL, OffsetArrays

    """
    `struct  RdwExcitationSettings  <:  AbstractProcessSettings`  ... defines a type for the details and parameters of computing electron-impact excitation lines.

        + electronEnergies        ::Array{Float64,1}             ... List of impact-energies of the incoming elecgtrons (in user-defined units).
        + includeBreit            ::Bool                         ... True if the Breit interaction is to be included, and false otherwise.
        + calcCollisionStrength   ::Bool                         ... True, if collision strength need to be calculated, and false otherwise.
        + printBefore             ::Bool                         ... True, if all energies and lines are printed before their evaluation.
        + lineSelection           ::LineSelection                ... Specifies the selected levels, if any.
        + energyShift             ::Float64                      ... An overall energy shift for all transitions |i> --> |f>.
        + maxKappa                ::Int64                        ... Maximum kappa value of partial waves to be included.
        + operator                ::AbstractEeInteraction   
            ... Interaction operator that is to be used for evaluating the e-e interaction amplitudes; allowed values are: 
                CoulombInteraction(), BreitInteraction(), ...
    """
    struct Settings  <:  AbstractProcessSettings
        electronEnergies          ::Array{Float64,1}
        includeBreit              ::Bool 
        calcCollisionStrength     ::Bool
        printBefore               ::Bool 
        lineSelection             ::LineSelection  
        energyShift               ::Float64    
        maxKappa                  ::Int64
        operator                  ::AbstractEeInteraction 
    end 


    """
    `RdwExcitation.Settings()`  ... constructor for the default values of electron-impact excitation line computations.
    """
    function Settings()
       Settings( Float64[], false, false, false, LineSelection(), 0., 0, CoulombInteraction())
    end


    # `Base.show(io::IO, settings::RdwExcitation.Settings)`  ... prepares a proper printout of the variable settings::RdwExcitation.Settings.
    function Base.show(io::IO, settings::RdwExcitation.Settings) 
        println(io, "electronEnergies:           $(settings.electronEnergies)  ")
        println(io, "includeBreit:               $(settings.includeBreit)  ")
        println(io, "calcCollisionStrength:      $(settings.calcCollisionStrength)  ")
        println(io, "printBefore:                $(settings.printBefore)  ")
        println(io, "lineSelection:              $(settings.lineSelection)  ")
        println(io, "energyShift:                $(settings.energyShift)  ")
        println(io, "maxKappa:                   $(settings.maxKappa)  ")
        println(io, "operator:                   $(settings.operator)  ")
    end


    """
    `struct  RdwExcitation.Channel`  
        ... defines a type for a electron-impact excitaiton channel to help characterize the incoming and outgoing (continuum) states of 
            many electron-states with a single free electron

        + initialKappa     ::Int64              ... partial-wave of the incoming free electron
        + finalKappa       ::Int64              ... partial-wave of the outgoing free electron
        + symmetry         ::LevelSymmetry      ... total angular momentum and parity of the scattering state
        + initialPhase     ::Float64            ... phase of the incoming partial wave
        + finalPhase       ::Float64            ... phase of the outgoing partial wave
        + amplitude        ::Complex{Float64}   ... Collision amplitude associated with the given channel.
    """
    struct  Channel
        initialKappa       ::Int64 
        finalKappa         ::Int64 
        symmetry           ::LevelSymmetry
        initialPhase       ::Float64
        finalPhase         ::Float64
        amplitude          ::Complex{Float64}
    end


    # `Base.show(io::IO, channel::RdwExcitation.Channel)`  ... prepares a proper printout of the variable channel::RdwExcitation.Channel.
    function Base.show(io::IO, channel::RdwExcitation.Channel) 
        println(io, "initialKappa:       $(channel.initialKappa)  ")
        println(io, "finalKappa:         $(channel.finalKappa)  ")
        println(io, "symmetry:           $(channel.symmetry)  ")
        println(io, "initialPhase:       $(channel.initialPhase)  ")
        println(io, "finalPhase:         $(channel.finalPhase)  ")
        println(io, "amplitude:          $(channel.amplitude)  ")
    end


    """
    `struct  RdwExcitation.Line`  
        ... defines a type for a electron-impact excitation line that may include the definition of channels and their corresponding amplitudes.

        + initialLevel           ::Level         ... initial- (bound-state) level
        + finalLevel             ::Level         ... final- (bound-state) level
        + initialElectronEnergy  ::Float64       ... energy of the incoming (initial-state) free-electron
        + finalElectronEnergy    ::Float64       ... energy of the outgoing (final-state) free-electron
        + crossSection           ::Float64       ... total cross section of this line
        + collisionStrength      ::Float64       ... total collision strength of this line
        + channels               ::Array{RdwExcitation.Channel,1}  ... List of RdwExcitation channels of this line.
    """
    struct  Line
        initialLevel             ::Level
        finalLevel               ::Level
        initialElectronEnergy    ::Float64
        finalElectronEnergy      ::Float64
        crossSection             ::Float64 
        collisionStrength        ::Float64 
        channels                 ::Array{RdwExcitation.Channel,1}
    end 


    """
    `RdwExcitation.Line()`  ... 'empty' constructor for an electron-impact excitation line between a specified initial and final level.
    """
    function Line()
        Line(Level(), Level(), 0., 0., 0., 0., RdwExcitation.Channel[] )
    end


    """
    `RdwExcitation.Line(initialLevel::Level, finalLevel::Level, crossSection::Float64)`  
        ... constructor for an electron-impact excitation line between a specified initial and final level.
    """
    function Line(initialLevel::Level, finalLevel::Level, crossSection::Float64)
        Line(initialLevel, finalLevel, 0., 0., crossSection, 0., RdwExcitation.Channel[] )
    end


    # `Base.show(io::IO, line::RdwExcitation.Line)`  ... prepares a proper printout of the variable line::RdwExcitation.Line.
    function Base.show(io::IO, line::RdwExcitation.Line) 
        println(io, "initialLevel:            $(line.initialLevel)  ")
        println(io, "finalLevel:              $(line.finalLevel)  ")
        println(io, "initialElectronEnergy:   $(line.initialElectronEnergy)  ")
        println(io, "finalElectronEnergy:     $(line.finalElectronEnergy)  ")
        println(io, "crossSection:            $(line.crossSection)  ")
        println(io, "collisionStrength:       $(line.collisionStrength)  ")
        println(io, "channels:                $(line.channels)  ")
    end


    """
    `RdwExcitation.amplitude(kind::AbstractEeInteraction, channel::RdwExcitation.Channel, cFinalLevel::Level, cInitialLevel::Level, 
                                grid::Radial.Grid; printout::Bool=true)`  
        ... to compute the kind in  CoulombInteraction(), BreitInteraction(), CoulombBreit() electron-impact interaction amplitude 
            <(alpha_f J_f, kappa_f) J_t || O^(e-e, kind) || (alpha_i J_i, kappa_i) J_t>  due to the interelectronic interaction for 
            the given final and initial (continuum) level. A value::ComplexF64 is returned.
    """
    function amplitude(kind::AbstractEeInteraction, channel::RdwExcitation.Channel, cFinalLevel::Level, cInitialLevel::Level, 
                       grid::Radial.Grid; printout::Bool=true)
        nf = length(cFinalLevel.basis.csfs);      fPartial = Subshell(9,channel.finalKappa)         
        ni = length(cInitialLevel.basis.csfs);    iPartial = Subshell(9,channel.initialKappa)    
        
        if  printout  printstyled("Compute ($kind) e-e matrix of dimension $nf x $ni in the final- and initial-state (continuum) bases " *
                                  "for the transition [$(cInitialLevel.index)- ...] " * 
                                  "and for partial waves $(string(fPartial)[2:end]),  $(string(iPartial)[2:end])... ", color=:light_green)    end
        matrix = zeros(ComplexF64, nf, ni)
        #
        ##x @show cInitialLevel.basis.subshells
        ##x @show cFinalLevel.basis.subshells
        if  cInitialLevel.basis.subshells == cFinalLevel.basis.subshells
            iLevel = cInitialLevel;   fLevel = cFinalLevel
        else
            subshells = Basics.merge(cInitialLevel.basis.subshells, cFinalLevel.basis.subshells)
            iLevel    = Level(cInitialLevel, subshells)
            fLevel    = Level(cFinalLevel,   subshells)
        end
        #
        if      kind in [ CoulombInteraction(), BreitInteraction(), CoulombBreit()]        ## pure V^Coulomb interaction
        #--------------------------------------------------------------------------
            for  r = 1:nf
                for  s = 1:ni
                    if  iLevel.basis.csfs[s].J != iLevel.J  ||  iLevel.basis.csfs[s].parity != iLevel.parity      continue    end 
                    subshellList = fLevel.basis.subshells
                    opa  = SpinAngular.TwoParticleOperator(0, plus, true)
                    wa   = SpinAngular.computeCoefficients(opa, fLevel.basis.csfs[r], iLevel.basis.csfs[s], subshellList)
                    #
                    me = 0.
                    for  coeff in wa
                        if   kind in [ CoulombInteraction(), CoulombBreit()]    
                            me = me + coeff.V * InteractionStrength.XL_Coulomb(coeff.nu, 
                                                    fLevel.basis.orbitals[coeff.a], fLevel.basis.orbitals[coeff.b],
                                                    iLevel.basis.orbitals[coeff.c], iLevel.basis.orbitals[coeff.d], grid)   end
                        if   kind in [ BreitInteraction(), CoulombBreit()]    
                            me = me + coeff.V * InteractionStrength.XL_Breit(coeff.nu, 
                                                    fLevel.basis.orbitals[coeff.a], fLevel.basis.orbitals[coeff.b],
                                                    iLevel.basis.orbitals[coeff.c], iLevel.basis.orbitals[coeff.d], grid)   end
                    end
                    matrix[r,s] = me
                end
            end 
            if  printout  printstyled("done. \n", color=:light_green)    end
            amplitude = transpose(fLevel.mc) * matrix * iLevel.mc 
            amplitude = im^( -1.0 * Basics.subshell_l(Subshell(102, channel.finalKappa)) )   * exp(  im*channel.finalPhase )   * 
                        im^Basics.subshell_l(Subshell(101, channel.initialKappa)) * exp(  im*channel.initialPhase ) * amplitude
            @show amplitude
            #
            #
         elseif  kind == "H-E"
        #--------------------
            amplitude = 0.;    error("stop a")
        else    error("stop b")
        end
        
        return( amplitude )
    end


    """
    `RdwExcitation.computeAmplitudesProperties(line::RdwExcitation.Line, nm::Nuclear.Model, grid::Radial.Grid, 
                                                  settings::RdwExcitation.Settings; printout::Bool=true)`  
        ... to compute all amplitudes and properties of the given line; a line::RdwExcitation.Line is returned for which the amplitudes and
            properties are now evaluated.
    """
    function  computeAmplitudesProperties(line::RdwExcitation.Line, nm::Nuclear.Model, grid::Radial.Grid, 
                                          settings::RdwExcitation.Settings; printout::Bool=true)
        newChannels = RdwExcitation.Channel[];   contSettings = Continuum.Settings(false, grid.NoPoints-50);   cross = 0.;   coll = 0.
        amplitude1 = zeros(ComplexF64,181);   coll1 = zeros(ComplexF64,181)
        
        # Define a common subshell list for both multiplets
        subshellList = Basics.generate("subshells: ordered list for two bases", line.finalLevel.basis, line.initialLevel.basis)
        Defaults.setDefaults("relativistic subshell list", subshellList; printout=false)
        
        # First determine a common set of continuum orbitals for the incoming and outgoing electron
        ciOrbitals = Dict{Subshell, Orbital}();     ciPhases = Dict{Subshell, Float64}()
        cfOrbitals = Dict{Subshell, Orbital}();     cfPhases = Dict{Subshell, Float64}()
        for channel in line.channels
            # Generate the continuum orbitals if they were not yet generated before
            iSubshell  = Subshell(101, channel.initialKappa)
            if  !haskey(ciOrbitals, iSubshell)
                newiLevel  = Basics.generateLevelWithSymmetryReducedBasis(line.initialLevel, subshellList)
                #ciOrbital, iPhase     = Continuum.generateOrbitalForLevel(line.initialElectronEnergy, iSubshell, newiLevel, nm, grid, contSettings)
                ciOrbital, iPhase     = generateFreeOrbital(line.initialElectronEnergy, iSubshell, newiLevel, nm, grid, contSettings)
                ciOrbitals[iSubshell] = ciOrbital;      ciPhases[iSubshell] = iPhase
            end
            
            fSubshell  = Subshell(102, channel.finalKappa)
            if  !haskey(cfOrbitals, fSubshell)
                newfLevel  = Basics.generateLevelWithSymmetryReducedBasis(line.finalLevel,   subshellList)
                #cfOrbital, fPhase     = Continuum.generateOrbitalForLevel(line.finalElectronEnergy,   fSubshell, newfLevel, nm, grid, contSettings)
                cfOrbital, fPhase     = generateFreeOrbital(line.finalElectronEnergy,  fSubshell, newfLevel, nm, grid, contSettings)
                cfOrbitals[fSubshell] = cfOrbital;      cfPhases[fSubshell] = fPhase
            end
       	end
        
       	for channel in line.channels
            # Generate two continuum orbitals
            newiLevel  = Basics.generateLevelWithSymmetryReducedBasis(line.initialLevel, subshellList)
            newfLevel  = Basics.generateLevelWithSymmetryReducedBasis(line.finalLevel,   subshellList)
            iSubshell  = Subshell(101, channel.initialKappa)
            fSubshell  = Subshell(102, channel.finalKappa)
            ciOrbital  = ciOrbitals[iSubshell];     iPhase = ciPhases[iSubshell]
            cfOrbital  = cfOrbitals[fSubshell];     fPhase = cfPhases[fSubshell]
            ##x ciOrbital, iPhase  = Continuum.generateOrbitalForLevel(line.initialElectronEnergy, iSubshell, newiLevel, nm, grid, contSettings)
            ##x cfOrbital, fPhase  = Continuum.generateOrbitalForLevel(line.finalElectronEnergy,   fSubshell, newfLevel, nm, grid, contSettings)
            newiLevel = Basics.generateLevelWithExtraElectron(ciOrbital, channel.symmetry, newiLevel)
            newiLevel = Basics.generateLevelWithExtraSubshell(fSubshell,   newiLevel)
            newfLevel = Basics.generateLevelWithExtraSubshell(iSubshell, newfLevel)
            newfLevel = Basics.generateLevelWithExtraElectron(cfOrbital, channel.symmetry, newfLevel)
            ##x @show newiLevel.basis.subshells
            ##x @show newfLevel.basis.subshells
            newChannel = RdwExcitation.Channel(channel.initialKappa, channel.finalKappa, channel.symmetry, iPhase, fPhase, 0.)
            #
            amplitude  = RdwExcitation.amplitude(settings.operator, newChannel, newfLevel, newiLevel, grid, printout=printout)

            coll       = coll + AngularMomentum.bracket([channel.symmetry.J]) * conj(amplitude) * amplitude

            push!( newChannels, RdwExcitation.Channel(newChannel.initialKappa, newChannel.finalKappa, newChannel.symmetry, iPhase, fPhase, amplitude) )
        end
	
        # Calculate the electron-impact excitation strength and cross section
        cross   = coll
        cross = Defaults.convertUnits("cross section: from atomic", coll.re)
        newLine = RdwExcitation.Line( line.initialLevel, line.finalLevel, line.initialElectronEnergy, line.finalElectronEnergy, cross, coll.re, newChannels)
        return( newLine )
    end

    """
    Calculate the amplitude
    """
    function getAmplitude(line::RdwExcitation.Line, iSubshell::Subshell, fSubshell::Subshell, 
        newiLevel::Level, newfLevel::Level, amplitude::ComplexF64)

        matx = OffsetArray(zeros(ComplexF64, 181), 0:180)

        Ji = line.initialLevel.J.num // line.initialLevel.J.den
        Jf = line.finalLevel.J.num // line.finalLevel.J.den

        for θ = 0:180
            sum = 0
            for Mi in -Ji:Ji
                for Mf in -Jf:Jf
                    for μa in -1//2:1//2
                        for μb in -1//2:1//2
                            sum = amplitude * calcCoeff( line, iSubshell, fSubshell, newiLevel, newfLevel, Mi, Mf, μa, μb, float(θ) )
                            matx[θ] = matx[θ] + conj(sum) * sum
                        end
                    end
                end
            end
        end

        return matx

    end

    """
    Calculate total coefficient by summing over Mi, Mf, μa, μb
    """
    function totCoeff( line::RdwExcitation.Line, iSubshell::Subshell, fSubshell::Subshell, 
        newiLevel::Level, newfLevel::Level, θ::Float64)
        sum = 0.0
        Ji = line.initialLevel.J.num // line.initialLevel.J.den
        Jf = line.finalLevel.J.num // line.finalLevel.J.den
        for Mi in -Ji:Ji
            for Mf in -Jf:Jf
                for μa in -1//2:1//2
                    for μb in -1//2:1//2
                        sum += calcCoeff( line, iSubshell, fSubshell, newiLevel, newfLevel, Mi, Mf, μa, μb, θ )
                    end
                end
            end
        end
        return sum
    end

    
    """
    This function calculated the coefficients in formula 1 
    Here the mla is 0 as incident is quantized along z-Axis
    which makes the ma = mla + μa = 0 + μa = μa

    """
    function calcCoeff( line::RdwExcitation.Line, iSubshell::Subshell, fSubshell::Subshell, 
                    newiLevel::Level, newfLevel::Level, Mi::Rational{Int64}, Mf::Rational{Int64},
                     μa::Rational{Int64}, μb::Rational{Int64}, θ::Float64 )

        Ea = line.initialElectronEnergy;    
        Eb = line.finalElectronEnergy
        c  = 1/Defaults.FINE_STRUCTURE_CONSTANT             # Speed of light in atomic unit
        la = Basics.subshell_l(iSubshell);  
        lb = Basics.subshell_l(fSubshell)
        ja = Basics.subshell_j(iSubshell).num // Basics.subshell_j(iSubshell).den
        jb = Basics.subshell_j(fSubshell).num // Basics.subshell_j(fSubshell).den
        #ma = ja
        Ji = line.initialLevel.J.num // line.initialLevel.J.den
        Jf = line.finalLevel.J.num // line.finalLevel.J.den
        Jt = newiLevel.J.num // newiLevel.J.den

        ## Mi = line.initialLevel.M.num // line.initialLevel.M.den
        ## Mf = line.finalLevel.M.num // line.finalLevel.M.den

        ##μa = 1//2; μb = 1//2; 

        cof = 0.5 * (π^(-3/2)) * (Ea + c^2 ) * (Eb + c^2 ) / (Ea * Eb) * sqrt((2la + 1)/4π)
        cof = cof * MyClebschGordan( la, 0, 1//2, μa, ja, μa ) 

        sum = 0.0
        for mb in -jb:jb           
            for mlb in -lb:lb
                sum2 = 0.0
                for Mt = -Jt:Jt
                    sum2 += MyClebschGordan(Ji, Mi, ja, μa, Jt, Mt) * MyClebschGordan(Jf, Mf, jb, mb, Jt, Mt)
                end
                sum2 *= MyClebschGordan( lb, mlb, 1//2, μb, jb, μb) * Ylmθ(lb, mlb, θ, 0.0)

                sum += sum2
            end
        end
        return sum
    end

    """
    Ylm(l,m,θ,ϕ) = GSL.sf_legendre_sphPlm(l, m, cos(θ)) * ℯ^(im*m*ϕ)
    This function Calculates the spherical harmonics as the default JAC can not above l > 4
    """
    function Ylmθ(l::Int64, m::Int64, θ::Float64, ϕ::Float64)
        if m >= 0
            return GSL.sf_legendre_sphPlm(l, m, cos(θ)) * ℯ^(im*m*ϕ)
        else
            m = -m
            return (-1)^m * conj(GSL.sf_legendre_sphPlm(l, m, cos(θ)) * ℯ^(im*m*ϕ))
        end
    end

    """
    My ClebschGordan as the JAC one gives a negative sign sometimes
    (-1)^(ja-jb+M)*sqrt(2j+1)*Wigner_3j
    """
    function MyClebschGordan(ja, ma, jb, mb, Jab, Mab)
        mab = - Basics.twice(Mab) / 2
        pp  = (ja - jb + Mab)
        cg  = (-1)^pp * sqrt(Basics.twice(Jab) + 1) * sf_coupling_3j(Basics.twice(ja), Basics.twice(jb), Basics.twice(Jab),
                                                                     Basics.twice(ma), Basics.twice(mb), Basics.twice(mab))
        return( cg )
    end

    """
    `RdwExcitation.computeLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet, nm::Nuclear.Model, grid::Radial.Grid, 
                                   settings::RdwExcitation.Settings; output=true)`  
        ... to compute the electron-impact excitation transition amplitudes and all properties as requested by the given settings. 
            A list of lines::Array{RdwExcitation.Lines} is returned.
    """
    function  computeLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet, nm::Nuclear.Model, grid::Radial.Grid, settings::RdwExcitation.Settings; 
                           output=true)
        println("")
        printstyled("RdwExcitation.computePathways(): The computation of electron-impact excitation cross sections starts now ... \n", color=:light_green)
        printstyled("--------------------------------------------------------------------------------------------------------------- \n", color=:light_green)
        println("")
        #
        lines = RdwExcitation.determineLines(finalMultiplet, initialMultiplet, settings)
        # Display all selected lines before the computations start
        if  settings.printBefore    RdwExcitation.displayLines(lines)    end
        # Calculate all amplitudes and requested properties
        newLines = RdwExcitation.Line[]
        for  line in lines
            newLine = RdwExcitation.computeAmplitudesProperties(line, nm, grid, settings) 
            push!( newLines, newLine)
        end
        # Print all results to screen
        RdwExcitation.displayResults(newLines)
        #
        if    output    return( newLines )
        else            return( nothing )
        end
    end


    """
    `RdwExcitation.computeMatrix(finalBasis::Basis, initialBasis::Basis, settings::RdwExcitation.Settings)`  
        ... to compute the transition matrix  (<finalContinuumCSF_r|| V(e-e) ||initialContinuumCSF_s>)  between the CSF_r from the 
            finalContinuumBasis and the CSF_s from the initialContinuumBasis. A (non-quadratic) matrix::Array{Float64,2} with 
            dimensions [length(finalContinuumBasis.csfs) x length(initialContinuumBasis.csfs)] is returned. Note that this transition 
            matrix is typically specific to just one Eimex channel due to the different energies, partial waves and overall symmetry 
            of the scattering states. **Not yet implemented !**
    """
    function computeMatrix(finalBasis::Basis, initialBasis::Basis, settings::RdwExcitation.Settings)   
        error("Not yet implemented.")
    end



    """
    `RdwExcitation.determineChannels(finalLevel::Level, initialLevel::Level, settings::RdwExcitation.Settings)`  
        ... to determine a list of electron-impact excitation Channels for a transitions from the initial to the final level and by 
            taking into account the particular settings of for this computation; an Array{RdwExcitation.Channel,1} is returned.
    """
    function determineChannels(finalLevel::Level, initialLevel::Level, settings::RdwExcitation.Settings)
        channels = RdwExcitation.Channel[];   
        symi = LevelSymmetry(initialLevel.J, initialLevel.parity);    symf = LevelSymmetry(finalLevel.J, finalLevel.parity) 
        #
        for  inKappa = -(settings.maxKappa+1):settings.maxKappa
            if  inKappa == 0    continue    end
            symtList = AngularMomentum.allowedTotalSymmetries(symi, inKappa)
            for  symt in symtList
                outKappaList = AngularMomentum.allowedKappaSymmetries(symt, symf)
                for  outKappa in outKappaList
                    #if abs(outKappa) > settings.maxKappa continue end  #alok
                    push!(channels, RdwExcitation.Channel( inKappa, outKappa, symt, 0., 0.,Complex(0.)) )
                end
            end
        end
        return( channels )  
    end


    """
    `RdwExcitation.determineLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet, settings::RdwExcitation.Settings)`  
        ... to determine a list of RdwExcitation.Line's for transitions between levels from the initial- and final-state multiplets, 
            and by taking into account the particular selections and settings for this computation; an Array{RdwExcitation.Line,1} is 
            returned. Apart from the level specification, all physical properties are set to zero during the initialization process.
    """
    function  determineLines(finalMultiplet::Multiplet, initialMultiplet::Multiplet, settings::RdwExcitation.Settings)
        energyShift = Defaults.convertUnits("energy: to atomic", settings.energyShift)
        lines = RdwExcitation.Line[]
        for  iLevel  in  initialMultiplet.levels
            for  fLevel  in  finalMultiplet.levels
                if  Basics.selectLevelPair(iLevel, fLevel, settings.lineSelection)
                    for  en in settings.electronEnergies
                        initialElectronEnergy  = Defaults.convertUnits("energy: to atomic", en)
                        finalElectronEnergy    = initialElectronEnergy - (fLevel.energy - iLevel.energy) + energyShift
                        if  finalElectronEnergy < 0    continue   end  
                        channels = RdwExcitation.determineChannels(fLevel, iLevel, settings) 
                        push!( lines, RdwExcitation.Line(iLevel, fLevel, initialElectronEnergy, finalElectronEnergy, 0., 0., channels) )
                    end
                end
            end
        end
        return( lines )
    end


    """
    `RdwExcitation.displayLines(lines::Array{RdwExcitation.Line,1})`  
        ... to display a list of lines and channels that have been selected due to the prior settings. A neat table of all 
            selected transitions and energies is printed but nothing is returned otherwise.
    """
    function  displayLines(lines::Array{RdwExcitation.Line,1})
        nx = 180
        println(" ")
        println("  Selected electron-impact ionization lines:")
        println(" ")
        println("  ", TableStrings.hLine(nx))
        sa = "  ";   sb = "  "
        sa = sa * TableStrings.center(18, "i-level-f"; na=2);                                sb = sb * TableStrings.hBlank(20)
        sa = sa * TableStrings.center(18, "i--J^P--f"; na=3);                                sb = sb * TableStrings.hBlank(22)
        sa = sa * TableStrings.center(12, "f--Energy--i"; na=4)               
        sb = sb * TableStrings.center(12, TableStrings.inUnits("energy"); na=5)
        sa = sa * TableStrings.center(12, "Energy e_in"; na=3);              
        sb = sb * TableStrings.center(12, TableStrings.inUnits("energy"); na=3)
        sa = sa * TableStrings.center(12, "Energy e_out"; na=4);              
        sb = sb * TableStrings.center(12, TableStrings.inUnits("energy"); na=4)
        sa = sa * TableStrings.flushleft(57, "List of partial waves and total symmetries"; na=4)  
        sb = sb * TableStrings.flushleft(57, "partial-in [total J^P] partial-out        "; na=4)
        println(sa);    println(sb);    println("  ", TableStrings.hLine(nx)) 
        #
        NoLines = 0   
        for  line in lines
            NoLines = NoLines + 1
            sa  = "  ";    isym = LevelSymmetry( line.initialLevel.J, line.initialLevel.parity)
                           fsym = LevelSymmetry( line.finalLevel.J,   line.finalLevel.parity)
            sa = sa * TableStrings.center(18, TableStrings.levels_if(line.initialLevel.index, line.finalLevel.index); na=2)
            sa = sa * TableStrings.center(18, TableStrings.symmetries_if(isym, fsym); na=4)
            en = line.finalLevel.energy - line.initialLevel.energy
            sa = sa * @sprintf("%.6e", Defaults.convertUnits("energy: from atomic", en))                          * "    "
            sa = sa * @sprintf("%.6e", Defaults.convertUnits("energy: from atomic", line.initialElectronEnergy))  * "    "
            sa = sa * @sprintf("%.6e", Defaults.convertUnits("energy: from atomic", line.finalElectronEnergy))    * "    "
            kappaInOutSymmetryList = Tuple{Int64,Int64,LevelSymmetry}[]
            for  i in 1:length(line.channels)
                push!( kappaInOutSymmetryList, (line.channels[i].initialKappa, line.channels[i].finalKappa, line.channels[i].symmetry) ) 
            end
            wa = TableStrings.kappaKappaSymmetryTupels(95, kappaInOutSymmetryList)
            sb = sa * wa[1];    println( sb )  
            for  i = 2:length(wa)
                NoLines = NoLines + 1
                sb = TableStrings.hBlank( length(sa) ) * wa[i];    println( sb )
            end
            # Avoid long table printouts
            if NoLines > 100  println("\n  RdwExcitation.displayLines():  A maximum of 100 lines are printed in this table. \n")
               break   
            end
        end
        println("  ", TableStrings.hLine(nx))
        #
        return( nothing )
    end


    """
    `RdwExcitation.displayResults(lines::Array{RdwExcitation.Line,1})`  
        ... to list all results, energies, cross sections, etc. of the selected lines. A neat table is printed but nothing is 
            returned otherwise.
    """
    function  displayResults(lines::Array{RdwExcitation.Line,1})
        nx = 131
        println(" ")
        println("  Electron-impact excitation cross sections:")
        println(" ")
        println("  ", TableStrings.hLine(nx))
        sa = "  ";   sb = "  "
        sa = "  ";   sb = "  "
        sa = sa * TableStrings.center(18, "i-level-f"; na=2);                                sb = sb * TableStrings.hBlank(20)
        sa = sa * TableStrings.center(18, "i--J^P--f"; na=3);                                sb = sb * TableStrings.hBlank(22)
        sa = sa * TableStrings.center(12, "f--Energy--i"; na=4)               
        sb = sb * TableStrings.center(12, TableStrings.inUnits("energy"); na=5)
        sa = sa * TableStrings.center(12, "Energy e_in"; na=3);              
        sb = sb * TableStrings.center(12, TableStrings.inUnits("energy"); na=3)
        sa = sa * TableStrings.center(12, "Energy e_out"; na=3);              
        sb = sb * TableStrings.center(12, TableStrings.inUnits("energy"); na=3)
        sa = sa * TableStrings.center(15, "Cross section"; na=3)      
        sb = sb * TableStrings.center(15, TableStrings.inUnits("cross section"); na=3)
        sa = sa * TableStrings.center(15, "Collision strength"; na=3)      
        sb = sb * TableStrings.center(15, " ";                  na=3)
        println(sa);    println(sb);    println("  ", TableStrings.hLine(nx)) 
        #   
        for  line in lines
            sa  = "  ";    isym = LevelSymmetry( line.initialLevel.J, line.initialLevel.parity)
                           fsym = LevelSymmetry( line.finalLevel.J,   line.finalLevel.parity)
            sa = sa * TableStrings.center(18, TableStrings.levels_if(line.initialLevel.index, line.finalLevel.index); na=2)
            sa = sa * TableStrings.center(18, TableStrings.symmetries_if(isym, fsym); na=4)
            en = line.finalLevel.energy - line.initialLevel.energy
            sa = sa * @sprintf("%.6e", Defaults.convertUnits("energy: from atomic", en))                          * "    "
            sa = sa * @sprintf("%.6e", Defaults.convertUnits("energy: from atomic", line.initialElectronEnergy))  * "    "
            sa = sa * @sprintf("%.6e", Defaults.convertUnits("energy: from atomic", line.finalElectronEnergy))    * "     "
            #sa = sa * @sprintf("%.6e", Defaults.convertUnits("cross section: from atomic", line.crossSection))    * "        "
            sa = sa * @sprintf("%.6e", line.crossSection)    * "        "
            sa = sa * @sprintf("%.6e", line.collisionStrength)                                                    * "    "
            println(sa)
        end
        println("  ", TableStrings.hLine(nx))
        #
        return( nothing )
    end

    """
    RdwExcitation.generateFreeOrbital()
    genrates distorted wavefunction for projectile electron
    returns cOrbital, iphase, cphase
    iphase is the inner phase shift
    cphase is the Coulomb phase shift
    """
    function generateFreeOrbital(energy::Float64, sh::Subshell, level::Level, nm::Nuclear.Model, grid::Radial.Grid, settings::Continuum.Settings)
        
        # Generate a (local) potential for the given level
        nuclearPotential  = Nuclear.nuclearPotential(nm, grid)
        wp = compute("radial potential: Kohn-Sham", grid, level)   
        #wp = computePotential(grid, level)

        pot = Basics.add(nuclearPotential, wp)
        #pot = - nuclearPotential.Zr + wp.Zr
        #pot = - ( nuclearPotential.Zr ./ grid.r) + wp

        # Obtaining static potential
        #pot = computePotential(grid, nm, level)

        # rV for input to dfree (radial.f)
        #pot = pot .* grid.r

        # Calling the radial fortran program
        p, q, iphase, cphase = dfree(grid.r, -pot.Zr, energy, sh.kappa)
        corbital = Orbital(sh, false, true, energy, p, q, zeros(Float64, grid.NoPoints), zeros(Float64, grid.NoPoints), grid)
        
        if level.basis.NoElectrons - convert(Int64,nm.Z) == 0
            phase = iphase
        else
            phase = iphase + cphase
        end

        return corbital, phase

    end
    

    """
    (Not being used)
    Compute static potential for Distorted wave
        changed from rV to only V; r will the multiplied latter
    """
    function computePotential(grid::Radial.Grid, nm::Nuclear.Model, level::Level)
        basis = level.basis;       npoints = grid.NoPoints
        rhot = zeros( npoints );   wb = zeros( npoints );   wx = zeros( npoints )
        # Compute the charge density of the core orbitals for the given level
        for  sh in basis.subshells
            orb  = basis.orbitals[sh]
            occ  = Basics.computeMeanSubshellOccupation(sh, [level])
            nrho = length(orb.P)
            for    i = 1:nrho   rhot[i] = rhot[i] + occ * (orb.P[i]^2 + orb.Q[i]^2)    end
        end
        # Define the integrant and take the integrals (without alpha)
        for i = 1:npoints
            for j = 1:npoints    rg = max( grid.r[i],  grid.r[j]);    wx[j] = (rhot[j]/rg)   end
            wb[i] = RadialIntegrals.V0(wx, npoints, grid::Radial.Grid)
        end

        # Generate a (local) potential for the given level
        nuclearPotential  = Nuclear.nuclearPotential(nm, grid)

        pot = - ( nuclearPotential.Zr ./ grid.r) + wb

        return( pot )
    end


    """
    Takes the user input radial grid points (R), R*V, Energy, kappa
    """
    function dfree(R0::Vector{Float64}, RV0::Vector{Float64}, E::Float64, k::Int64)
    
        npts=length(R0)   # No of radila points
        r0=zeros(npts+1); r0[2:npts+1]=R0
        rV0=zeros(npts+1); rV0[2:npts+1]=RV0
        phase=Ref{Float64}(0.0)
        delta=Ref{Float64}(0.0)
        r=zeros(Float64,25000)
        p=zeros(Float64,25000)
        q=zeros(Float64,25000)

        ccall((:mydfree, "/home/alok/.julia/packages/JAC/yvKzj/radial/mod_dfree.so"), Cvoid, 
                (Ref{Int64},Ptr{Float64}, Ptr{Float64}, Ref{Float64}, Ref{Int64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ref{Float64}, Ref{Float64}), 
                npts, r0, rV0, E, k, r, p, q, phase, delta)
        
        print("phase = ", phase.x, "  delta = ", delta.x)

        return p[2:npts+1], q[2:npts+1], phase.x, delta.x
    end


end # module
