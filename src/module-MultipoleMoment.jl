
"""
`module  JAC.MultipoleMoment`  
    ... a submodel of JAC that contains all methods for computing (electric) dipole D_z amplitudes 
        between two bound-state levels.
"""
module MultipoleMoment

    using Printf, JAC, ..Basics, ..Defaults, ..ManyElectron, ..Radial, ..InteractionStrength

    """
    `MultipoleMoment.amplitude(mp::EmMultipole, gauge::EmGauge, omega::Float64, finalLevel::Level, initialLevel::Level, 
                                   grid::Radial.Grid; display::Bool=false)`  
        ... to compute the multipole-moment amplitude  
            <alpha_f J_f || Q^(M) (omega, gauge) || alpha_i J_i>  for the given final and initial level. All the amplitudes are 
            calculated by means of Johnson's (2007) multipole-transition amplitudes from which also all constrains w.r.t the 
            allowed gauges are inherited. A value::ComplexF64 is returned.
    """
    function amplitude(mp::EmMultipole, gauge::EmGauge, omega::Float64, finalLevel::Level, initialLevel::Level, 
                       grid::Radial.Grid; display::Bool=false)

        amplitude = MultipoleMoment.transitionAmplitude(mp, gauge, omega, finalLevel, initialLevel, grid; display=false)
        #
        if  display
            sa = @sprintf("%.5e", amplitude.re) * "  " * @sprintf("%.5e", amplitude.im)
            println("    < level=$(finalLevel.index) [J=$(finalLevel.J)$(string(finalLevel.parity))] || Q^($mp) ($omega a.u., $gauge) ||" *
                    " $(initialLevel.index) [$(initialLevel.J)$(string(initialLevel.parity))] >  = " * sa)
            printSummary, iostream = Defaults.getDefaults("summary flag/stream")
            if  printSummary
                println(iostream, "    Multipole-moment amplitude:     " *
                                  " < level=$(finalLevel.index) [J=$(finalLevel.J)$(string(finalLevel.parity))] || Q^($mp) ($omega a.u., $gauge) ||" *
                                  " $(initialLevel.index) [$(initialLevel.J)$(string(initialLevel.parity))] >  = " * sa)
            end
        end
        
        return( amplitude )
    end


    """
    `MultipoleMoment.dipoleAmplitude(finalLevel::Level, initialLevel::Level, grid::Radial.Grid; display::Bool=false)`  
         ... to compute the dipole amplitude   <(alpha_f J_f, kappa) J_i || D || alpha_i J_i>  for the given final 
            and initial level. A value::ComplexF64 is returned.
    """
    function dipoleAmplitude(finalLevel::Level, initialLevel::Level, grid::Radial.Grid; display::Bool=false)
        
        if     finalLevel.parity == initialLevel.parity   
            amplitude = 0.0 + 0.0 * im
        else
            nf = length(finalLevel.basis.csfs);    ni = length(initialLevel.basis.csfs)
            printstyled("Compute dipole matrix of dimension $nf x $ni in the final- and initial-state bases " *
                        "for the transition [$(initialLevel.index)- $(finalLevel.index)] ... ", color=:light_green)
            matrix = zeros(ComplexF64, nf, ni)
            #
            for  r = 1:nf
                for  s = 1:ni
                    ##x wa = Basics.compute("angular coefficients: 1-p, Grasp92", 0, 1, finalLevel.basis.csfs[r], initialLevel.basis.csfs[s])
                    # Calculate the spin-angular coefficients
                    if  Defaults.saRatip()
                        waR = Basics.compute("angular coefficients: 1-p, Grasp92", 0, 1, finalLevel.basis.csfs[r], initialLevel.basis.csfs[s])
                        wa  = waR       
                    end
                    if  Defaults.saGG()
                        subshellList = initialLevel.basis.subshells
                        opa = SpinAngular.OneParticleOperator(1, plus, true)
                        waG = SpinAngular.computeCoefficients(opa, finalLevel.basis.csfs[r], initialLevel.basis.csfs[s], subshellList) 
                        wa  = waG
                    end
                    if  Defaults.saRatip() && Defaults.saGG() && true
                        if  length(waR) != 0     println("\n>> Angular coeffients from GRASP/MCT   = $waR ")    end
                        if  length(waG) != 0     println(  ">> Angular coeffients from SpinAngular = $waG ")    end
                    end
                    #
                    for  coeff in wa
                        ja = Basics.subshell_2j(finalLevel.basis.orbitals[coeff.a].subshell)
                        jb = Basics.subshell_2j(initialLevel.basis.orbitals[coeff.b].subshell)
                        tamp  = JAC.InteractionStrength.dipole(finalLevel.basis.orbitals[coeff.a], initialLevel.basis.orbitals[coeff.b], grid)
                        matrix[r,s] = matrix[r,s] + coeff.T * tamp  
                    end
                end
            end
            printstyled("done. \n", color=:light_green)
            amplitude::ComplexF64 = transpose(finalLevel.mc) * matrix * initialLevel.mc 
        end
        #
        if  display
            sa = @sprintf("%.5e", amplitude.re) * "  " * @sprintf("%.5e", amplitude.im)
            println("    < level=$(finalLevel.index) [J=$(finalLevel.J)$(string(finalLevel.parity))] || D ||" *
                    " $(initialLevel.index) [$(initialLevel.J)$(string(initialLevel.parity))] >  = " * sa)
            printSummary, iostream = Defaults.getDefaults("summary flag/stream")
            if  printSummary
                println(iostream, "    Dipole amplitude:  < level=$(finalLevel.index) [J=$(finalLevel.J)$(string(finalLevel.parity))] || D ||" *
                                  " $(initialLevel.index) [$(initialLevel.J)$(string(initialLevel.parity))] >  = " * sa)
            end
        end
        
        return( amplitude )
    end


    """
    `JAC.MultipoleMoment.transitionAmplitude(mp::EmMultipole, gauge::EmGauge, omega::Float64, finalLevel::Level, initialLevel::Level, 
                                             grid::Radial.Grid; display::Bool=false)`  
        ... to compute Johnson's (2007) multipole-transition amplitude  <alpha_f J_f || T^(M, absorption) (omega, gauge) || alpha_i J_i>  
            for the given final and initial level. Allowed gauges are {Velocity, Length} for the electric-multipole transition 
            amplitudes and Magnetic for the magnetic-multipole transition amplitudes, respectively. A value::ComplexF64 is 
            returned.
    """
    function transitionAmplitude(mp::EmMultipole, gauge::EmGauge, omega::Float64, finalLevel::Level, initialLevel::Level, 
                                 grid::Radial.Grid; display::Bool=false)
        nf = length(finalLevel.basis.csfs);    ni = length(initialLevel.basis.csfs)
        printstyled("Compute multipole-moment transition matrix of dimension $nf x $ni in the final- and initial-state bases " *
                    "for the transition [$(initialLevel.index)- $(finalLevel.index)] ... ", color=:light_green)
        matrix = zeros(ComplexF64, nf, ni)
        #
        for  r = 1:nf
            for  s = 1:ni
                wa = compute("angular coefficients: 1-p, Grasp92", 0, mp.L, finalLevel.basis.csfs[r], initialLevel.basis.csfs[s])
                for  coeff in wa
                    ##x ja = Basics.subshell_2j(finalLevel.basis.orbitals[coeff.a].subshell)
                    ##x jb = Basics.subshell_2j(initialLevel.basis.orbitals[coeff.b].subshell)
                    tamp  = JAC.InteractionStrength.multipoleTransition(mp, gauge, omega, finalLevel.basis.orbitals[coeff.a], 
                                                                        initialLevel.basis.orbitals[coeff.b], grid)
                    matrix[r,s] = matrix[r,s] + coeff.T * tamp  
                end
            end
        end
        printstyled("done. \n", color=:light_green)
        amplitude = transpose(finalLevel.mc) * matrix * initialLevel.mc 
        #
        if  display
            sa = @sprintf("%.5e", amplitude.re) * "  " * @sprintf("%.5e", amplitude.im)
            println("    < level=$(finalLevel.index) [J=$(finalLevel.J)$(string(finalLevel.parity))] ||" *
                    " T^($mp, absorption) ($omega a.u., $gauge) ||" *
                    " $(initialLevel.index) [$(initialLevel.J)$(string(initialLevel.parity))] >  = " * sa)
            printSummary, iostream = Defaults.getDefaults("summary flag/stream")
            if  printSummary
                println(iostream, "    Multipole-transition amplitude:      < level=$(finalLevel.index) [J=$(finalLevel.J)$(string(finalLevel.parity))] ||" *
                                  " T^($mp, absorption) ($omega a.u., $gauge) ||" *
                                  " $(initialLevel.index) [$(initialLevel.J)$(string(initialLevel.parity))] >  = " * sa)
            end
        end
        
        return( amplitude )
    end

end # module
