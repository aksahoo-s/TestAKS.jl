
"""
`module JAC.LSjj`  
    ... a submodel of JAC that contains methods and (numerical) values for performing the jj-LS transformation of atomic
        levels; this transformation is mainly based on global data lists which are only accessible within this module.
"""
module  LSjj

    using  Printf,  ..AngularMomentum, ..Basics, ..Defaults, ..ManyElectron, ..TableStrings
    export CsfNR, BasisNR, LevelNR, MultipletNR

    
    abstract type  AbstractOpenShell                        end
    struct         ZeroOpenShell    <:  AbstractOpenShell   end
    struct         OneOpenShell     <:  AbstractOpenShell   end
    struct         TwoOpenShells    <:  AbstractOpenShell   end
    struct         ThreeOpenShells  <:  AbstractOpenShell   end
    
    
    """
    `struct   LSjj.ShellStateNR`  
        ... defines a struct for the nonrelativistic antisymmetric shell states in the seniority scheme; this struct can be accessed 
            only internally and, therefore, only the standard constructor is supported.

        + shell        ::Shell           ... to refer to a particular shell. 
        + occ          ::Int64           ... occupation of the shell
        + w            ::Int64           ... w quantum number of shell state
        + QQ           ::Int64           ... 2*Q of the shell state
        + LL           ::Int64           ... 2*L of the shell state
        + SS           ::Int64           ... 2*S of the shell state
    """
    struct  ShellStateNR
        shell        ::Shell
        occ          ::Int64
        w            ::Int64
        QQ           ::Int64
        LL           ::Int64
        SS           ::Int64
    end     


    # `Base.show(io::IO, state::ShellStateNR)`  ... prepares a proper printout of the variable state::ShellStateNR.
    function Base.show(io::IO, state::ShellStateNR) 
        print(io, string(state) )
    end


    # `Base.string(state::ShellStateNR)`  ... provides a proper printout of the variable state::ShellStateNR.
    function Base.string(state::ShellStateNR)
        sa = "ShellState: [" * string(state.shell) * "^$(state.occ), w=$(state.w), 2*Q=$(state.QQ), 2*L=$(state.LL), 2*S=$(state.SS)]"
    end


    """
    `LSjj.shellStateString(shell::String, occ::Int64, w::Int64, Q::AngularJ64, 
                           L::AngularJ64, S::AngularJ64, LX::AngularJ64, SX::AngularJ64)`  
        ... to provide a string of a given shell state in the form !!To be worked out !! '[2p_1/2^occ]_(seniorityNr, J_sub), X=Xo' ... .
    """
    function shellStateString(shell::String, occ::Int64, w::Int64, Q::AngularJ64, L::AngularJ64, S::AngularJ64, LX::AngularJ64, SX::AngularJ64)
        error("stop a: To be adapted yet. ")
        sa = "[" * subshell * "^$occ]_($seniorityNr, " * string(Jsub) * ") X=" * string(X)
        return( sa )
    end

    
    """
    `struct  CsfNR`  
        ... defines a type for a nonrelativistic configuration state function (CSF) in terms of its orbitals (sequence of 
            orbitals), their occupation as well as the angular momenta, seniority and the coupling of the corresponding 
            antisymmetric subshell states.
    
        + J                ::AngularJ64           ... Total angular momentum J.
        + parity           ::Parity               ... Total parity.
        + occupation       ::Array{Int64,1}       ... Occupation of the orbitals with regard to the specified shell list.
        + w                ::Array{Int64,1}       ... List of w quantum numbers of the antisymmetric subshell states.
        + QQ               ::Array{Int64,1}       ... List of total quasispin 2*Q, likely = seniority values of the 
                                                      antisymmetric subshell states.
        + shellL           ::Array{AngularJ64,1}  ... List of L-values of the antisymmetric subshell states.
        + shellS           ::Array{AngularJ64,1}  ... List of S-values of the antisymmetric subshell states.
        + shellLX          ::Array{AngularJ64,1}  ... intermediate L-values in the coupling of the antisymmetric shell states.
        + shellSX          ::Array{AngularJ64,1}  ... intermediate S-values in the coupling of the antisymmetric shell states.
    """ 
    struct  CsfNR
        J                  ::AngularJ64
        parity             ::Parity
        occupation         ::Array{Int64,1}
        w                  ::Array{Int64,1} 
        QQ                 ::Array{Int64,1}
        shellL             ::Array{AngularJ64,1}
        shellS             ::Array{AngularJ64,1}
        shellLX            ::Array{AngularJ64,1} 
        shellSX            ::Array{AngularJ64,1}
    end
    
    
    """
    `Base.show(io::IO, csf::CsfNR)`  ... prepares a proper printout of the variable csf::CsfNR.
    """
    function Base.show(io::IO, csf::CsfNR) 
        println(io, "J:             $(csf.J)  ")
        println(io, "parity:        $(csf.parity)  ")
        println(io, "occupation:    $(csf.occupation)  ")
        println(io, "w:             $(csf.w)  ")
        println(io, "QQ:            $(csf.QQ)  ")
        println(io, "shellL:        $(csf.shellL)  ")
        println(io, "shellS:        $(csf.shellS)  ")
        println(io, "shellLX:       $(csf.shellLX)  ")
        println(io, "shellSX:       $(csf.shellSX)  ")
    end
    

    """
    `Base.:(==)(csfa::CsfNR, csfb::CsfNR)`  
        ... compares (recursively) two nonrelativistic CSFs and return true if all subfields are equal, and false otherwise.
    """
    function  Base.:(==)(csfa::CsfNR, csfb::CsfNR)
        # It is assumed that both csf's are defined with regard to the same shell list.
    
        if  length(csfa.occupation) != length(csfb.occupation)      return( false )    end
        if  csfa.J  !=  csfb.J  ||  csfa.parity  !=  csfb.parity	return( false )    end
    	if  csfa.occupation      !=  csfb.occupation			    return( false )    end
        if  csfa.w               !=  csfb.w			                return( false )    end
        if  csfa.QQ              !=  csfb.QQ			            return( false )    end
        if  csfa.shellL          !=  csfb.shellL			        return( false )    end
        if  csfa.shellS          !=  csfb.shellS			        return( false )    end
        if  csfa.shellLX         !=  csfb.shellLX			        return( false )    end
        if  csfa.shellSX         !=  csfb.shellSX			        return( false )    end
    	
        return( true )
    end
    
    
    """
    `struct  LSjj.BasisNR`  
        ... defines a type for a nonrelativistic atomic basis, including the (full) specification of the configuration space 
            but without radial orbitals.

        + NoElectrons    ::Int64              ... No. of electrons.
        + shells         ::Array{Shell,1}     ... Explicitly given shell list for this basis.
        + csfs           ::Array{CsfNR,1}     ... List of nonrelativistic LS-coupled CSF.	    
    """
    struct  BasisNR
        NoElectrons	     ::Int64 
        shells	         ::Array{Shell,1}
        csfs             ::Array{CsfNR,1}	      
    end 
    
    
    """
    `LSjj.BasisNR()`  ... constructor for an 'empty' instance::BasisNR.
    """
    function BasisNR()
    	  BasisNR(0, Shell[], CsfNR[] )
    end
    
    
    """
    `Base.show(io::IO, basis::BasisNR)`  ... prepares a proper printout of the variable basis::BasisNR.
    """
    function Base.show(io::IO, basis::BasisNR) 
        println(io, "NoElectrons:        $(basis.NoElectrons)  ")
        println(io, "shells:             $(basis.shells)  ")
        println(io, "csfs:               $(basis.csfs)  ")
    end
    

    """
    `struct  LSjj.LevelNR`  
        ... defines a type for an atomic level in terms of its quantum number, energy and with regard to a nonrelativistic basis.
    
        + J            ::AngularJ64       ... Total angular momentum J.
        + parity       ::Parity           ... Parity of the level.
        + index        ::Int64            ... Index of this level in its original multiplet, or 0.
        + energy       ::Float64          ... energy
        + basis        ::BasisNR          ... nonrelativistic atomic basis
        + mc           ::Vector{Float64}  ... Vector of mixing coefficients w.r.t basis.
    """
    struct  LevelNR
        J              ::AngularJ64
        parity         ::Parity
        index          ::Int64
        energy         ::Float64
        basis          ::BasisNR
        mc             ::Vector{Float64}
    end 
    
    
    # `Base.show(io::IO, level::LevelNR)`  ... prepares a proper printout of the variable level::LevelNR.
    function Base.show(io::IO, level::LevelNR) 
        println(io, "Level: L = $(level.L), S = $(level.S), J = $(level.J), parity = $(level.parity), index = $(level.index) ")
        println(io, "energy:         $(level.energy)  ")
        println(io, "basis:           (level.basis)  ")
        println(io, "mc:             $(level.mc)  ")
    end
    
    
    """
    `struct  LSjj.MultipletNR`  
        ... defines a type for an ordered list of atomic levelsin a nonrelativistic basis.
    
        + name       ::String	         ... A name associated to the multiplet.
        + levels     ::Array{LevelNR,1}  ... A list of levels (pointers).
    """
    struct  MultipletNR
        name         ::String
        levels       ::Array{LevelNR,1}
    end 
    
    
    # `Base.show(io::IO, multiplet::MultipletNR)`  ... prepares a proper printout of the variable multiplet::MultipletNR.
    function Base.show(io::IO, multiplet::MultipletNR) 
        println(io, "name:        $(multiplet.name)  ")
        println(io, "levels:      $(multiplet.levels)  ")
    end


    """
    `LSjj.expandLevelsIntoLS(multiplet::Multiplet, settings::ManyElectron.LSjjSettings)  
        ... expand and print all (selected) levels from the multiplet into their LS representation. The request is controlled
            by settings.makeIt  and the details of this expansion by further parameters given in these settings.
            nothing is returned.
    """
    function expandLevelsIntoLS(multiplet::Multiplet, settings::ManyElectron.LSjjSettings)
        # Return immediately if not expansion is to be made.
        if  !(settings.makeIt)   return( nothing )    end
        
        # Determine all selected levels by their indiced; at present, simply all levels are transformed at present
        indexList = Int64[]
        for  levelR in multiplet.levels    push!(indexList, levelR.index)     end
        
        # Make the jj-LS expansion of all selected levels
        if  Basics.isStandardSubshellList(multiplet.levels[1].basis)
            shellList = Basics.extractNonrelativisticShellList(multiplet.levels[1].basis.subshells) 
            confList  = Basics.extractNonrelativisticConfigurations(multiplet.levels[1].basis)
            csfsNR    = LSjj.generateNonrelativisticCsfList(confList, shellList)
            basisNR   = BasisNR(multiplet.levels[1].basis.NoElectrons, shellList, csfsNR)
            ncsfs     = length(basisNR.csfs)
            printstyled("\n  LSjj.expandLevelsIntoLS():: Relativistic basis with $(length(multiplet.levels[1].basis.csfs)) CsfR " *
                        "will be transformed into a nonrelativistic basis with $ncsfs CsfNR. \n", color=:light_green)
                    
            # Create a dictionary of (empty) eigenvectors for the selected levels
            mcVectors = Dict{Int64, Array{Float64,1}}()
            for index in indexList    mcVectors = Base.merge( mcVectors, Dict(index => zeros(ncsfs)) )    end
            
            # Expand in terms, if requested, all relativistic CSF into the nonrelativistic basis
            for  r = 1:length(multiplet.levels[1].basis.csfs)
                csfR = multiplet.levels[1].basis.csfs[r]
                # Cycle over this CSF if it has too low weight (not yet)
                # Determine the number of open shells in CsfR
                conf       = Basics.extractNonrelativisticConfigurationFromCsfR(csfR, multiplet.levels[1].basis)
                openShells = Basics.extractNoOpenShells(conf)
                mcCsfR     = LSjj.expandCsfRintoNonrelativisticBasis(openShells, csfR, multiplet.levels[1].basis, basisNR)
                # Cycle through all selected levels    
                for  levelR  in  multiplet.levels  
                    if  levelR.index  in  indexList
                        mcLevel = levelR.mc[r]
                        mcVectors[levelR.index] = mcVectors[levelR.index] + mcLevel * mcCsfR
                    end
                end
            end
            
            # Use the dictionary of expansion coefficients to define the nonrelativistic levels and multiplet
            levelsNR = LevelNR[]
            for  levelR  in  multiplet.levels  
                if  levelR.index  in  indexList
                    levelNR = LevelNR( levelR.J, levelR.parity, levelR.index, levelR.energy, basisNR, mcVectors[levelR.index])
                    push!( levelsNR, levelNR)
                end
            end
            multipletNR = MultipletNR( "LS-expanded " * multiplet.name, levelsNR)
            
            # Check the consistency of the LS expansion of the selected levels
            LSjj.checkLSexpansion(multipletNR)
            
            # Print the results to screen and elsewhere
            LSjj.displayLSexpansion(stdout, multiplet, multipletNR)
            printSummary, iostream = Defaults.getDefaults("summary flag/stream")
            if  printSummary   LSjj.displayLSexpansion(iostream, multiplet, multipletNR)    end
            
        else  
            @warn "LSjj.expandLevelsIntoLS():: Inappropriate basis without a standard list of subshells;" *
                "no jj-LS epxansion is supported in this case."
        end
        
        return( nothing )
    end


    """
    `LSjj.checkLSexpansion(multipletNR::MultipletNR)`
        ... to check the consistency and normalization of the transformed eigenvectors; nothing is returned.
    """
    function checkLSexpansion(multipletNR::MultipletNR) 
        basis = multipletNR.levels[1].basis
        for  level  in  multipletNR.levels
            weight = 0.
            for  r = 1:length(basis.csfs)    
                weight = weight + level.mc[r]*level.mc[r]
                if  level.mc[r]^2 > 0.  &&  level.J != basis.csfs[r].J   
                    println("Level $(level.index):: Total angular momenta $(level.J)  !=  $(basis.csfs[r].J) for r = $r")
                end
            end
            #
            if  weight < 0.99   println("Level $(level.index):: Total weight $weight  < 1.0 ")      end
        end
        return( nothing )
    end



    """
    `LSjj.displayLSexpansion(stream::IO, multiplet::Multiplet, multipletNR::MultipletNR)`
        ... to display the LS expansion of the selected (relativistic) levels from multiplet in a neat format;
            nothing is returned.
    """
    function displayLSexpansion(stream::IO, multiplet::Multiplet, multipletNR::MultipletNR)
        
        println(stream, " ")
        println(stream, "  LS-expansion of selected atomic levels:")
        println(stream, " ")
        
        # Determine all selected levels by their indiced; at present, simply all levels are transformed at present
        indexList = Int64[]
        for  levelR in multiplet.levels    push!(indexList, levelR.index)     end
        
        for  levelR  in  multiplet.levels
            if  levelR.index in indexList
                # First, find the levelR in multipletNR.levels
                for  levelNR in multipletNR.levels
                    if  levelNR.index == levelR.index
                        # Extract the weights of all (nonrelativistic) configurations for the selected level
                        weights = LSjj.extractConfigurationWeightsOfLevelNR(levelNR)
                        # Find the configuration with the largest weight
                        wa = 0.;    conf = Configuration("1s")
                        for (k,v) in  weights
                            if  v > wa      wa = v;     conf = k    end
                        end
                        # Print the level details and the weight of the 'most important' nonrelativistic configuration
                        println(stream, "    Level   $(levelNR.index)      $(LevelSymmetry(levelNR.J,levelNR.parity))   " *
                                        @sprintf("%.8e", Defaults.convertUnits("energy: from atomic", levelNR.energy))  * "  " *
                                        TableStrings.inUnits("energy") * "  has weight " * @sprintf("%.4e", wa) * " of $conf" )
                        # Next, print also the contribution of the major nonrelativistic csf
                        wb = abs.(levelNR.mc);   wb = wb.^2;   idx = sortperm(wb, rev = true)
                        for  ix = 1:min(3, length(idx))
                            wx = wb[ idx[ix] ]
                            if   wx > 1.0e-4    println(stream, "       $(ix))   " * @sprintf("%.4e", wx) * "   of  " *
                                                        LSjj.shortString(levelNR.basis.csfs[ idx[ix] ], levelNR.basis) )    end
                        end
                        println("")
                    end
                end
            end
        end
        
        return( nothing )
    end


    """
    `LSjj.expandCsfRintoNonrelativisticBasis(openShells::ZeroOpenShell, csfR::CsfR, basisR::Basis, basisNR::BasisNR)`
        ... to expand a relativistic CsfR without any nonrelativistic open shell into the given nonrelativistic basis; 
            a mcVector::Array{Float64,1} is returned with the same lengths as length(basisNR.csfs). It is checked that the 
            expansion is complete, and a warning is issued if this is not the case.
    """
    function expandCsfRintoNonrelativisticBasis(openShells::ZeroOpenShell, csfR::CsfR, basisR::Basis, basisNR::BasisNR)
        if  csfR.J != AngularJ64(0)  ||   csfR.parity != Basics.plus    error("stop a")   end
        mcVector = Float64[]
        confR     = Basics.extractNonrelativisticConfigurationFromCsfR(csfR,  basisR)
        for  csf in basisNR.csfs
            confNR = LSjj.extractConfigurationFromCsfNR(csf, basisNR)   
            if      confR == confNR    push!( mcVector, 1.0)
            else                       push!( mcVector, 0.0)
            end
        end
        return( mcVector )
    end


    """
    `LSjj.expandCsfRintoNonrelativisticBasis(openShells::OneOpenShell, csfR::CsfR, basisR::Basis, basisNR::BasisNR)`
        ... to expand a relativistic CsfR with one nonrelativistic open shell into the given nonrelativistic basis; 
            a mcVector::Array{Float64,1} is returned with the same lengths as length(basisNR.csfs). It is checked that the 
            expansion is complete, and a warning is issued if this is not the case.
    """
    function expandCsfRintoNonrelativisticBasis(openShells::OneOpenShell, csfR::CsfR, basisR::Basis, basisNR::BasisNR)
        mcVector      = Float64[]
        # Find the open subshells and determine all the quantum numbers
        rQN           = Basics.extractOpenShellQNfromCsfR(csfR, basisR)
        rShells, rOCC = Basics.extractShellOccupationFromCsfR(csfR, basisR)
        
        if  length( keys(rQN) ) != 1    error("stop a")     end
        if  rShells != basisNR.shells   error("stop b")     end
        
        for  (rsh, rv)  in  rQN
            rQNm = rv[1];   rQNp = rv[2]    # Get QN of the relativistic subshells
            for  s = 1:length(basisNR.csfs)
                csfNR = basisNR.csfs[s]
                if  csfR.J != csfNR.J   
                    push!( mcVector, 0.)
                else
                    me = 0.
                    if  rOCC == csfNR.occupation
                        nrQN = LSjj.extractOpenShellQNfromCsfNR(csfNR, basisNR)
                        if  length( keys(nrQN) ) != 1    error("stop c")     end
                        for  (nrsh, nrv)  in  nrQN
                            if  rsh != nrsh              error("stop d")     end
                            # Determine the jj-LS overlap coefficient from the tabulation; no re-coupling coefficients need to be
                            # considered for a single nonrelativistic open shell
                            if  rsh.l == 0  twojm = 1    else    twojm = 2*rsh.l - 1  end
                            twojp = 2*rsh.l + 1 
                            QQm = Int64( (twojm+1)/2 - rQNm[3]);     QQp = Int64( (twojp+1)/2 - rQNp[3])
                            me  = LSjj.getLSjjCoefficient(rsh.l, nrv[2], 
                                    LS_jj_qn(nrv[3], nrv[4], nrv[5], nrv[6], rQNp[5], rQNm[2], QQm, rQNm[4], QQp, rQNp[4]) )
                        end
                    end
                    push!( mcVector, me)
                end
            end
        end
        #
        return( mcVector )
    end


    """
    `LSjj.expandCsfRintoNonrelativisticBasis(openShells::TwoOpenShells, csfR::CsfR, basisR::Basis, basisNR::BasisNR)`
        ... to expand a relativistic CsfR with two nonrelativistic open shells into the given nonrelativistic basis; 
            a mcVector::Array{Float64,1} is returned with the same lengths as length(basisNR.csfs). It is checked that the 
            expansion is complete, and a warning is issued if this is not the case.
    """
    function expandCsfRintoNonrelativisticBasis(openShells::TwoOpenShells, csfR::CsfR, basisR::Basis, basisNR::BasisNR)
        mcVector      = Float64[]
        # Find the open subshells and determine all the quantum numbers
        rQN           = Basics.extractOpenShellQNfromCsfR(csfR, basisR)
        rShells, rOCC = Basics.extractShellOccupationFromCsfR(csfR, basisR)
        rKeys         = keys(rQN)
        
        if  length(rKeys) != 2          error("stop a")     end
        if  rShells != basisNR.shells   error("stop b")     end
        
        # Initialize all quantum numbers that occur in loops below
        rsh1 = rsh2 = nrsh1 = nrsh2 = Shell("9m")   
        Nm1 = JJm1 = JJp1 = XXp1 = QQm1 = QQp1 = Nm2 = JJm2 = JJp2 = XXm2 = QQm2 = QQp2 = 0  
        Jm1 = Jp1 = Xp1 = Jm2 = Jp2 = Xm2 = L1 = S1 = L2 = S2 = L12 = S12 = J = AngularJ64(0)
        N1 = w1 = QQ1 = LL1 = SS1 = N2 = w2 = QQ2 = LL2 = SS2 =  0 
        
        # Set open-shell QN for the given csfR
        ns = 0
        for  rsh in rShells
            if  rsh in rKeys
                ns = ns + 1;    rv = rQN[rsh];    rQNm = rv[1];    rQNp = rv[2]
                twojp = 2*rsh.l + 1;        if  rsh.l == 0  twojm = 1    else    twojm = 2*rsh.l - 1  end
                #
                if      ns == 1  Nm1 = rQNm[2];    JJm1 = rQNm[4];       JJp1 = rQNp[4];         XXp1 = rQNp[5];    rsh1 = rsh
                                QQm1 = Int64( (twojm+1)/2 - rQNm[3]);   QQp1 = Int64( (twojp+1)/2 - rQNp[3])
                                Jm1  = AngularJ64(JJm1//2);             Jp1  = AngularJ64(JJp1//2);    Xp1  = AngularJ64(XXp1//2)
                elseif  ns == 2  Nm2 = rQNm[2];    JJm2 = rQNm[4];       JJp2 = rQNp[4];         XXm2 = rQNm[5];    rsh2 = rsh
                                QQm2 = Int64( (twojm+1)/2 - rQNm[3]);   QQp2 = Int64( (twojp+1)/2 - rQNp[3])
                                Jm2  = AngularJ64(JJm2//2);             Jp2  = AngularJ64(JJp2//2);    Xm2  = AngularJ64(XXm2//2)
                else    error("stop c")     
                end
            end
        end
        
        # Now cycle over all nonrelativistic csfNR in the given basis
        for  s = 1:length(basisNR.csfs)
            csfNR = basisNR.csfs[s]
            if  csfR.J != csfNR.J   
                push!( mcVector, 0.)
            else
                me = 0.
                if  rOCC == csfNR.occupation
                    nrQN   = LSjj.extractOpenShellQNfromCsfNR(csfNR, basisNR)
                    nrKeys = keys(nrQN)
                    # Set open-shell QN for the given csfNR in the loop
                    ns = 0
                    for  rsh in rShells
                        if  rsh in nrKeys
                            ns = ns + 1;    nrv = nrQN[rsh]
                            if      ns == 1     N1 = nrv[2];     w1 = nrv[3];   QQ1 = nrv[4];    LL1 = nrv[5];    SS1 = nrv[6];    nrsh1 = rsh
                                                L1 = AngularJ64(LL1//2);        S1  = AngularJ64(SS1//2)
                            elseif  ns == 2     N2 = nrv[2];     w2 = nrv[3];   QQ2 = nrv[4];    LL2 = nrv[5];    SS2 = nrv[6];    nrsh2 = rsh
                                                L2 = AngularJ64(LL2//2);        S2  = AngularJ64(SS2//2)
                                                L12 = AngularJ64(nrv[7]//2);    S12 = AngularJ64(nrv[8]//2);      J = csfR.J
                            else    error("stop d")     
                            end
                        end
                    end
                    #
                    if  rsh1 != nrsh1   ||   rsh2 != nrsh2      error("stop e")     end
                    
                    # Cycle over all intermediate angular momenta T1, T2
                    T1List = oplus(L1, S1);     T2List = oplus(L2, S2)
                    for  T1 in T1List
                        for  T2 in T2List
                            wa1 = LSjj.getLSjjCoefficient(rsh1.l, N1, LS_jj_qn(w1, QQ1, LL1, SS1, Basics.twice(T1), Nm1, QQm1, JJm1, QQp1, JJp1) )
                            wa2 = LSjj.getLSjjCoefficient(rsh2.l, N2, LS_jj_qn(w2, QQ2, LL2, SS2, Basics.twice(T2), Nm2, QQm2, JJm2, QQp2, JJp2) )
                            if AngularMomentum.triangularDelta(T1,Jm1,Jp1) == 1  &&  T1 == Xp1
                                me  = me + sqrt(AngularMomentum.bracket([T1, T2, L12, S12])) * AngularMomentum.Wigner_9j(L1,S1,T1, L2,S2,T2, L12,S12,J) *
                                        ## phase from Valeria (2020) which seem to have no effect upon the expansion
                                        ## AngularMomentum.phaseFactor([L1, 1, S1, -1, T1,  1, L2, 1, S2, -1, T2, 1, L12, 1, S12, -1, J]) *   ## phase from Valeria (2020)
                                        AngularMomentum.phaseFactor([Jm2, 1, Jp2, 1, T1, 1, J]) * 
                                        sqrt(AngularMomentum.bracket([T2, Xm2])) * AngularMomentum.Wigner_6j(T1, Jm2, Xm2, Jp2, J, T2) * wa1 * wa2
                            end
                        end
                    end
                end
                push!( mcVector, me)
            end
        end
        #
        return( mcVector )
    end


    """
    `LSjj.expandCsfRintoNonrelativisticBasis(openShells::ThreeOpenShells, csfR::CsfR, basisR::Basis, basisNR::BasisNR)`
        ... to expand a relativistic CsfR with three nonrelativistic open shells into the given nonrelativistic basis; 
            a mcVector::Array{Float64,1} is returned with the same lengths as length(basisNR.csfs). It is checked that the 
            expansion is complete, and a warning is issued if this is not the case.
    """
    function expandCsfRintoNonrelativisticBasis(openShells::ThreeOpenShells, csfR::CsfR, basisR::Basis, basisNR::BasisNR)
        mcVector      = Float64[]
        # Find the open subshells and determine all the quantum numbers
        rQN           = Basics.extractOpenShellQNfromCsfR(csfR, basisR)
        rShells, rOCC = Basics.extractShellOccupationFromCsfR(csfR, basisR)
        rKeys         = keys(rQN)
        
        if  length(rKeys) != 3          error("stop a")     end
        if  rShells != basisNR.shells   error("stop b")     end
        
        # Initialize all quantum numbers that occur in loops below
        rsh1 = rsh2 = rsh3 = nrsh1 = nrsh2 = nrsh3 = Shell("9m")   
        Nm1 = JJm1 = JJp1 = XXp1 = QQm1 = QQp1 = Nm2 = JJm2 = JJp2 = XXm2 = QQm2 = QQp2 = Nm3 = JJm3 = JJp3 = XXm3 = QQm3 = QQp3 = 0  
        Jm1 = Jp1 = Xp1 = Jm2 = Jp2 = Xm2 = Jm3 = Jp3 = Xm3 = L1 = S1 = L2 = S2 = L3 = S3 = L12 = S12 = T12 = L123 = S123 = J = AngularJ64(0)
        N1 = w1 = QQ1 = LL1 = SS1 = N2 = w2 = QQ2 = LL2 = SS2 = N3 = w3 = QQ3 = LL3 = SS3 =  0 
        
        # Set open-shell QN for the given csfR
        ns = 0
        for  rsh in rShells
            if  rsh in rKeys
                ns = ns + 1;    rv = rQN[rsh];    rQNm = rv[1];    rQNp = rv[2]
                twojp = 2*rsh.l + 1;        if  rsh.l == 0  twojm = 1    else    twojm = 2*rsh.l - 1  end
                #
                if      ns == 1  Nm1 = rQNm[2];    JJm1 = rQNm[4];       JJp1 = rQNp[4];         XXp1 = rQNp[5];    rsh1 = rsh
                                QQm1 = Int64( (twojm+1)/2 - rQNm[3]);   QQp1 = Int64( (twojp+1)/2 - rQNp[3])
                                Jm1  = AngularJ64(JJm1//2);             Jp1  = AngularJ64(JJp1//2);    Xp1  = AngularJ64(XXp1//2)
                elseif  ns == 2  Nm2 = rQNm[2];    JJm2 = rQNm[4];       JJp2 = rQNp[4];         XXm2 = rQNm[5];    rsh2 = rsh
                                QQm2 = Int64( (twojm+1)/2 - rQNm[3]);   QQp2 = Int64( (twojp+1)/2 - rQNp[3])
                                Jm2  = AngularJ64(JJm2//2);             Jp2  = AngularJ64(JJp2//2);    Xm2  = AngularJ64(XXm2//2)
                elseif  ns == 3  Nm3 = rQNm[2];    JJm3 = rQNm[4];       JJp3 = rQNp[4];         XXm3 = rQNm[5];    rsh3 = rsh
                                QQm3 = Int64( (twojm+1)/2 - rQNm[3]);   QQp3 = Int64( (twojp+1)/2 - rQNp[3])
                                Jm3  = AngularJ64(JJm3//2);             Jp3  = AngularJ64(JJp3//2);    Xm3  = AngularJ64(XXm3//2)
                else    error("stop c")     
                end
            end
        end
        
        # Now cycle over all nonrelativistic csfNR in the given basis
        for  s = 1:length(basisNR.csfs)
            csfNR = basisNR.csfs[s]
            if  csfR.J != csfNR.J   
                push!( mcVector, 0.)
            else
                me = 0.
                if  rOCC == csfNR.occupation
                    nrQN   = LSjj.extractOpenShellQNfromCsfNR(csfNR, basisNR)
                    nrKeys = keys(nrQN)
                    # Set open-shell QN for the given csfNR in the loop
                    ns = 0
                    for  rsh in rShells
                        if  rsh in nrKeys
                            ns = ns + 1;    nrv = nrQN[rsh]
                            if      ns == 1     N1 = nrv[2];     w1 = nrv[3];   QQ1 = nrv[4];    LL1 = nrv[5];    SS1 = nrv[6];    nrsh1 = rsh
                                                L1 = AngularJ64(LL1//2);        S1  = AngularJ64(SS1//2)
                            elseif  ns == 2     N2 = nrv[2];     w2 = nrv[3];   QQ2 = nrv[4];    LL2 = nrv[5];    SS2 = nrv[6];    nrsh2 = rsh
                                                L2 = AngularJ64(LL2//2);        S2  = AngularJ64(SS2//2)
                                                L12 = AngularJ64(nrv[7]//2);    S12 = AngularJ64(nrv[8]//2);
                            elseif  ns == 3     N3 = nrv[2];     w3 = nrv[3];   QQ3 = nrv[4];    LL3 = nrv[5];    SS3 = nrv[6];    nrsh3 = rsh
                                                L3 = AngularJ64(LL3//2);        S3  = AngularJ64(SS3//2)
                                                L123 = AngularJ64(nrv[7]//2);   S123 = AngularJ64(nrv[8]//2);     J = csfR.J
                            else    error("stop d")     
                            end
                        end
                    end
                    #
                    if  rsh1 != nrsh1   ||   rsh2 != nrsh2   ||   rsh3 != nrsh3      error("stop e")     end
                    
                    # Cycle over all intermediate angular momenta T1, T2, T12, T3
                    T1List = oplus(L1, S1);     T2List = oplus(L2, S2);     T3List = oplus(L3, S3);     T12List = oplus(L12, S12)
                    for  T1 in T1List
                        for  T2 in T2List
                            for  T3 in T3List
                                wa1 = LSjj.getLSjjCoefficient(rsh1.l, N1, LS_jj_qn(w1, QQ1, LL1, SS1, Basics.twice(T1), Nm1, QQm1, JJm1, QQp1, JJp1) )
                                wa2 = LSjj.getLSjjCoefficient(rsh2.l, N2, LS_jj_qn(w2, QQ2, LL2, SS2, Basics.twice(T2), Nm2, QQm2, JJm2, QQp2, JJp2) )
                                wa3 = LSjj.getLSjjCoefficient(rsh3.l, N3, LS_jj_qn(w3, QQ3, LL3, SS3, Basics.twice(T3), Nm3, QQm3, JJm3, QQp3, JJp3) )
                                for  T12 in T12List
                                    ##x println("T1 = $T1   T2 = $T2   T3 = $T3   wa1 = $wa1   wa2 = $wa2   wa3 = $wa3")
                                    if AngularMomentum.triangularDelta(T1,Jm1,Jp1) == 1  &&  T1 == Xp1  &&
                                    AngularMomentum.triangularDelta(T2,Jm2,Jp2) == 1  &&  AngularMomentum.triangularDelta(T3,Jm3,Jp3) == 1
                                        me  = me + sqrt(AngularMomentum.bracket([T1, T2, T3, T12, L12, L123, S12, S123])) * 
                                                AngularMomentum.Wigner_9j(L2,L1,L12, S2,S1,S12, T2,T1,T12) * AngularMomentum.Wigner_9j(J,L123,S123, T3,L3,S3, T12,L12,S12) *
                                                AngularMomentum.phaseFactor([L1, 1, S1, -1, T1,  1, L2, 1, S2, 1, L3, -1, L123, 1, S3, -1, S123, -1, J, -1, T12, -1, T12]) * 
                                                AngularMomentum.phaseFactor([Xp1, 1, Xp1, 1, Jm2, 1, Jp2, 1, T2, 1, T12, 1, T12, 1, Jm3, 1, Jp3, 1, T3]) * 
                                                sqrt(AngularMomentum.bracket([T2, Xm2, T3, Xm3])) * AngularMomentum.Wigner_6j(Jm2, Jp2, T2, T12, Xp1, Xm2) * 
                                                AngularMomentum.Wigner_6j(Jm3, Jp3, T3, J, T12, Xm3) *  wa1 * wa2 * wa3
                                    end
                                end
                            end
                        end
                    end
                end
                push!( mcVector, me)
            end
        end
        #
        return( mcVector )
    end


    """
    `LSjj.extractConfigurationFromCsfNR(csfNR::CsfNR, basisNR::BasisNR)`
        ... to extract the nonrelativistic configuration from the given csfNR, if this is part of basisNR.
            A nonrelativistic conf::Configuration is returned.
    """
    function extractConfigurationFromCsfNR(csfNR::CsfNR, basisNR::BasisNR)
        shells = Dict{Shell,Int64}()
        for  s = 1:length(basisNR.shells)
            shell = basisNR.shells[s];    occ = csfNR.occupation[s]
            if  occ > 0     shells = Base.merge( shells, Dict( shell => occ))   end
        end
        conf   = Configuration( shells, basisNR.NoElectrons)
        return( conf )
    end


    """
    `LSjj.extractConfigurationWeightsOfLevelNR(levelNR::LevelNR)`
        ... to extract the nonrelativistic configurations and their (total) weights that contribute to the given levelNR;
            a dictionary weights::Dict{Configuration,Float64} is returned.
    """
    function extractConfigurationWeightsOfLevelNR(levelNR::LevelNR)
        weights = Dict{String,Float64}()
        for  s = 1:length(levelNR.basis.csfs)
            csf  = levelNR.basis.csfs[s];    weight = abs(levelNR.mc[s])^2
            conf = Base.string(LSjj.extractConfigurationFromCsfNR(csf, levelNR.basis))
            if      haskey(weights, conf)    weights[conf] = weights[conf] + weight
            else    weights = Base.merge( weights, Dict( conf => weight))
            end
        end
        
        return( weights )
    end


    """
    `LSjj.extractOpenShellQNfromCsfNR(csfNR::CsfNR, basisNR::BasisNR)`  
        ... extracts the quantum numbers of the nonrelativistic open shells from the given csfNR. Here, we assume that csfNR is 
            defined in the given basis. An nrQN::Dict{Shell,NTuple{8,Int64}} = Dict(shell => qn ) is returned with 
            qn = (idx, N, w, QQ, LL, SS, LLx, SSx) which contains all the shell quantum numbers and their coupling in csfNR.
    """
    function extractOpenShellQNfromCsfNR(csfNR::CsfNR, basisNR::BasisNR)
        wa = Dict{Shell,NTuple{8,Int64}}()
        for  s = 1:length(basisNR.shells)
            shell = basisNR.shells[s];   occ = csfNR.occupation[s]
            # Now collect the quantum numbers into tuples
            if  0 < occ < 2* (2*shell.l + 1)
                wa = Base.merge( wa, Dict(shell => (s, csfNR.occupation[s], csfNR.w[s], csfNR.QQ[s], 
                                                    Basics.twice(csfNR.shellL[s]), Basics.twice(csfNR.shellS[s]), 
                                                    Basics.twice(csfNR.shellLX[s]), Basics.twice(csfNR.shellSX[s])) ) )
            end
        end
        
        return( wa )
    end


    """
    `LSjj.generateNonrelativisticCsfList(confList::Array{Configuration,1}, shellList::Array{Shell,1})`
        ... to construct from a given list of configuration all possible (nonrelativistic) CSF with regard to the 
            given shellList; a unique list::Array{CsfNR,1} is returned.  
    """
    function generateNonrelativisticCsfList(confList::Array{Configuration,1}, shellList::Array{Shell,1}) 
        
        csfList = CsfNR[];    previousCsfs = CsfNR[]

        # Cycle through all configurations 
        for conf  in  confList
            first = true;    parity  = Basics.determineParity(conf)
            for  shell  in shellList
                if   shell  in  keys(conf.shells)    occ = conf.shells[shell]    else    occ = 0    end
                if   first
                    stateList   = LSjj.provideShellStates(shell, occ)
                    currentCsfs = CsfNR[]
                    for  state in stateList
                        push!( currentCsfs, CsfNR( AngularJ64(0), parity, [state.occ], [state.w], [state.QQ],
                                                [AngularJ64(state.LL//2)], [AngularJ64(state.SS//2)],
                                                [AngularJ64(state.LL//2)], [AngularJ64(state.SS//2)] ) )
                    end
                    previousCsfs = copy(currentCsfs)
                    first        = false
                else
                    # Now support also all couplings of the subshell states with the CSFs that were built-up so far
                    stateList   = LSjj.provideShellStates(shell, occ)
                    currentCsfs = CsfNR[]
                    ##x for  i = 1:length(previousCsfs)    println("generate-aa: ", previousCsfs[i])    end
                    for  csf in  previousCsfs
                        for  state in stateList
                            occupation = deepcopy(csf.occupation);    w = deepcopy(csf.w);    QQ = deepcopy(csf.QQ);    
                            shellL     = deepcopy(csf.shellL);        shellS  = deepcopy(csf.shellS)
                            push!(occupation, state.occ);   push!(w, state.w);   push!(QQ, state.QQ);   
                            push!(shellL, AngularJ64(state.LL//2) );  push!(shellS, AngularJ64(state.SS//2) )  
                            newLXList = oplus( csf.shellLX[end], AngularJ64(state.LL//2) )
                            newSXList = oplus( csf.shellSX[end], AngularJ64(state.SS//2) )
                            for  newLX in newLXList
                                for  newSX in newSXList
                                    shellLX = deepcopy(csf.shellLX);   push!(shellLX, newLX) 
                                    shellSX = deepcopy(csf.shellSX);   push!(shellSX, newSX) 
                                    push!( currentCsfs, CsfNR( AngularJ64(0), parity, occupation, w, QQ, shellL, shellS, shellLX, shellSX) ) 
                                end
                            end
                        end
                    end
                    previousCsfs = copy(currentCsfs)
                end
            end
            append!( csfList, previousCsfs)
        end
        currentCsfs = CsfNR[]
        for  csf in  csfList   #### previousCsfs
            newJList = oplus( csf.shellLX[end], csf.shellSX[end] )
            ##x println("abc: LX = $(csf.shellLX[end])   SX = $(csf.shellSX[end])   newJList = $newJList")
            for  newJ in newJList
                push!( currentCsfs, CsfNR( newJ, csf.parity, csf.occupation, csf.w, csf.QQ, csf.shellL, csf.shellS, csf.shellLX, csf.shellSX) )
            end
        end
        
        ##x println("Final Csfs = $currentCsfs ")
        return( currentCsfs )
    end



    """
    `LSjj.provideShellStates(sh::Shell, occ::Int64)` 
        ... to provide all antisymmetric ShellStates within the nonrelativistic quasi-spin scheme for the given 
            shell and occupation; these shell states are listed by Gaigalas et al., CPC 135 (2002) 219.
            An Array{ShellStateNR,1} is returned.
    """
    function provideShellStates(sh::Shell, occ::Int64)
        wb = LSjj.ShellStateNR[]

        #                                                                                w QQ LL SS
        if       sh.l == 0
            if       occ == 0  || occ == 2      push!( wb, LSjj.ShellStateNR( sh, occ, 0, 1, 0, 0) );       return(wb)
            elseif   occ == 1                   push!( wb, LSjj.ShellStateNR( sh, occ, 0, 0, 0, 1) );       return(wb)
            end
        elseif   sh.l == 1
            if       occ == 0  || occ == 6      push!( wb, LSjj.ShellStateNR( sh, occ, 0, 3, 0, 0) );       return(wb)
            elseif   occ == 1  || occ == 5      push!( wb, LSjj.ShellStateNR( sh, occ, 0, 2, 2, 1) );       return(wb)
            elseif   occ == 2  || occ == 4      push!( wb, LSjj.ShellStateNR( sh, occ, 0, 1, 2, 2) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 0, 3, 0, 0) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 0, 1, 4, 0) );       return(wb)
            elseif   occ == 3                   push!( wb, LSjj.ShellStateNR( sh, occ, 0, 0, 0, 3) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 0, 2, 2, 1) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 0, 0, 4, 1) );       return(wb)
            end
        elseif   sh.l == 2
            if       occ == 0  || occ == 10     push!( wb, LSjj.ShellStateNR( sh, occ, 0, 5, 0, 0) );       return(wb)
            elseif   occ == 1  || occ == 9      push!( wb, LSjj.ShellStateNR( sh, occ, 0, 4, 4, 1) );       return(wb)
            elseif   occ == 2  || occ == 8      push!( wb, LSjj.ShellStateNR( sh, occ, 0, 3, 2, 2) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 0, 3, 6, 2) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 0, 5, 0, 0) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 0, 3, 4, 0) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 0, 3, 8, 0) );       return(wb)
            elseif   occ == 3  || occ == 7      push!( wb, LSjj.ShellStateNR( sh, occ, 0, 2, 2, 3) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 0, 2, 6, 3) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 0, 2, 2, 1) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 0, 4, 4, 1) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 0, 2, 4, 1) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 0, 2, 6, 1) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 0, 2, 8, 1) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 0, 2,10, 1) );       return(wb)
            elseif   occ == 4  || occ == 6      push!( wb, LSjj.ShellStateNR( sh, occ, 0, 1, 4, 4) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 0, 3, 2, 2) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 0, 1, 2, 2) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 0, 1, 4, 2) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 0, 3, 6, 2) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 0, 1, 6, 2) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 0, 1, 8, 2) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 0, 1,10, 2) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 0, 5, 0, 0) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 0, 1, 0, 0) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 0, 3, 4, 0) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 0, 1, 4, 0) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 0, 1, 6, 0) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 0, 3, 8, 0) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 0, 1, 8, 0) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 0, 1,12, 0) );       return(wb)
            elseif   occ == 5                   push!( wb, LSjj.ShellStateNR( sh, occ, 0, 0, 0, 5) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 0, 2, 2, 3) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 0, 0, 4, 3) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 0, 2, 6, 3) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 0, 0, 8, 3) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 0, 0, 0, 1) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 0, 2, 2, 1) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 0, 4, 4, 1) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 0, 2, 4, 1) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 0, 0, 4, 1) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 0, 2, 6, 1) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 0, 0, 6, 1) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 0, 2, 8, 1) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 0, 0, 8, 1) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 0, 2,10, 1) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 0, 0,12, 1) );       return(wb)
            else     error("stop c")
            end
        elseif   sh.l == 3
            if       occ == 0  || occ == 14     push!( wb, LSjj.ShellStateNR( sh, occ, 1, 7, 0, 0) );       return(wb)
            elseif   occ == 1  || occ == 13     push!( wb, LSjj.ShellStateNR( sh, occ, 1, 6, 6, 1) );       return(wb)
            elseif   occ == 2  || occ == 12     push!( wb, LSjj.ShellStateNR( sh, occ, 1, 5, 2, 2) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 1, 5, 6, 2) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 1, 5,10, 2) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 1, 7, 0, 0) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 1, 5, 4, 0) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 1, 5, 8, 0) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 1, 5,12, 0) );       return(wb)
            elseif   occ == 3  || occ == 11     push!( wb, LSjj.ShellStateNR( sh, occ, 1, 4, 0, 3) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 1, 4, 4, 3) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 1, 4, 6, 3) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 1, 4, 8, 3) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 1, 4,12, 3) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 1, 4, 2, 1) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 1, 4, 4, 1) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 2, 4, 4, 1) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 1, 6, 6, 1) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 2, 4, 6, 1) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 1, 4, 8, 1) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 2, 4, 8, 1) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 1, 4,10, 1) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 2, 4,10, 1) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 1, 4,12, 1) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 1, 4,14, 1) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 1, 4,16, 1) );       return(wb)
            elseif   occ == 4  || occ == 10     push!( wb, LSjj.ShellStateNR( sh, occ, 0, 3, 0, 4) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 1, 3, 4, 4) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 1, 3, 6, 4) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 1, 3, 8, 4) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 1, 3,12, 4) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 1, 5, 2, 2) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 2, 3, 2, 2) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 3, 3, 2, 2) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 1, 3, 4, 2) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 2, 3, 4, 2) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 1, 5, 6, 2) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 2, 3, 6, 2) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 3, 3, 6, 2) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 4, 3, 6, 2) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 1, 3, 8, 2) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 2, 3, 8, 2) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 3, 3, 8, 2) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 1, 5,10, 2) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 2, 3,10, 2) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 3, 3,10, 2) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 4, 3,10, 2) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 1, 3,12, 2) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 2, 3,12, 2) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 1, 3,14, 2) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 2, 3,14, 2) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 1, 3,16, 2) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 1, 3,18, 2) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 1, 7, 0, 0) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 2, 3, 0, 0) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 1, 5, 4, 0) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 2, 3, 4, 0) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 3, 3, 4, 0) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 4, 3, 4, 0) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 1, 3, 6, 0) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 1, 5, 8, 0) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 2, 3, 8, 0) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 3, 3, 8, 0) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 4, 3, 8, 0) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 1, 3,10, 0) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 2, 3,10, 0) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 1, 5,12, 0) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 2, 3,12, 0) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 3, 3,12, 0) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 1, 3,14, 0) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 1, 3,16, 0) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 2, 3,16, 0) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 1, 3,20, 0) );       return(wb)
            elseif   occ == 5  || occ ==  9     push!( wb, LSjj.ShellStateNR( sh, occ, 0, 2, 2, 5) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 0, 2, 6, 5) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 0, 2,10, 5) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 1, 4, 0, 3) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 1, 2, 2, 3) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 2, 2, 2, 3) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 1, 4, 4, 3) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 2, 2, 4, 3) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 3, 2, 4, 3) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 1, 4, 6, 3) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 2, 2, 6, 3) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 3, 2, 6, 3) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 4, 2, 6, 3) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 1, 4, 8, 3) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 2, 2, 8, 3) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 3, 2, 8, 3) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 4, 2, 8, 3) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 1, 2,10, 3) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 2, 2,10, 3) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 3, 2,10, 3) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 1, 4,12, 3) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 2, 2,12, 3) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 3, 2,12, 3) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 1, 2,14, 3) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 2, 2,14, 3) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 1, 2,16, 3) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 0, 2,18, 3) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 1, 4, 2, 1) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 2, 2, 2, 1) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 3, 2, 2, 1) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 4, 2, 2, 1) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 1, 4, 4, 1) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 2, 4, 4, 1) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 3, 2, 4, 1) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 4, 2, 4, 1) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 5, 2, 4, 1) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 1, 6, 6, 1) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 2, 4, 6, 1) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 3, 2, 6, 1) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 4, 2, 6, 1) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 5, 2, 6, 1) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 6, 2, 6, 1) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 7, 2, 6, 1) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 1, 4, 8, 1) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 2, 4, 8, 1) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 3, 2, 8, 1) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 4, 2, 8, 1) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 5, 2, 8, 1) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 6, 2, 8, 1) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 1, 4,10, 1) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 2, 4,10, 1) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 3, 2,10, 1) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 4, 2,10, 1) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 5, 2,10, 1) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 6, 2,10, 1) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 7, 2,10, 1) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 1, 4,12, 1) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 2, 2,12, 1) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 3, 2,12, 1) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 4, 2,12, 1) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 5, 2,12, 1) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 1, 4,14, 1) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 2, 2,14, 1) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 3, 2,14, 1) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 4, 2,14, 1) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 5, 2,14, 1) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 1, 4,16, 1) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 2, 2,16, 1) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 3, 2,16, 1) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 1, 2,18, 1) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 2, 2,18, 1) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 1, 2,20, 1) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 0, 2,22, 1) );       return(wb)

            elseif   occ == 6  || occ ==  8     push!( wb, LSjj.ShellStateNR( sh, occ, 0, 1, 6, 6) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 0, 3, 0, 4) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 0, 1, 2, 4) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 1, 3, 4, 4) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 2, 1, 4, 4) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 3, 1, 4, 4) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 1, 3, 6, 4) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 2, 1, 6, 4) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 1, 3, 8, 4) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 2, 1, 8, 4) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 3, 1, 8, 4) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 1, 1,10, 4) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 2, 1,10, 4) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 1, 3,12, 4) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 2, 1,12, 4) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 0, 1,14, 4) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 0, 1,16, 4) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 1, 5, 2, 2) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 2, 3, 2, 2) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 3, 3, 2, 2) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 4, 1, 2, 2) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 5, 1, 2, 2) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 6, 1, 2, 2) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 1, 3, 4, 2) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 2, 3, 4, 2) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 3, 1, 4, 2) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 4, 1, 4, 2) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 5, 1, 4, 2) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 1, 5, 6, 2) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 2, 3, 6, 2) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 3, 3, 6, 2) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 4, 3, 6, 2) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 5, 1, 6, 2) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 6, 1, 6, 2) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 7, 1, 6, 2) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 8, 1, 6, 2) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 9, 1, 6, 2) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 1, 3, 8, 2) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 2, 3, 8, 2) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 3, 3, 8, 2) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 4, 1, 8, 2) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 5, 1, 8, 2) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 6, 1, 8, 2) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 7, 1, 8, 8) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 1, 5,10, 2) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 2, 3,10, 2) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 3, 3,10, 2) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 4, 3,10, 2) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 5, 1,10, 2) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 6, 1,10, 2) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 7, 1,10, 2) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 8, 1,10, 2) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 9, 1,10, 2) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 1, 3,12, 2) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 2, 3,12, 2) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 3, 1,12, 2) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 4, 1,12, 2) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 5, 1,12, 2) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 6, 1,12, 2) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 1, 3,14, 2) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 2, 3,14, 2) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 3, 1,14, 2) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 4, 1,14, 2) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 5, 1,14, 2) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 6, 1,14, 2) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 1, 3,16, 2) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 2, 1,16, 2) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 3, 1,16, 2) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 1, 3,18, 2) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 2, 1,18, 2) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 3, 1,18, 2) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 0, 1,20, 2) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 0, 1,22, 2) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 1, 7, 0, 0) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 2, 3, 0, 0) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 3, 1, 0, 0) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 4, 1, 0, 0) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 0, 1, 2, 0) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 1, 5, 4, 0) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 2, 3, 4, 0) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 3, 3, 4, 0) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 4, 3, 4, 0) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 5, 1, 4, 0) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 6, 1, 4, 0) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 1, 3, 6, 0) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 2, 1, 6, 0) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 3, 1, 6, 0) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 4, 1, 6, 0) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 1, 5, 8, 0) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 2, 3, 8, 0) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 3, 3, 8, 0) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 4, 3, 8, 0) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 5, 1, 8, 0) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 6, 1, 8, 0) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 7, 1, 8, 0) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 8, 1, 8, 0) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 1, 3,10, 0) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 2, 3,10, 0) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 3, 1,10, 0) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 4, 1,10, 0) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 1, 5,12, 0) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 2, 3,12, 0) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 3, 3,12, 0) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 4, 1,12, 0) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 5, 1,12, 0) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 6, 1,12, 0) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 7, 1,12, 0) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 1, 3,14, 0) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 2, 1,14, 0) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 3, 1,14, 0) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 1, 3,16, 0) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 2, 3,16, 0) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 3, 1,16, 0) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 4, 1,16, 0) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 1, 1,18, 0) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 2, 1,18, 0) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 1, 3,20, 0) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 2, 1,20, 0) )
                                                push!( wb, LSjj.ShellStateNR( sh, occ, 0, 1,24, 0) );       return(wb)
            else     error("stop d")
            end

        # For an extension to other subshell-occupations, see the Racah-manual; Table ...
        else  error("stop a")
        end
    end


    """
    `LSjj.shortString(csfNR::CsfNR, basisNR::BasisNR)`
        ... to compose a single-line string of the nonrelativistic csfNR, if this is part of basisNR.
            Only the open shells are printed out explicity; an sa::String is returned.
    """
    function shortString(csfNR::CsfNR, basisNR::BasisNR)
        function  letter(La::AngularJ64) 
            L  = Int64( Basics.twice(La)/2 )
            wa = [ "P", "D", "F", "G", "H", "I", "K", "L", "M", "N", "O", "Q"]
            if      L == 0      wb = "S"    
            elseif  L in 1:12    wb = wa[L]
            else    error("stop a")
            end
            return( wb )
        end
        #
        function  multiplicity(S::AngularJ64) 
            return( string(Basics.twice(S) + 1) )
        end
        
        sa = "  "
        for  s = 1:length(basisNR.shells)
            shell = basisNR.shells[s];    occ = csfNR.occupation[s]
            if  occ == 0  ||  occ == 2(shell.l + 1)     break   end
            sa = sa * "[" * string(shell) * "^" * string(occ) * " ^" * multiplicity(csfNR.shellS[s]) * letter(csfNR.shellL[s]) * 
                    "] " * " ^" * multiplicity(csfNR.shellSX[s]) * letter(csfNR.shellLX[s])
        end
        sa = sa * " ^" * multiplicity(csfNR.shellSX[end]) * letter(csfNR.shellLX[end]) * "_" * string(csfNR.J)
        return( sa )
    end



#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################


    """
    `struct  LS_jj_qn`  ... defines a struct for the generalized quantum numbers (qn) of a LS_jj matrix element.

        + w           ::Int64    ... w quantum number
        + QQ          ::Int64    ... subshell total quasispin 2*Q
        + LL          ::Int64    ... subshell total angular momentum 2*L
        + SS          ::Int64    ... subshell total angular momentum 2*S
        + JJ          ::Int64    ... subshell total angular momentum 2*J
        + Nm          ::Int64    ... occupation Nm
        + Qm          ::Int64    ... subshell quantum number 2*Qm
        + Jm          ::Int64    ... subshell quantum number 2*Jm
        + Qp          ::Int64    ... subshell quantum number 2*Qp
        + Jp          ::Int64    ... subshell quantum number 2*Jp 
    """  
    struct  LS_jj_qn
        w             ::Int64
        QQ            ::Int64
        LL            ::Int64
        SS            ::Int64
        JJ            ::Int64
        Nm            ::Int64
        Qm            ::Int64
        Jm            ::Int64
        Qp            ::Int64
        Jp            ::Int64
    end


    """
    `struct  LS_jj_me`  ... defines a struct for the generalized quantum numbers (qn) of a LS_jj matrix element.

        + qn          ::LS_jj_qn    ... list of quantum numbers of this LS_jj matrix element
        + factor      ::Int64       ... integer factor
        + nom         ::Int64       ... nominator
        + denom       ::Int64       ... denominator
    """  
    struct  LS_jj_me
        qn            ::LS_jj_qn
        factor        ::Int64 
        nom           ::Int64
        denom         ::Int64
    end


    """
    `LSjj.getLSjjCoefficient(l::Int64, N::Int64, qn::LS_jj_qn)`  
        ... to return the value of the LS-jj transformation matrix element for a given set of quantum numbers.
            Note that all (generalized) angular momentum quantum numbers except of l must be given twice the original numbers, 
            i.e. for the quantum numbers Q, L, S, J, Qm, Jm, Qp, Jp.
    """
    function getLSjjCoefficient(l::Int64, N::Int64, qn::LS_jj_qn)
        global  LS_jj_p_3,  LS_jj_p_4,  LS_jj_p_5,  LS_jj_p_6,  LS_jj_d_3,  LS_jj_d_4,  LS_jj_d_5,  LS_jj_d_6,  LS_jj_d_7,
                LS_jj_d_8,  LS_jj_d_9,  LS_jj_d_10,  LS_jj_f_3,  LS_jj_f_4,  LS_jj_f_5,  LS_jj_f_6
        ##x println("getLSjjCoefficient: aa  l = $l   N = $N   qn.w = $(qn.w)")
        #
        wa = 0.
        # Deal separately with the occupations N = 0, 1, 2  before the matrix elements are taken from some predefined lists
        if      N == 0
        elseif  N == 1
            if  qn.Jm == qn.Jp     error("stop a")                 end
            if  qn.JJ == qn.Jm  ||  qn.JJ == qn.Jp   wa = 1.0      end
        elseif  N == 2
            Nm = qn.Nm;    Np = N - Nm;     l2 = l+l;   LL = qn.LL;   SS = qn.SS;   JJ = qn.JJ
            if          Nm == 2     &&  Np == 0     
                if  l2 == 0     twoj1 = twoj2 = 1   else    twoj1 = twoj2 = l2 - 1                end
            elseif      Nm == 1     &&  Np == 1
                if  l2 == 0     twoj1 = twoj2 = 1   else    twoj1 = l2 - 1;   twoj2 = l2 + 1      end
            elseif      Nm == 0     &&  Np == 2
                twoj1 = twoj2 = l2 + 1
            else        error("stop b")                             
            end
            # Compute matrix element
            if          rem(LL+SS,4) != 0
            elseif      rem(LL+SS+JJ,2) != 0
            elseif      rem(twoj1+twoj2+JJ,2) != 0
            else
                if      twoj1 == twoj2   && rem(JJ,4) != 0
                elseif  twoj1 == twoj2      wa = (twoj1+1) * sqrt((LL+1)*(SS+1))
                else    wa = sqrt( 2*(twoj1+1)*(twoj2+1)*(LL+1)*(SS+1) )
                end
                wa = wa * AngularMomentum.Wigner_9j(l,l, LL/2, 0.5, 0.5, SS/2, twoj1/2, twoj2/2, JJ/2)
            end
            # Now take matrix elements from the table
        elseif  l == 0
        elseif  l == 1
            if      N == 3   ## Use data from the array LS_jj_p_3
                for  me  in  LS_jj_p_3    if  qn == me.qn    wa = me.factor * sqrt( me.nom/me.denom );    break    end    end
            elseif  N == 4   ## Use data from the array LS_jj_p_4
                for  me  in  LS_jj_p_4    if  qn == me.qn    wa = me.factor * sqrt( me.nom/me.denom );    break    end    end
            elseif  N == 5   ## Use data from the array LS_jj_p_5
                ##x println("getLSjjCoefficient: bb")
                for  me  in  LS_jj_p_5    
                    ##x println("** qn = $qn       me.qn = $(me.qn)")
                    if  qn == me.qn    wa = me.factor * sqrt( me.nom/me.denom );    break    end    
                end
            elseif  N == 6   ## Use data from the array LS_jj_p_6
                for  me  in  LS_jj_p_6    if  qn == me.qn    wa = me.factor * sqrt( me.nom/me.denom );    break    end    end
            else    error("stop c")
            end
        elseif  l == 2
            if      N == 3   ## Use data from the array LS_jj_d_3
                for  me  in  LS_jj_d_3    if  qn == me.qn    wa = me.factor * sqrt( me.nom/me.denom );    break    end    end
            elseif  N == 4   ## Use data from the array LS_jj_d_4
                for  me  in  LS_jj_d_4    if  qn == me.qn    wa = me.factor * sqrt( me.nom/me.denom );    break    end    end
            elseif  N == 5   ## Use data from the array LS_jj_d_5
                for  me  in  LS_jj_d_5    if  qn == me.qn    wa = me.factor * sqrt( me.nom/me.denom );    break    end    end
            elseif  N == 6   ## Use data from the array LS_jj_d_6
                for  me  in  LS_jj_d_6    if  qn == me.qn    wa = me.factor * sqrt( me.nom/me.denom );    break    end    end
            elseif  N == 7   ## Use data from the array LS_jj_d_7
                for  me  in  LS_jj_d_7    if  qn == me.qn    wa = me.factor * sqrt( me.nom/me.denom );    break    end    end
            elseif  N == 8   ## Use data from the array LS_jj_d_8
                for  me  in  LS_jj_d_8    if  qn == me.qn    wa = me.factor * sqrt( me.nom/me.denom );    break    end    end
            elseif  N == 9   ## Use data from the array LS_jj_d_9
                for  me  in  LS_jj_d_9    if  qn == me.qn    wa = me.factor * sqrt( me.nom/me.denom );    break    end    end
            elseif  N == 10  ## Use data from the array LS_jj_d_10
                for  me  in  LS_jj_d_10   if  qn == me.qn    wa = me.factor * sqrt( me.nom/me.denom );    break    end    end
            else    error("stop d")
            end
        elseif  l == 3
            if      N == 3   ## Use data from the array LS_jj_f_3
                for  me  in  LS_jj_f_3    if  qn == me.qn    wa = me.factor * sqrt( me.nom/me.denom );    break    end    end
            elseif  N == 4   ## Use data from the array LS_jj_f_4
                for  me  in  LS_jj_f_4    if  qn == me.qn    wa = me.factor * sqrt( me.nom/me.denom );    break    end    end
            elseif  N == 5   ## Use data from the array LS_jj_f_5
                for  me  in  LS_jj_f_5    if  qn == me.qn    wa = me.factor * sqrt( me.nom/me.denom );    break    end    end
            elseif  N == 6   ## Use data from the array LS_jj_f_6
                for  me  in  LS_jj_f_6    if  qn == me.qn    wa = me.factor * sqrt( me.nom/me.denom );    break    end    end
            else    error("stop f")
            end
        else        error("stop g")
        end
        ##x println("getLSjjCoefficient: cc     wa = $wa")
        
        return(wa)
    end

    ## p^3  shell;  11 ME in total
    const   LS_jj_p_3 = [ LS_jj_me( LS_jj_qn(0, 2, 2, 1, 1, 1, 0, 1, 2, 0),  1, 1, 1),
                          LS_jj_me( LS_jj_qn(0, 0, 0, 3, 3, 0, 1, 0, 1, 3), -1, 2, 9),
                          LS_jj_me( LS_jj_qn(0, 0, 0, 3, 3, 1, 0, 1, 0, 4),  1, 5, 9),
                          LS_jj_me( LS_jj_qn(0, 0, 0, 3, 3, 2, 1, 0, 1, 3),  1, 2, 9),
                          LS_jj_me( LS_jj_qn(0, 2, 2, 1, 3, 0, 1, 0, 1, 3),  1, 1, 2),
                          LS_jj_me( LS_jj_qn(0, 2, 2, 1, 3, 2, 1, 0, 1, 3),  1, 1, 2),
                          LS_jj_me( LS_jj_qn(0, 0, 4, 1, 3, 0, 1, 0, 1, 3),  1, 5,18),
                          LS_jj_me( LS_jj_qn(0, 0, 4, 1, 3, 1, 0, 1, 0, 4),  1, 4, 9),
                          LS_jj_me( LS_jj_qn(0, 0, 4, 1, 3, 2, 1, 0, 1, 3), -1, 5,18),
                          LS_jj_me( LS_jj_qn(0, 0, 4, 1, 5, 1, 0, 1, 0, 4),  1, 1, 1)  ]
    ## p^4  shell;  9 ME in total
    const   LS_jj_p_4 = [ LS_jj_me( LS_jj_qn(0, 3, 0, 0, 0, 0, 1, 0, 2, 0),  1, 1, 3),
                          LS_jj_me( LS_jj_qn(0, 3, 0, 0, 0, 2, 1, 0, 2, 0),  1, 2, 3),
                          LS_jj_me( LS_jj_qn(0, 1, 2, 2, 0, 0, 1, 0, 2, 0), -1, 2, 3),
                          LS_jj_me( LS_jj_qn(0, 1, 2, 2, 0, 2, 1, 0, 2, 0),  1, 1, 3),
                          LS_jj_me( LS_jj_qn(0, 1, 2, 2, 2, 1, 0, 1, 1, 3),  1, 1, 1),
                          LS_jj_me( LS_jj_qn(0, 1, 2, 2, 4, 1, 0, 1, 1, 3),  1, 1, 3),
                          LS_jj_me( LS_jj_qn(0, 1, 2, 2, 4, 2, 1, 0, 0, 4),  1, 2, 3),
                          LS_jj_me( LS_jj_qn(0, 1, 4, 0, 4, 1, 0, 1, 1, 3), -1, 2, 3),
                          LS_jj_me( LS_jj_qn(0, 1, 4, 0, 4, 2, 1, 0, 0, 4),  1, 1, 3)  ]
    ## p^5  shell;  2 ME in total
    const   LS_jj_p_5 = [ LS_jj_me( LS_jj_qn(0, 2, 2, 1, 1, 1, 0, 1, 2, 0),  1, 1, 1),
                          LS_jj_me( LS_jj_qn(0, 2, 2, 1, 3, 2, 1, 0, 1, 3),  1, 1, 1)  ]
    ## p^6  shell;  1 ME in total
    const   LS_jj_p_6 = [ LS_jj_me( LS_jj_qn(0, 3, 0, 0, 0, 2, 1, 0, 2, 0),  1, 1, 1)  ]

    
    ## d^3  shell;  73 ME in total
    const   LS_jj_d_3 = [ LS_jj_me( LS_jj_qn(0, 2, 2, 3, 1, 1, 1, 3, 1, 4), -1,  7, 15), ## 1
                          LS_jj_me( LS_jj_qn(0, 2, 2, 3, 1, 2, 0, 4, 2, 5), -1,  8, 15), ## 2
                          LS_jj_me( LS_jj_qn(0, 2, 2, 1, 1, 1, 1, 3, 1, 4), -1,  8, 15), ## 3
                          LS_jj_me( LS_jj_qn(0, 2, 2, 1, 1, 2, 0, 4, 2, 5),  1,  7, 15), ## 4
                          LS_jj_me( LS_jj_qn(0, 2, 2, 3, 3, 0, 2, 0, 0, 3), -1, 21,125), ## 5
                          LS_jj_me( LS_jj_qn(0, 2, 2, 3, 3, 1, 1, 3, 3, 0),  1,  8,375), ## 6
                          LS_jj_me( LS_jj_qn(0, 2, 2, 3, 3, 1, 1, 3, 1, 4), -1, 28, 75), ## 7
                          LS_jj_me( LS_jj_qn(0, 2, 2, 3, 3, 2, 0, 4, 2, 5), -1, 28, 75), ## 8
                          LS_jj_me( LS_jj_qn(0, 2, 2, 3, 3, 3, 1, 3, 3, 0), -1,  8,125), ## 9
                          LS_jj_me( LS_jj_qn(0, 2, 2, 1, 3, 0, 2, 0, 0, 3), -1, 12, 25), ## 10
                          LS_jj_me( LS_jj_qn(0, 2, 2, 1, 3, 1, 1, 3, 3, 0), -1,  7,150), ## 11
                          LS_jj_me( LS_jj_qn(0, 2, 2, 1, 3, 1, 1, 3, 1, 4), -1,  1, 15), ## 12
                          LS_jj_me( LS_jj_qn(0, 2, 2, 1, 3, 2, 0, 4, 2, 5),  1,  4, 15), ## 13
                          LS_jj_me( LS_jj_qn(0, 2, 2, 1, 3, 3, 1, 3, 3, 0),  1,  7, 50), ## 14
                          LS_jj_me( LS_jj_qn(0, 4, 4, 1, 3, 1, 1, 3, 3, 0),  1,  3,  4), ## 16
                          LS_jj_me( LS_jj_qn(0, 4, 4, 1, 3, 3, 1, 3, 3, 0),  1,  1,  4), ## 19
                          LS_jj_me( LS_jj_qn(0, 2, 4, 1, 3, 0, 2, 0, 0, 3), -1,  8, 25), ## 20
                          LS_jj_me( LS_jj_qn(0, 2, 4, 1, 3, 1, 1, 3, 3, 0),  1,  7,100), ## 21
                          LS_jj_me( LS_jj_qn(0, 2, 4, 1, 3, 1, 1, 3, 1, 4),  1,  2,  5), ## 22
                          LS_jj_me( LS_jj_qn(0, 2, 4, 1, 3, 3, 1, 3, 3, 0), -1, 21,100), ## 24
                          LS_jj_me( LS_jj_qn(0, 2, 6, 3, 3, 0, 2, 0, 0, 3),  1,  4,125), ## 25
                          LS_jj_me( LS_jj_qn(0, 2, 6, 3, 3, 1, 1, 3, 3, 0),  1, 14,125), ## 26
                          LS_jj_me( LS_jj_qn(0, 2, 6, 3, 3, 1, 1, 3, 1, 4), -1,  4, 25), ## 27
                          LS_jj_me( LS_jj_qn(0, 2, 6, 3, 3, 2, 0, 4, 2, 5),  1,  9, 25), ## 28
                          LS_jj_me( LS_jj_qn(0, 2, 6, 3, 3, 3, 1, 3, 3, 0), -1, 42,125), ## 29
                          LS_jj_me( LS_jj_qn(0, 2, 2, 3, 5, 0, 2, 0, 2, 5),  1, 24,125), ## 30
                          LS_jj_me( LS_jj_qn(0, 2, 2, 3, 5, 1, 1, 3, 1, 4), -1,  1, 25), ## 31
                          LS_jj_me( LS_jj_qn(0, 2, 2, 3, 5, 1, 1, 3, 1, 8), -1, 72,125), ## 32
                          LS_jj_me( LS_jj_qn(0, 2, 2, 3, 5, 2, 2, 0, 2, 5), -1, 24,125), ## 33
                          LS_jj_me( LS_jj_qn(0, 4, 4, 1, 5, 0, 2, 0, 2, 5),  1,  1,  2), ## 35
                          LS_jj_me( LS_jj_qn(0, 4, 4, 1, 5, 2, 2, 0, 2, 5),  1,  1,  2), ## 38
                          LS_jj_me( LS_jj_qn(0, 2, 4, 1, 5, 0, 2, 0, 2, 5), -1,  7,150), ## 40
                          LS_jj_me( LS_jj_qn(0, 2, 4, 1, 5, 1, 1, 3, 1, 4),  1, 16, 35), ## 41
                          LS_jj_me( LS_jj_qn(0, 2, 4, 1, 5, 1, 1, 3, 1, 8), -1, 32,175), ## 42
                          LS_jj_me( LS_jj_qn(0, 2, 4, 1, 5, 2, 2, 0, 2, 5),  1,  7,150), ## 43
                          LS_jj_me( LS_jj_qn(0, 2, 4, 1, 5, 2, 0, 4, 2, 5),  1,  4, 15), ## 44
                          LS_jj_me( LS_jj_qn(0, 2, 6, 3, 5, 0, 2, 0, 2, 5),  1, 28,375), ## 45
                          LS_jj_me( LS_jj_qn(0, 2, 6, 3, 5, 1, 1, 3, 1, 4), -1,  8,175), ## 46
                          LS_jj_me( LS_jj_qn(0, 2, 6, 3, 5, 1, 1, 3, 1, 8),  1,121,875), ## 47
                          LS_jj_me( LS_jj_qn(0, 2, 6, 3, 5, 2, 2, 0, 2, 5), -1, 28,375), ## 48
                          LS_jj_me( LS_jj_qn(0, 2, 6, 3, 5, 2, 0, 4, 2, 5),  1,  2,  3), ## 49
                          LS_jj_me( LS_jj_qn(0, 2, 6, 1, 5, 0, 2, 0, 2, 5),  1, 14, 75), ## 50
                          LS_jj_me( LS_jj_qn(0, 2, 6, 1, 5, 1, 1, 3, 1, 4),  1, 16, 35), ## 51
                          LS_jj_me( LS_jj_qn(0, 2, 6, 1, 5, 1, 1, 3, 1, 8),  1, 18,175), ## 52
                          LS_jj_me( LS_jj_qn(0, 2, 6, 1, 5, 2, 2, 0, 2, 5), -1, 14, 75), ## 53
                          LS_jj_me( LS_jj_qn(0, 2, 6, 1, 5, 2, 0, 4, 2, 5), -1,  1, 15), ## 54
                          LS_jj_me( LS_jj_qn(0, 2, 6, 3, 7, 1, 1, 3, 1, 4), -1,  3, 35), ## 55
                          LS_jj_me( LS_jj_qn(0, 2, 6, 3, 7, 1, 1, 3, 1, 8),  1,  5,  7), ## 56
                          LS_jj_me( LS_jj_qn(0, 2, 6, 3, 7, 2, 0, 4, 2, 5),  1,  1,  5), ## 57
                          LS_jj_me( LS_jj_qn(0, 2, 6, 1, 7, 1, 1, 3, 1, 4),  1,121,280), ## 58
                          LS_jj_me( LS_jj_qn(0, 2, 6, 1, 7, 1, 1, 3, 1, 8),  1, 15, 56), ## 59
                          LS_jj_me( LS_jj_qn(0, 2, 6, 1, 7, 2, 0, 4, 2, 5), -1,  3, 10), ## 60
                          LS_jj_me( LS_jj_qn(0, 2, 8, 1, 7, 1, 1, 3, 1, 4), -1, 27, 56), ## 61
                          LS_jj_me( LS_jj_qn(0, 2, 8, 1, 7, 1, 1, 3, 1, 8),  1,  1, 56), ## 62
                          LS_jj_me( LS_jj_qn(0, 2, 8, 1, 7, 2, 0, 4, 2, 5), -1,  1,  2), ## 63
                          LS_jj_me( LS_jj_qn(0, 2, 6, 3, 9, 0, 2, 0, 0, 9), -1, 12, 25), ## 64
                          LS_jj_me( LS_jj_qn(0, 2, 6, 3, 9, 1, 1, 3, 1, 8),  1, 11, 25), ## 65
                          LS_jj_me( LS_jj_qn(0, 2, 6, 3, 9, 2, 0, 4, 2, 5),  1,  2, 25), ## 66
                          LS_jj_me( LS_jj_qn(0, 2, 8, 1, 9, 0, 2, 0, 0, 9), -1, 54,125), ## 67
                          LS_jj_me( LS_jj_qn(0, 2, 8, 1, 9, 1, 1, 3, 1, 8), -1, 22,125), ## 68
                          LS_jj_me( LS_jj_qn(0, 2, 8, 1, 9, 2, 0, 4, 2, 5), -1, 49,125), ## 69
                          LS_jj_me( LS_jj_qn(0, 2,10, 1, 9, 0, 2, 0, 0, 9), -1, 11,125), ## 70
                          LS_jj_me( LS_jj_qn(0, 2,10, 1, 9, 1, 1, 3, 1, 8), -1, 48,125), ## 71
                          LS_jj_me( LS_jj_qn(0, 2,10, 1, 9, 2, 0, 4, 2, 5),  1, 66,125), ## 72
                          LS_jj_me( LS_jj_qn(0, 2,10, 1,11, 1, 1, 3, 1, 8), -1,  1,  1) ]## 73
    ## d^4  shell;  198 ME in total 					    	       
    const   LS_jj_d_4 = [ LS_jj_me( LS_jj_qn(0, 5, 0, 0, 0, 0, 2, 0, 3, 0),  1,  3,  10), ## 1
                          LS_jj_me( LS_jj_qn(0, 5, 0, 0, 0, 2, 2, 0, 3, 0),  1,  3,   5), ## 3
                          LS_jj_me( LS_jj_qn(0, 5, 0, 0, 0, 4, 2, 0, 3, 0),  1,  1,  10), ## 5
                          LS_jj_me( LS_jj_qn(0, 1, 0, 0, 0, 0, 2, 0, 3, 0), -1,  7, 250), ## 6
                          LS_jj_me( LS_jj_qn(0, 1, 0, 0, 0, 1, 1, 3, 0, 3), -1, 64, 125), ## 7
                          LS_jj_me( LS_jj_qn(0, 1, 0, 0, 0, 2, 2, 0, 3, 0),  1,  7, 125), ## 8
                          LS_jj_me( LS_jj_qn(0, 1, 0, 0, 0, 2, 0, 4, 1, 4),  1,  8,  25), ## 9
                          LS_jj_me( LS_jj_qn(0, 1, 0, 0, 0, 4, 2, 0, 3, 0), -1, 21, 250), ## 10
                          LS_jj_me( LS_jj_qn(0, 3, 2, 2, 0, 0, 2, 0, 3, 0), -1,  8,  15), ## 11
                          LS_jj_me( LS_jj_qn(0, 3, 2, 2, 0, 2, 2, 0, 3, 0),  1,  1,  15), ## 13
                          LS_jj_me( LS_jj_qn(0, 3, 2, 2, 0, 4, 2, 0, 3, 0),  1,  2,   5), ## 15
                          LS_jj_me( LS_jj_qn(0, 1, 2, 2, 0, 0, 2, 0, 3, 0), -1, 28, 375), ## 16
                          LS_jj_me( LS_jj_qn(0, 1, 2, 2, 0, 1, 1, 3, 0, 3),  1, 54, 125), ## 17
                          LS_jj_me( LS_jj_qn(0, 1, 2, 2, 0, 2, 2, 0, 3, 0),  1, 56, 375), ## 18
                          LS_jj_me( LS_jj_qn(0, 1, 2, 2, 0, 2, 0, 4, 1, 4),  1,  3,  25), ## 19
                          LS_jj_me( LS_jj_qn(0, 1, 2, 2, 0, 4, 2, 0, 3, 0), -1, 28, 125), ## 20
                          LS_jj_me( LS_jj_qn(0, 1, 4, 4, 0, 0, 2, 0, 3, 0), -1,  8, 125), ## 21
                          LS_jj_me( LS_jj_qn(0, 1, 4, 4, 0, 1, 1, 3, 0, 3), -1,  7, 125), ## 22
                          LS_jj_me( LS_jj_qn(0, 1, 4, 4, 0, 2, 2, 0, 3, 0),  1, 16, 125), ## 23
                          LS_jj_me( LS_jj_qn(0, 1, 4, 4, 0, 2, 0, 4, 1, 4), -1, 14,  25), ## 24
                          LS_jj_me( LS_jj_qn(0, 1, 4, 4, 0, 4, 2, 0, 3, 0), -1, 24, 125), ## 25
                          LS_jj_me( LS_jj_qn(0, 3, 2, 2, 2, 1, 1, 3, 2, 5),  1,  2,   3), ## 26
                          LS_jj_me( LS_jj_qn(0, 3, 2, 2, 2, 3, 1, 3, 2, 5),  1,  1,   3), ## 29
                          LS_jj_me( LS_jj_qn(0, 1, 2, 2, 2, 1, 1, 3, 2, 5),  1,  7,  75), ## 30
                          LS_jj_me( LS_jj_qn(0, 1, 2, 2, 2, 1, 1, 3, 0, 3),  1,  9,  25), ## 31
                          LS_jj_me( LS_jj_qn(0, 1, 2, 2, 2, 2, 0, 4, 1, 4),  1,  9,  25), ## 32
                          LS_jj_me( LS_jj_qn(0, 1, 2, 2, 2, 3, 1, 3, 2, 5), -1, 14,  75), ## 33
                          LS_jj_me( LS_jj_qn(0, 1, 4, 4, 2, 1, 1, 3, 2, 5),  1,  6,  25), ## 34
                          LS_jj_me( LS_jj_qn(0, 1, 4, 4, 2, 1, 1, 3, 0, 3), -1,  7,  50), ## 35
                          LS_jj_me( LS_jj_qn(0, 1, 4, 4, 2, 2, 0, 4, 1, 4), -1,  7,  50), ## 36
                          LS_jj_me( LS_jj_qn(0, 1, 4, 4, 2, 3, 1, 3, 2, 5), -1, 12,  25), ## 37
                          LS_jj_me( LS_jj_qn(0, 1, 4, 2, 2, 1, 1, 3, 0, 3),  1,  1,   2), ## 39
                          LS_jj_me( LS_jj_qn(0, 1, 4, 2, 2, 2, 0, 4, 1, 4), -1,  1,   2), ## 40
                          LS_jj_me( LS_jj_qn(0, 3, 2, 2, 4, 0, 2, 0, 1, 4),  1, 56, 375), ## 42
                          LS_jj_me( LS_jj_qn(0, 3, 2, 2, 4, 1, 1, 3, 2, 5),  1, 42, 125), ## 43
                          LS_jj_me( LS_jj_qn(0, 3, 2, 2, 4, 2, 2, 0, 1, 4),  1,112, 375), ## 45
                          LS_jj_me( LS_jj_qn(0, 3, 2, 2, 4, 2, 0, 4, 3, 0), -1,  6, 125), ## 46
                          LS_jj_me( LS_jj_qn(0, 3, 2, 2, 4, 3, 1, 3, 2, 5),  1, 21, 125), ## 49
                          LS_jj_me( LS_jj_qn(0, 1, 2, 2, 4, 0, 2, 0, 1, 4), -1, 64, 375), ## 50
                          LS_jj_me( LS_jj_qn(0, 1, 2, 2, 4, 1, 1, 3, 2, 5), -1,  3, 125), ## 51
                          LS_jj_me( LS_jj_qn(0, 1, 2, 2, 4, 1, 1, 3, 0, 3),  1, 27, 125), ## 52
                          LS_jj_me( LS_jj_qn(0, 1, 2, 2, 4, 2, 2, 0, 1, 4),  1, 32, 375), ## 53
                          LS_jj_me( LS_jj_qn(0, 1, 2, 2, 4, 2, 0, 4, 1, 4),  1,  3,  35), ## 55
                          LS_jj_me( LS_jj_qn(0, 1, 2, 2, 4, 2, 0, 4, 1, 8),  1,324, 875), ## 56
                          LS_jj_me( LS_jj_qn(0, 1, 2, 2, 4, 3, 1, 3, 2, 5),  1,  6, 125), ## 57
                          LS_jj_me( LS_jj_qn(0, 1, 4, 4, 4, 0, 2, 0, 1, 4),  1,  2,  25), ## 58
                          LS_jj_me( LS_jj_qn(0, 1, 4, 4, 4, 1, 1, 3, 2, 5),  1,  2,  25), ## 59
                          LS_jj_me( LS_jj_qn(0, 1, 4, 4, 4, 1, 1, 3, 0, 3), -1,  1,  50), ## 60
                          LS_jj_me( LS_jj_qn(0, 1, 4, 4, 4, 2, 2, 0, 1, 4), -1,  1,  25), ## 61
                          LS_jj_me( LS_jj_qn(0, 1, 4, 4, 4, 2, 0, 4, 1, 4), -1,  1,  14), ## 63
                          LS_jj_me( LS_jj_qn(0, 1, 4, 4, 4, 2, 0, 4, 1, 8),  1, 96, 175), ## 64
                          LS_jj_me( LS_jj_qn(0, 1, 4, 4, 4, 3, 1, 3, 2, 5), -1,  4,  25), ## 65
                          LS_jj_me( LS_jj_qn(0, 1, 4, 2, 4, 0, 2, 0, 1, 4), -1,  6,  25), ## 66
                          LS_jj_me( LS_jj_qn(0, 1, 4, 2, 4, 1, 1, 3, 2, 5),  1,  8,  75), ## 67
                          LS_jj_me( LS_jj_qn(0, 1, 4, 2, 4, 1, 1, 3, 0, 3),  1,  3,  50), ## 68
                          LS_jj_me( LS_jj_qn(0, 1, 4, 2, 4, 2, 2, 0, 1, 4),  1,  3,  25), ## 69
                          LS_jj_me( LS_jj_qn(0, 1, 4, 2, 4, 2, 0, 4, 1, 4), -1,  3,  14), ## 71
                          LS_jj_me( LS_jj_qn(0, 1, 4, 2, 4, 2, 0, 4, 1, 8), -1,  8, 175), ## 72
                          LS_jj_me( LS_jj_qn(0, 1, 4, 2, 4, 3, 1, 3, 2, 5), -1, 16,  75), ## 73
                          LS_jj_me( LS_jj_qn(0, 3, 4, 0, 4, 0, 2, 0, 1, 4),  1,  4,  25), ## 74
                          LS_jj_me( LS_jj_qn(0, 3, 4, 0, 4, 1, 1, 3, 2, 5), -1,  4,  25), ## 75
                          LS_jj_me( LS_jj_qn(0, 3, 4, 0, 4, 2, 2, 0, 1, 4),  1,  8,  25), ## 77
                          LS_jj_me( LS_jj_qn(0, 3, 4, 0, 4, 2, 0, 4, 3, 0),  1,  7,  25), ## 78
                          LS_jj_me( LS_jj_qn(0, 3, 4, 0, 4, 3, 1, 3, 2, 5), -1,  2,  25), ## 81
                          LS_jj_me( LS_jj_qn(0, 1, 4, 0, 4, 0, 2, 0, 1, 4), -1,  2,  25), ## 82
                          LS_jj_me( LS_jj_qn(0, 1, 4, 0, 4, 1, 1, 3, 2, 5), -1,  2,  25), ## 83
                          LS_jj_me( LS_jj_qn(0, 1, 4, 0, 4, 1, 1, 3, 0, 3), -1,  8,  25), ## 84
                          LS_jj_me( LS_jj_qn(0, 1, 4, 0, 4, 2, 2, 0, 1, 4),  1,  1,  25), ## 85
                          LS_jj_me( LS_jj_qn(0, 1, 4, 0, 4, 2, 0, 4, 1, 4), -1,  2,   7), ## 87
                          LS_jj_me( LS_jj_qn(0, 1, 4, 0, 4, 2, 0, 4, 1, 8),  1,  6, 175), ## 88
                          LS_jj_me( LS_jj_qn(0, 1, 4, 0, 4, 3, 1, 3, 2, 5),  1,  4,  25), ## 89
                          LS_jj_me( LS_jj_qn(0, 3, 6, 2, 4, 0, 2, 0, 1, 4), -1,  3, 125), ## 90
                          LS_jj_me( LS_jj_qn(0, 3, 6, 2, 4, 1, 1, 3, 2, 5),  1, 64, 375), ## 91
                          LS_jj_me( LS_jj_qn(0, 3, 6, 2, 4, 2, 2, 0, 1, 4), -1,  6, 125), ## 93
                          LS_jj_me( LS_jj_qn(0, 3, 6, 2, 4, 2, 0, 4, 3, 0),  1, 84, 125), ## 94
                          LS_jj_me( LS_jj_qn(0, 3, 6, 2, 4, 3, 1, 3, 2, 5),  1, 32, 375), ## 97
                          LS_jj_me( LS_jj_qn(0, 1, 6, 2, 4, 0, 2, 0, 1, 4), -1, 12, 125), ## 98
                          LS_jj_me( LS_jj_qn(0, 1, 6, 2, 4, 1, 1, 3, 0, 3), -1, 48, 125), ## 100
                          LS_jj_me( LS_jj_qn(0, 1, 6, 2, 4, 2, 2, 0, 1, 4),  1,  6, 125), ## 101
                          LS_jj_me( LS_jj_qn(0, 1, 6, 2, 4, 2, 0, 4, 1, 4),  1, 12,  35), ## 103
                          LS_jj_me( LS_jj_qn(0, 1, 6, 2, 4, 2, 0, 4, 1, 8),  1,  1, 875), ## 104
                          LS_jj_me( LS_jj_qn(0, 1, 6, 2, 4, 3, 1, 3, 2, 5), -1, 32, 375), ## 105
                          LS_jj_me( LS_jj_qn(0, 1, 4, 4, 6, 1, 1, 3, 2, 5),  1,  1,  25), ## 106
                          LS_jj_me( LS_jj_qn(0, 1, 4, 4, 6, 1, 1, 3, 0, 3),  1,  2, 175), ## 107
                          LS_jj_me( LS_jj_qn(0, 1, 4, 4, 6, 1, 1, 3, 0, 9), -1,  3,   7), ## 108
                          LS_jj_me( LS_jj_qn(0, 1, 4, 4, 6, 2, 0, 4, 1, 4), -1,  2, 175), ## 109
                          LS_jj_me( LS_jj_qn(0, 1, 4, 4, 6, 2, 0, 4, 1, 8),  1,  3,   7), ## 110
                          LS_jj_me( LS_jj_qn(0, 1, 4, 4, 6, 3, 1, 3, 2, 5), -1,  2,  25), ## 111
                          LS_jj_me( LS_jj_qn(0, 1, 4, 2, 6, 1, 1, 3, 2, 5),  1,  2,  21), ## 112
                          LS_jj_me( LS_jj_qn(0, 1, 4, 2, 6, 1, 1, 3, 0, 3),  1,  3,  49), ## 113
                          LS_jj_me( LS_jj_qn(0, 1, 4, 2, 6, 1, 1, 3, 0, 9),  1, 18,  49), ## 114
                          LS_jj_me( LS_jj_qn(0, 1, 4, 2, 6, 2, 0, 4, 1, 4), -1, 12,  49), ## 115
                          LS_jj_me( LS_jj_qn(0, 1, 4, 2, 6, 2, 0, 4, 1, 8),  1,  2,  49), ## 116
                          LS_jj_me( LS_jj_qn(0, 1, 4, 2, 6, 3, 1, 3, 2, 5), -1,  4,  21), ## 117
                          LS_jj_me( LS_jj_qn(0, 3, 6, 2, 6, 1, 1, 3, 2, 5),  1,  2,   3), ## 118
                          LS_jj_me( LS_jj_qn(0, 3, 6, 2, 6, 3, 1, 3, 2, 5),  1,  1,   3), ## 123
                          LS_jj_me( LS_jj_qn(0, 1, 6, 2, 6, 1, 1, 3, 2, 5),  1,  1, 600), ## 124
                          LS_jj_me( LS_jj_qn(0, 1, 6, 2, 6, 1, 1, 3, 0, 3), -1, 48, 175), ## 125
                          LS_jj_me( LS_jj_qn(0, 1, 6, 2, 6, 1, 1, 3, 0, 9),  1,  9,  56), ## 126
                          LS_jj_me( LS_jj_qn(0, 1, 6, 2, 6, 2, 0, 4, 1, 4),  1, 48, 175), ## 127
                          LS_jj_me( LS_jj_qn(0, 1, 6, 2, 6, 2, 0, 4, 1, 8),  1,  2,   7), ## 128
                          LS_jj_me( LS_jj_qn(0, 1, 6, 2, 6, 3, 1, 3, 2, 5), -1,  1, 300), ## 129
                          LS_jj_me( LS_jj_qn(0, 1, 6, 0, 6, 1, 1, 3, 2, 5), -1,  1,  10), ## 130
                          LS_jj_me( LS_jj_qn(0, 1, 6, 0, 6, 1, 1, 3, 0, 3),  1, 16,  35), ## 131
                          LS_jj_me( LS_jj_qn(0, 1, 6, 0, 6, 1, 1, 3, 0, 9),  1,  3,  70), ## 132
                          LS_jj_me( LS_jj_qn(0, 1, 6, 0, 6, 2, 0, 4, 1, 4),  1,  1,  35), ## 133
                          LS_jj_me( LS_jj_qn(0, 1, 6, 0, 6, 2, 0, 4, 1, 8),  1,  6,  35), ## 134
                          LS_jj_me( LS_jj_qn(0, 1, 6, 0, 6, 3, 1, 3, 2, 5),  1,  1,   5), ## 135
                          LS_jj_me( LS_jj_qn(0, 1, 8, 2, 6, 1, 1, 3, 2, 5),  1, 27, 280), ## 136
                          LS_jj_me( LS_jj_qn(0, 1, 8, 2, 6, 1, 1, 3, 0, 3),  1, 48, 245), ## 137
                          LS_jj_me( LS_jj_qn(0, 1, 8, 2, 6, 1, 1, 3, 0, 9), -1,  1,1960), ## 138
                          LS_jj_me( LS_jj_qn(0, 1, 8, 2, 6, 2, 0, 4, 1, 4),  1,108, 245), ## 139
                          LS_jj_me( LS_jj_qn(0, 1, 8, 2, 6, 2, 0, 4, 1, 8), -1, 18, 245), ## 140
                          LS_jj_me( LS_jj_qn(0, 1, 8, 2, 6, 3, 1, 3, 2, 5), -1, 27, 140), ## 141
                          LS_jj_me( LS_jj_qn(0, 1, 4, 4, 8, 0, 2, 0, 1, 8),  1, 24, 125), ## 142
                          LS_jj_me( LS_jj_qn(0, 1, 4, 4, 8, 1, 1, 3, 2, 5),  1,  3, 125), ## 143
                          LS_jj_me( LS_jj_qn(0, 1, 4, 4, 8, 1, 1, 3, 0, 9), -1, 11,  25), ## 144
                          LS_jj_me( LS_jj_qn(0, 1, 4, 4, 8, 2, 2, 0, 1, 8), -1, 12, 125), ## 145
                          LS_jj_me( LS_jj_qn(0, 1, 4, 4, 8, 2, 0, 4, 1, 4),  1,  2, 175), ## 146
                          LS_jj_me( LS_jj_qn(0, 1, 4, 4, 8, 2, 0, 4, 1, 8),  1, 33, 175), ## 147
                          LS_jj_me( LS_jj_qn(0, 1, 4, 4, 8, 3, 1, 3, 2, 5), -1,  6, 125), ## 148
                          LS_jj_me( LS_jj_qn(0, 3, 6, 2, 8, 0, 2, 0, 1, 8),  1,  4,  15), ## 149
                          LS_jj_me( LS_jj_qn(0, 3, 6, 2, 8, 1, 1, 3, 2, 5),  1,  2,  15), ## 150
                          LS_jj_me( LS_jj_qn(0, 3, 6, 2, 8, 2, 2, 0, 1, 8),  1,  8,  15), ## 152
                          LS_jj_me( LS_jj_qn(0, 3, 6, 2, 8, 3, 1, 3, 2, 5),  1,  1,  15), ## 155
                          LS_jj_me( LS_jj_qn(0, 1, 6, 2, 8, 0, 2, 0, 1, 8), -1, 64, 375), ## 156
                          LS_jj_me( LS_jj_qn(0, 1, 6, 2, 8, 1, 1, 3, 2, 5),  1,289,3000), ## 157
                          LS_jj_me( LS_jj_qn(0, 1, 6, 2, 8, 1, 1, 3, 0, 9),  1, 11, 200), ## 158
                          LS_jj_me( LS_jj_qn(0, 1, 6, 2, 8, 2, 2, 0, 1, 8),  1, 32, 375), ## 159
                          LS_jj_me( LS_jj_qn(0, 1, 6, 2, 8, 2, 0, 4, 1, 4),  1,  4, 175), ## 160
                          LS_jj_me( LS_jj_qn(0, 1, 6, 2, 8, 2, 0, 4, 1, 8),  1, 66, 175), ## 161
                          LS_jj_me( LS_jj_qn(0, 1, 6, 2, 8, 3, 1, 3, 2, 5), -1,289,1500), ## 162
                          LS_jj_me( LS_jj_qn(0, 1, 8, 2, 8, 0, 2, 0, 1, 8),  1, 16, 125), ## 163
                          LS_jj_me( LS_jj_qn(0, 1, 8, 2, 8, 1, 1, 3, 2, 5),  1, 81,1000), ## 164
                          LS_jj_me( LS_jj_qn(0, 1, 8, 2, 8, 1, 1, 3, 0, 9),  1, 33, 200), ## 165
                          LS_jj_me( LS_jj_qn(0, 1, 8, 2, 8, 2, 2, 0, 1, 8), -1,  8, 125), ## 166
                          LS_jj_me( LS_jj_qn(0, 1, 8, 2, 8, 2, 0, 4, 1, 4),  1, 48, 175), ## 167
                          LS_jj_me( LS_jj_qn(0, 1, 8, 2, 8, 2, 0, 4, 1, 8), -1, 22, 175), ## 168
                          LS_jj_me( LS_jj_qn(0, 1, 8, 2, 8, 3, 1, 3, 2, 5), -1, 81, 500), ## 169
                          LS_jj_me( LS_jj_qn(0, 3, 8, 0, 8, 0, 2, 0, 1, 8),  1,  1,  15), ## 170
                          LS_jj_me( LS_jj_qn(0, 3, 8, 0, 8, 1, 1, 3, 2, 5), -1,  8,  15), ## 171
                          LS_jj_me( LS_jj_qn(0, 3, 8, 0, 8, 2, 2, 0, 1, 8),  1,  2,  15), ## 173
                          LS_jj_me( LS_jj_qn(0, 3, 8, 0, 8, 3, 1, 3, 2, 5), -1,  4,  15), ## 176
                          LS_jj_me( LS_jj_qn(0, 1, 8, 0, 8, 0, 2, 0, 1, 8), -1, 44, 375), ## 177
                          LS_jj_me( LS_jj_qn(0, 1, 8, 0, 8, 1, 1, 3, 2, 5), -1, 11, 750), ## 178
                          LS_jj_me( LS_jj_qn(0, 1, 8, 0, 8, 1, 1, 3, 0, 9), -1,  9,  50), ## 179
                          LS_jj_me( LS_jj_qn(0, 1, 8, 0, 8, 2, 2, 0, 1, 8),  1, 22, 375), ## 180
                          LS_jj_me( LS_jj_qn(0, 1, 8, 0, 8, 2, 0, 4, 1, 4),  1, 99, 175), ## 181
                          LS_jj_me( LS_jj_qn(0, 1, 8, 0, 8, 2, 0, 4, 1, 8), -1,  6, 175), ## 182
                          LS_jj_me( LS_jj_qn(0, 1, 8, 0, 8, 3, 1, 3, 2, 5),  1, 11, 375), ## 183
                          LS_jj_me( LS_jj_qn(0, 1,10, 2, 8, 0, 2, 0, 1, 8),  1, 22, 375), ## 184
                          LS_jj_me( LS_jj_qn(0, 1,10, 2, 8, 1, 1, 3, 2, 5), -1, 44, 375), ## 185
                          LS_jj_me( LS_jj_qn(0, 1,10, 2, 8, 1, 1, 3, 0, 9),  1,  4,  25), ## 186
                          LS_jj_me( LS_jj_qn(0, 1,10, 2, 8, 2, 2, 0, 1, 8), -1, 11, 375), ## 187
                          LS_jj_me( LS_jj_qn(0, 1,10, 2, 8, 2, 0, 4, 1, 4),  1, 22, 175), ## 188
                          LS_jj_me( LS_jj_qn(0, 1,10, 2, 8, 2, 0, 4, 1, 8),  1, 48, 175), ## 189
                          LS_jj_me( LS_jj_qn(0, 1,10, 2, 8, 3, 1, 3, 2, 5),  1, 88, 375), ## 190
                          LS_jj_me( LS_jj_qn(0, 1, 8, 2,10, 1, 1, 3, 0, 9),  1, 16,  25), ## 191
                          LS_jj_me( LS_jj_qn(0, 1, 8, 2,10, 2, 0, 4, 1, 8), -1,  9,  25), ## 192
                          LS_jj_me( LS_jj_qn(0, 1,10, 2,10, 1, 1, 3, 0, 9),  1,  9,  25), ## 193
                          LS_jj_me( LS_jj_qn(0, 1,10, 2,10, 2, 0, 4, 1, 8),  1, 16,  25), ## 194
                          LS_jj_me( LS_jj_qn(0, 1,10, 2,12, 1, 1, 3, 0, 9),  1,  3,   5), ## 195
                          LS_jj_me( LS_jj_qn(0, 1,10, 2,12, 2, 0, 4, 1, 8),  1,  2,   5), ## 196
                          LS_jj_me( LS_jj_qn(0, 1,12, 0,12, 1, 1, 3, 0, 9),  1,  2,   5), ## 197
                          LS_jj_me( LS_jj_qn(0, 1,12, 0,12, 2, 0, 4, 1, 8), -1,  3,   5) ]## 198
    ## d^5  shell;  249 ME in total					    		
    const   LS_jj_d_5 = [ LS_jj_me( LS_jj_qn(0, 0, 0, 1, 1, 1, 1, 3, 1, 4),  1,   2,	5), ## 1
                          LS_jj_me( LS_jj_qn(0, 0, 0, 1, 1, 2, 0, 4, 0, 3),  1,   1,	5), ## 3
                          LS_jj_me( LS_jj_qn(0, 0, 0, 1, 1, 3, 1, 3, 1, 4), -1,   2,	5), ## 4
                          LS_jj_me( LS_jj_qn(0, 2, 2, 3, 1, 1, 1, 3, 1, 4), -1,   7,   30), ## 5
                          LS_jj_me( LS_jj_qn(0, 2, 2, 3, 1, 2, 0, 4, 2, 5), -1,   8,   15), ## 6
                          LS_jj_me( LS_jj_qn(0, 2, 2, 3, 1, 3, 1, 3, 1, 4), -1,   7,   30), ## 8
                          LS_jj_me( LS_jj_qn(0, 2, 2, 1, 1, 1, 1, 3, 1, 4), -1,   4,   15), ## 9
                          LS_jj_me( LS_jj_qn(0, 2, 2, 1, 1, 2, 0, 4, 2, 5),  1,   7,   15), ## 10
                          LS_jj_me( LS_jj_qn(0, 2, 2, 1, 1, 3, 1, 3, 1, 4), -1,   4,   15), ## 12
                          LS_jj_me( LS_jj_qn(0, 0, 4, 3, 1, 1, 1, 3, 1, 4),  1,   1,   10), ## 13
                          LS_jj_me( LS_jj_qn(0, 0, 4, 3, 1, 2, 0, 4, 0, 3), -1,   4,	5), ## 15
                          LS_jj_me( LS_jj_qn(0, 0, 4, 3, 1, 3, 1, 3, 1, 4), -1,   1,   10), ## 16
                          LS_jj_me( LS_jj_qn(0, 2, 2, 3, 3, 1, 1, 3, 3, 0),  1,  16,  375), ## 17
                          LS_jj_me( LS_jj_qn(0, 2, 2, 3, 3, 1, 1, 3, 1, 4), -1,  14,   75), ## 18
                          LS_jj_me( LS_jj_qn(0, 2, 2, 3, 3, 2, 2, 0, 0, 3), -1,  21,  125), ## 19
                          LS_jj_me( LS_jj_qn(0, 2, 2, 3, 3, 2, 0, 4, 2, 5), -1,  28,   75), ## 20
                          LS_jj_me( LS_jj_qn(0, 2, 2, 3, 3, 3, 1, 3, 3, 0), -1,  16,  375), ## 22
                          LS_jj_me( LS_jj_qn(0, 2, 2, 3, 3, 3, 1, 3, 1, 4), -1,  14,   75), ## 23
                          LS_jj_me( LS_jj_qn(0, 2, 2, 1, 3, 1, 1, 3, 3, 0), -1,   7,   75), ## 24
                          LS_jj_me( LS_jj_qn(0, 2, 2, 1, 3, 1, 1, 3, 1, 4), -1,   1,   30), ## 25
                          LS_jj_me( LS_jj_qn(0, 2, 2, 1, 3, 2, 2, 0, 0, 3), -1,  12,   25), ## 26
                          LS_jj_me( LS_jj_qn(0, 2, 2, 1, 3, 2, 0, 4, 2, 5),  1,   4,   15), ## 27
                          LS_jj_me( LS_jj_qn(0, 2, 2, 1, 3, 3, 1, 3, 3, 0),  1,   7,   75), ## 29
                          LS_jj_me( LS_jj_qn(0, 2, 2, 1, 3, 3, 1, 3, 1, 4), -1,   1,   30), ## 30
                          LS_jj_me( LS_jj_qn(0, 4, 4, 1, 3, 1, 1, 3, 3, 0),  1,   1,	2), ## 31
                          LS_jj_me( LS_jj_qn(0, 4, 4, 1, 3, 3, 1, 3, 3, 0),  1,   1,	2), ## 36
                          LS_jj_me( LS_jj_qn(0, 2, 4, 1, 3, 1, 1, 3, 3, 0),  1,   7,   50), ## 38
                          LS_jj_me( LS_jj_qn(0, 2, 4, 1, 3, 1, 1, 3, 1, 4),  1,   1,    5), ## 39
                          LS_jj_me( LS_jj_qn(0, 2, 4, 1, 3, 2, 2, 0, 0, 3), -1,   8,   25), ## 40
                          LS_jj_me( LS_jj_qn(0, 2, 4, 1, 3, 3, 1, 3, 3, 0), -1,   7,   50), ## 43
                          LS_jj_me( LS_jj_qn(0, 2, 4, 1, 3, 3, 1, 3, 1, 4),  1,   1,    5), ## 44
                          LS_jj_me( LS_jj_qn(0, 0, 4, 3, 3, 1, 1, 3, 1, 4),  1,   2,    5), ## 46
                          LS_jj_me( LS_jj_qn(0, 0, 4, 3, 3, 2, 0, 4, 0, 3), -1,   1,    5), ## 49
                          LS_jj_me( LS_jj_qn(0, 0, 4, 3, 3, 3, 1, 3, 1, 4), -1,   2,    5), ## 51
                          LS_jj_me( LS_jj_qn(0, 0, 4, 1, 3, 1, 1, 3, 1, 4), -1,   1,   10), ## 53
                          LS_jj_me( LS_jj_qn(0, 0, 4, 1, 3, 2, 0, 4, 0, 3), -1,   4,    5), ## 56
                          LS_jj_me( LS_jj_qn(0, 0, 4, 1, 3, 3, 1, 3, 1, 4),  1,   1,   10), ## 58
                          LS_jj_me( LS_jj_qn(0, 2, 6, 3, 3, 1, 1, 3, 3, 0),  1,  28,  125), ## 59
                          LS_jj_me( LS_jj_qn(0, 2, 6, 3, 3, 1, 1, 3, 1, 4), -1,   2,   25), ## 60
                          LS_jj_me( LS_jj_qn(0, 2, 6, 3, 3, 2, 2, 0, 0, 3),  1,   4,  125), ## 61
                          LS_jj_me( LS_jj_qn(0, 2, 6, 3, 3, 2, 0, 4, 2, 5),  1,   9,   25), ## 62
                          LS_jj_me( LS_jj_qn(0, 2, 6, 3, 3, 3, 1, 3, 3, 0), -1,  28,  125), ## 64
                          LS_jj_me( LS_jj_qn(0, 2, 6, 3, 3, 3, 1, 3, 1, 4), -1,   2,   25), ## 65
                          LS_jj_me( LS_jj_qn(0, 0, 0, 5, 5, 0, 2, 0, 2, 5), -1,  24,  625), ## 66
                          LS_jj_me( LS_jj_qn(0, 0, 0, 5, 5, 1, 1, 3, 1, 4),  1,   2,  125), ## 67
                          LS_jj_me( LS_jj_qn(0, 0, 0, 5, 5, 1, 1, 3, 1, 8),  1, 144,  625), ## 68
                          LS_jj_me( LS_jj_qn(0, 0, 0, 5, 5, 2, 2, 0, 2, 5),  1,  24,  625), ## 69
                          LS_jj_me( LS_jj_qn(0, 0, 0, 5, 5, 2, 0, 4, 0, 3), -1,   1,  125), ## 71
                          LS_jj_me( LS_jj_qn(0, 0, 0, 5, 5, 2, 0, 4, 0, 9), -1,  48,  125), ## 72
                          LS_jj_me( LS_jj_qn(0, 0, 0, 5, 5, 3, 1, 3, 1, 4), -1,   2,  125), ## 73
                          LS_jj_me( LS_jj_qn(0, 0, 0, 5, 5, 3, 1, 3, 1, 8), -1, 144,  625), ## 74
                          LS_jj_me( LS_jj_qn(0, 0, 0, 5, 5, 4, 2, 0, 2, 5), -1,  24,  625), ## 75
                          LS_jj_me( LS_jj_qn(0, 2, 2, 3, 5, 0, 2, 0, 2, 5),  1,  24,  125), ## 76
                          LS_jj_me( LS_jj_qn(0, 2, 2, 3, 5, 1, 1, 3, 1, 4), -1,   1,   50), ## 77
                          LS_jj_me( LS_jj_qn(0, 2, 2, 3, 5, 1, 1, 3, 1, 8), -1,  36,  125), ## 78
                          LS_jj_me( LS_jj_qn(0, 2, 2, 3, 5, 3, 1, 3, 1, 4), -1,   1,   50), ## 83
                          LS_jj_me( LS_jj_qn(0, 2, 2, 3, 5, 3, 1, 3, 1, 8), -1,  36,  125), ## 84
                          LS_jj_me( LS_jj_qn(0, 2, 2, 3, 5, 4, 2, 0, 2, 5), -1,  24,  125), ## 85
                          LS_jj_me( LS_jj_qn(0, 4, 4, 1, 5, 0, 2, 0, 2, 5),  1,   1,    6), ## 86
                          LS_jj_me( LS_jj_qn(0, 4, 4, 1, 5, 2, 2, 0, 2, 5),  1,   2,    3), ## 89
                          LS_jj_me( LS_jj_qn(0, 4, 4, 1, 5, 4, 2, 0, 2, 5),  1,   1,    6), ## 95
                          LS_jj_me( LS_jj_qn(0, 2, 4, 1, 5, 0, 2, 0, 2, 5), -1,   7,  150), ## 96
                          LS_jj_me( LS_jj_qn(0, 2, 4, 1, 5, 1, 1, 3, 1, 4),  1,   8,   35), ## 97
                          LS_jj_me( LS_jj_qn(0, 2, 4, 1, 5, 1, 1, 3, 1, 8), -1,  16,  175), ## 98
                          LS_jj_me( LS_jj_qn(0, 2, 4, 1, 5, 2, 0, 4, 2, 5),  1,   4,   15), ## 100
                          LS_jj_me( LS_jj_qn(0, 2, 4, 1, 5, 3, 1, 3, 1, 4),  1,   8,   35), ## 103
                          LS_jj_me( LS_jj_qn(0, 2, 4, 1, 5, 3, 1, 3, 1, 8), -1,  16,  175), ## 104
                          LS_jj_me( LS_jj_qn(0, 2, 4, 1, 5, 4, 2, 0, 2, 5),  1,   7,  150), ## 105
                          LS_jj_me( LS_jj_qn(0, 0, 4, 3, 5, 0, 2, 0, 2, 5), -1,   8,   75), ## 106
                          LS_jj_me( LS_jj_qn(0, 0, 4, 3, 5, 1, 1, 3, 1, 4),  1,  81,  490), ## 107
                          LS_jj_me( LS_jj_qn(0, 0, 4, 3, 5, 1, 1, 3, 1, 8), -1,   4, 1225), ## 108
                          LS_jj_me( LS_jj_qn(0, 0, 4, 3, 5, 2, 2, 0, 2, 5),  1,   8,   75), ## 109
                          LS_jj_me( LS_jj_qn(0, 0, 4, 3, 5, 2, 0, 4, 0, 3), -1,  36,  245), ## 111
                          LS_jj_me( LS_jj_qn(0, 0, 4, 3, 5, 2, 0, 4, 0, 9),  1,  48,  245), ## 112
                          LS_jj_me( LS_jj_qn(0, 0, 4, 3, 5, 3, 1, 3, 1, 4), -1,  81,  490), ## 113
                          LS_jj_me( LS_jj_qn(0, 0, 4, 3, 5, 3, 1, 3, 1, 8),  1,   4, 1225), ## 114
                          LS_jj_me( LS_jj_qn(0, 0, 4, 3, 5, 4, 2, 0, 2, 5), -1,   8,   75), ## 115
                          LS_jj_me( LS_jj_qn(0, 0, 4, 1, 5, 0, 2, 0, 2, 5),  1,   7,   75), ## 116
                          LS_jj_me( LS_jj_qn(0, 0, 4, 1, 5, 1, 1, 3, 1, 4),  1,   4,   35), ## 117
                          LS_jj_me( LS_jj_qn(0, 0, 4, 1, 5, 1, 1, 3, 1, 8), -1,   8,  175), ## 118
                          LS_jj_me( LS_jj_qn(0, 0, 4, 1, 5, 2, 2, 0, 2, 5), -1,   7,   75), ## 119
                          LS_jj_me( LS_jj_qn(0, 0, 4, 1, 5, 2, 0, 4, 0, 3), -1,   8,   35), ## 121
                          LS_jj_me( LS_jj_qn(0, 0, 4, 1, 5, 2, 0, 4, 0, 9), -1,   6,   35), ## 122
                          LS_jj_me( LS_jj_qn(0, 0, 4, 1, 5, 3, 1, 3, 1, 4), -1,   4,   35), ## 123
                          LS_jj_me( LS_jj_qn(0, 0, 4, 1, 5, 3, 1, 3, 1, 8),  1,   8,  175), ## 124
                          LS_jj_me( LS_jj_qn(0, 0, 4, 1, 5, 4, 2, 0, 2, 5),  1,   7,   75), ## 125
                          LS_jj_me( LS_jj_qn(0, 2, 6, 3, 5, 0, 2, 0, 2, 5),  1,  28,  375), ## 126
                          LS_jj_me( LS_jj_qn(0, 2, 6, 3, 5, 1, 1, 3, 1, 4), -1,   4,  175), ## 127
                          LS_jj_me( LS_jj_qn(0, 2, 6, 3, 5, 1, 1, 3, 1, 8),  1, 121, 1750), ## 128
                          LS_jj_me( LS_jj_qn(0, 2, 6, 3, 5, 2, 0, 4, 2, 5),  1,   2,    3), ## 130
                          LS_jj_me( LS_jj_qn(0, 2, 6, 3, 5, 3, 1, 3, 1, 4), -1,   4,  175), ## 133
                          LS_jj_me( LS_jj_qn(0, 2, 6, 3, 5, 3, 1, 3, 1, 8),  1, 121, 1750), ## 134
                          LS_jj_me( LS_jj_qn(0, 2, 6, 3, 5, 4, 2, 0, 2, 5), -1,  28,  375), ## 135
                          LS_jj_me( LS_jj_qn(0, 2, 6, 1, 5, 0, 2, 0, 2, 5),  1,  14,   75), ## 136
                          LS_jj_me( LS_jj_qn(0, 2, 6, 1, 5, 1, 1, 3, 1, 4),  1,   8,   35), ## 137
                          LS_jj_me( LS_jj_qn(0, 2, 6, 1, 5, 1, 1, 3, 1, 8),  1,   9,  175), ## 138
                          LS_jj_me( LS_jj_qn(0, 2, 6, 1, 5, 2, 0, 4, 2, 5), -1,   1,   15), ## 140
                          LS_jj_me( LS_jj_qn(0, 2, 6, 1, 5, 3, 1, 3, 1, 4),  1,   8,   35), ## 143
                          LS_jj_me( LS_jj_qn(0, 2, 6, 1, 5, 3, 1, 3, 1, 8),  1,   9,  175), ## 144
                          LS_jj_me( LS_jj_qn(0, 2, 6, 1, 5, 4, 2, 0, 2, 5), -1,  14,   75), ## 145
                          LS_jj_me( LS_jj_qn(0, 0, 6, 1, 5, 0, 2, 0, 2, 5),  1,  14,  375), ## 146
                          LS_jj_me( LS_jj_qn(0, 0, 6, 1, 5, 1, 1, 3, 1, 4), -1,   8,  175), ## 147
                          LS_jj_me( LS_jj_qn(0, 0, 6, 1, 5, 1, 1, 3, 1, 8),  1, 121,  875), ## 148
                          LS_jj_me( LS_jj_qn(0, 0, 6, 1, 5, 2, 2, 0, 2, 5), -1,  14,  375), ## 149
                          LS_jj_me( LS_jj_qn(0, 0, 6, 1, 5, 2, 0, 4, 0, 3), -1,  64,  175), ## 151
                          LS_jj_me( LS_jj_qn(0, 0, 6, 1, 5, 2, 0, 4, 0, 9),  1,  27,  175), ## 152
                          LS_jj_me( LS_jj_qn(0, 0, 6, 1, 5, 3, 1, 3, 1, 4),  1,   8,  175), ## 153
                          LS_jj_me( LS_jj_qn(0, 0, 6, 1, 5, 3, 1, 3, 1, 8), -1, 121,  875), ## 154
                          LS_jj_me( LS_jj_qn(0, 0, 6, 1, 5, 4, 2, 0, 2, 5),  1,  14,  375), ## 155
                          LS_jj_me( LS_jj_qn(0, 0, 8, 3, 5, 0, 2, 0, 2, 5), -1,  36,  625), ## 156
                          LS_jj_me( LS_jj_qn(0, 0, 8, 3, 5, 1, 1, 3, 1, 4), -1, 972, 6125), ## 157
                          LS_jj_me( LS_jj_qn(0, 0, 8, 3, 5, 1, 1, 3, 1, 8), -1,5043,61250), ## 158
                          LS_jj_me( LS_jj_qn(0, 0, 8, 3, 5, 2, 2, 0, 2, 5),  1,  36,  625), ## 159
                          LS_jj_me( LS_jj_qn(0, 0, 8, 3, 5, 2, 0, 4, 0, 3), -1,1536, 6125), ## 161
                          LS_jj_me( LS_jj_qn(0, 0, 8, 3, 5, 2, 0, 4, 0, 9), -1, 578, 6125), ## 162
                          LS_jj_me( LS_jj_qn(0, 0, 8, 3, 5, 3, 1, 3, 1, 4),  1, 972, 6125), ## 163
                          LS_jj_me( LS_jj_qn(0, 0, 8, 3, 5, 3, 1, 3, 1, 8),  1,5043,61250), ## 164
                          LS_jj_me( LS_jj_qn(0, 0, 8, 3, 5, 4, 2, 0, 2, 5), -1,  36,  625), ## 165
                          LS_jj_me( LS_jj_qn(0, 0, 4, 3, 7, 1, 1, 3, 1, 4), -1,   2,  245), ## 166
                          LS_jj_me( LS_jj_qn(0, 0, 4, 3, 7, 1, 1, 3, 1, 8), -1,  54,  245), ## 167
                          LS_jj_me( LS_jj_qn(0, 0, 4, 3, 7, 2, 0, 4, 0, 3),  1,   1,  245), ## 169
                          LS_jj_me( LS_jj_qn(0, 0, 4, 3, 7, 2, 0, 4, 0, 9),  1, 132,  245), ## 170
                          LS_jj_me( LS_jj_qn(0, 0, 4, 3, 7, 3, 1, 3, 1, 4),  1,   2,  245), ## 171
                          LS_jj_me( LS_jj_qn(0, 0, 4, 3, 7, 3, 1, 3, 1, 8),  1,  54,  245), ## 172
                          LS_jj_me( LS_jj_qn(0, 2, 6, 3, 7, 1, 1, 3, 1, 4), -1,   3,   70), ## 173
                          LS_jj_me( LS_jj_qn(0, 2, 6, 3, 7, 1, 1, 3, 1, 8),  1,   5,   14), ## 174
                          LS_jj_me( LS_jj_qn(0, 2, 6, 3, 7, 2, 0, 4, 2, 5),  1,   1,    5), ## 175
                          LS_jj_me( LS_jj_qn(0, 2, 6, 3, 7, 3, 1, 3, 1, 4), -1,   3,   70), ## 178
                          LS_jj_me( LS_jj_qn(0, 2, 6, 3, 7, 3, 1, 3, 1, 8),  1,   5,   14), ## 179
                          LS_jj_me( LS_jj_qn(0, 2, 6, 1, 7, 1, 1, 3, 1, 4),  1, 121,  560), ## 180
                          LS_jj_me( LS_jj_qn(0, 2, 6, 1, 7, 1, 1, 3, 1, 8),  1,  15,  112), ## 181
                          LS_jj_me( LS_jj_qn(0, 2, 6, 1, 7, 2, 0, 4, 2, 5), -1,   3,   10), ## 182
                          LS_jj_me( LS_jj_qn(0, 2, 6, 1, 7, 3, 1, 3, 1, 4),  1, 121,  560), ## 185
                          LS_jj_me( LS_jj_qn(0, 2, 6, 1, 7, 3, 1, 3, 1, 8),  1,  15,  112), ## 186
                          LS_jj_me( LS_jj_qn(0, 0, 6, 1, 7, 1, 1, 3, 1, 4), -1,  25,  112), ## 187
                          LS_jj_me( LS_jj_qn(0, 0, 6, 1, 7, 1, 1, 3, 1, 8),  1, 243, 2800), ## 188
                          LS_jj_me( LS_jj_qn(0, 0, 6, 1, 7, 2, 0, 4, 0, 3), -1,   2,    7), ## 190
                          LS_jj_me( LS_jj_qn(0, 0, 6, 1, 7, 2, 0, 4, 0, 9),  1,  33,  350), ## 191
                          LS_jj_me( LS_jj_qn(0, 0, 6, 1, 7, 3, 1, 3, 1, 4),  1,  25,  112), ## 192
                          LS_jj_me( LS_jj_qn(0, 0, 6, 1, 7, 3, 1, 3, 1, 8), -1, 243, 2800), ## 193
                          LS_jj_me( LS_jj_qn(0, 2, 8, 1, 7, 1, 1, 3, 1, 4), -1,  27,  112), ## 194
                          LS_jj_me( LS_jj_qn(0, 2, 8, 1, 7, 1, 1, 3, 1, 8),  1,   1,  112), ## 195
                          LS_jj_me( LS_jj_qn(0, 2, 8, 1, 7, 2, 0, 4, 2, 5), -1,   1,    2), ## 196
                          LS_jj_me( LS_jj_qn(0, 2, 8, 1, 7, 3, 1, 3, 1, 4), -1,  27,  112), ## 199
                          LS_jj_me( LS_jj_qn(0, 2, 8, 1, 7, 3, 1, 3, 1, 8),  1,   1,  112), ## 200
                          LS_jj_me( LS_jj_qn(0, 0, 8, 3, 7, 1, 1, 3, 1, 4), -1,   9,   98), ## 201
                          LS_jj_me( LS_jj_qn(0, 0, 8, 3, 7, 1, 1, 3, 1, 8), -1,1369, 7350), ## 202
                          LS_jj_me( LS_jj_qn(0, 0, 8, 3, 7, 2, 0, 4, 0, 3), -1,   4,   49), ## 204
                          LS_jj_me( LS_jj_qn(0, 0, 8, 3, 7, 2, 0, 4, 0, 9), -1,1331, 3675), ## 205
                          LS_jj_me( LS_jj_qn(0, 0, 8, 3, 7, 3, 1, 3, 1, 4),  1,   9,   98), ## 206
                          LS_jj_me( LS_jj_qn(0, 0, 8, 3, 7, 3, 1, 3, 1, 8),  1,1369, 7350), ## 207
                          LS_jj_me( LS_jj_qn(0, 0, 8, 1, 7, 1, 1, 3, 1, 4), -1,  99,  560), ## 208
                          LS_jj_me( LS_jj_qn(0, 0, 8, 1, 7, 1, 1, 3, 1, 8),  1,  11, 1680), ## 209
                          LS_jj_me( LS_jj_qn(0, 0, 8, 1, 7, 2, 0, 4, 0, 3),  1,  22,   35), ## 211
                          LS_jj_me( LS_jj_qn(0, 0, 8, 1, 7, 2, 0, 4, 0, 9), -1,   1,  210), ## 212
                          LS_jj_me( LS_jj_qn(0, 0, 8, 1, 7, 3, 1, 3, 1, 4),  1,  99,  560), ## 213
                          LS_jj_me( LS_jj_qn(0, 0, 8, 1, 7, 3, 1, 3, 1, 8), -1,  11, 1680), ## 214
                          LS_jj_me( LS_jj_qn(0, 2, 6, 3, 9, 1, 1, 3, 1, 8),  1,  11,   50), ## 215
                          LS_jj_me( LS_jj_qn(0, 2, 6, 3, 9, 2, 2, 0, 0, 9), -1,  12,   25), ## 216
                          LS_jj_me( LS_jj_qn(0, 2, 6, 3, 9, 2, 0, 4, 2, 5),  1,   2,   25), ## 217
                          LS_jj_me( LS_jj_qn(0, 2, 6, 3, 9, 3, 1, 3, 1, 8),  1,  11,   50), ## 219
                          LS_jj_me( LS_jj_qn(0, 2, 8, 1, 9, 1, 1, 3, 1, 8), -1,  11,  125), ## 220
                          LS_jj_me( LS_jj_qn(0, 2, 8, 1, 9, 2, 2, 0, 0, 9), -1,  54,  125), ## 221
                          LS_jj_me( LS_jj_qn(0, 2, 8, 1, 9, 2, 0, 4, 2, 5), -1,  49,  125), ## 222
                          LS_jj_me( LS_jj_qn(0, 2, 8, 1, 9, 3, 1, 3, 1, 8), -1,  11,  125), ## 224
                          LS_jj_me( LS_jj_qn(0, 0, 8, 3, 9, 1, 1, 3, 1, 8), -1,   1,    6), ## 225
                          LS_jj_me( LS_jj_qn(0, 0, 8, 3, 9, 2, 0, 4, 0, 9), -1,   2,    3), ## 228
                          LS_jj_me( LS_jj_qn(0, 0, 8, 3, 9, 3, 1, 3, 1, 8),  1,   1,    6), ## 229
                          LS_jj_me( LS_jj_qn(0, 0, 8, 1, 9, 1, 1, 3, 1, 8),  1,   1,    3), ## 230
                          LS_jj_me( LS_jj_qn(0, 0, 8, 1, 9, 2, 0, 4, 0, 9), -1,   1,    3), ## 233
                          LS_jj_me( LS_jj_qn(0, 0, 8, 1, 9, 3, 1, 3, 1, 8), -1,   1,    3), ## 234
                          LS_jj_me( LS_jj_qn(0, 2,10, 1, 9, 1, 1, 3, 1, 8), -1,  24,  125), ## 235
                          LS_jj_me( LS_jj_qn(0, 2,10, 1, 9, 2, 2, 0, 0, 9), -1,  11,  125), ## 236
                          LS_jj_me( LS_jj_qn(0, 2,10, 1, 9, 2, 0, 4, 2, 5),  1,  66,  125), ## 237
                          LS_jj_me( LS_jj_qn(0, 2,10, 1, 9, 3, 1, 3, 1, 8), -1,  24,  125), ## 239
                          LS_jj_me( LS_jj_qn(0, 0, 8, 3,11, 1, 1, 3, 1, 8), -1,   6,   25), ## 240
                          LS_jj_me( LS_jj_qn(0, 0, 8, 3,11, 2, 0, 4, 0, 9), -1,  13,   25), ## 241
                          LS_jj_me( LS_jj_qn(0, 0, 8, 3,11, 3, 1, 3, 1, 8),  1,   6,   25), ## 242
                          LS_jj_me( LS_jj_qn(0, 2,10, 1,11, 1, 1, 3, 1, 8), -1,   1,    2), ## 243
                          LS_jj_me( LS_jj_qn(0, 2,10, 1,11, 3, 1, 3, 1, 8), -1,   1,    2), ## 245
                          LS_jj_me( LS_jj_qn(0, 0,12, 1,11, 1, 1, 3, 1, 8),  1,  13,   50), ## 246
                          LS_jj_me( LS_jj_qn(0, 0,12, 1,11, 2, 0, 4, 0, 9), -1,  12,   25), ## 247
                          LS_jj_me( LS_jj_qn(0, 0,12, 1,11, 3, 1, 3, 1, 8), -1,  13,   50), ## 248
                          LS_jj_me( LS_jj_qn(0, 0,12, 1,13, 2, 0, 4, 0, 9), -1,   1,    1) ]## 249
    ## d^6  shell;  198 ME in total                                                       
    const   LS_jj_d_6 = [ LS_jj_me( LS_jj_qn(0, 5, 0, 0, 0, 0, 2, 0, 3, 0),  1,  1,  10), ## 1
                          LS_jj_me( LS_jj_qn(0, 5, 0, 0, 0, 2, 2, 0, 3, 0),  1,  3,   5), ## 2
                          LS_jj_me( LS_jj_qn(0, 5, 0, 0, 0, 4, 2, 0, 3, 0),  1,  3,  10), ## 5
                          LS_jj_me( LS_jj_qn(0, 1, 0, 0, 0, 0, 2, 0, 3, 0), -1, 21, 250), ## 6
                          LS_jj_me( LS_jj_qn(0, 1, 0, 0, 0, 2, 2, 0, 3, 0),  1,  7, 125), ## 7
                          LS_jj_me( LS_jj_qn(0, 1, 0, 0, 0, 2, 0, 4, 1, 4),  1,  8,  25), ## 8
                          LS_jj_me( LS_jj_qn(0, 1, 0, 0, 0, 3, 1, 3, 0, 3), -1, 64, 125), ## 9
                          LS_jj_me( LS_jj_qn(0, 1, 0, 0, 0, 4, 2, 0, 3, 0), -1,  7, 250), ## 10
                          LS_jj_me( LS_jj_qn(0, 3, 2, 2, 0, 0, 2, 0, 3, 0), -1,  2,   5), ## 11
                          LS_jj_me( LS_jj_qn(0, 3, 2, 2, 0, 2, 2, 0, 3, 0), -1,  1,  15), ## 12
                          LS_jj_me( LS_jj_qn(0, 3, 2, 2, 0, 4, 2, 0, 3, 0),  1,  8,  15), ## 15
                          LS_jj_me( LS_jj_qn(0, 1, 2, 2, 0, 0, 2, 0, 3, 0), -1, 28, 125), ## 16
                          LS_jj_me( LS_jj_qn(0, 1, 2, 2, 0, 2, 2, 0, 3, 0),  1, 56, 375), ## 17
                          LS_jj_me( LS_jj_qn(0, 1, 2, 2, 0, 2, 0, 4, 1, 4),  1,  3,  25), ## 18
                          LS_jj_me( LS_jj_qn(0, 1, 2, 2, 0, 3, 1, 3, 0, 3),  1, 54, 125), ## 19
                          LS_jj_me( LS_jj_qn(0, 1, 2, 2, 0, 4, 2, 0, 3, 0), -1, 28, 375), ## 20
                          LS_jj_me( LS_jj_qn(0, 1, 4, 4, 0, 0, 2, 0, 3, 0), -1, 24, 125), ## 21
                          LS_jj_me( LS_jj_qn(0, 1, 4, 4, 0, 2, 2, 0, 3, 0),  1, 16, 125), ## 22
                          LS_jj_me( LS_jj_qn(0, 1, 4, 4, 0, 2, 0, 4, 1, 4), -1, 14,  25), ## 23
                          LS_jj_me( LS_jj_qn(0, 1, 4, 4, 0, 3, 1, 3, 0, 3), -1,  7, 125), ## 24
                          LS_jj_me( LS_jj_qn(0, 1, 4, 4, 0, 4, 2, 0, 3, 0), -1,  8, 125), ## 25
                          LS_jj_me( LS_jj_qn(0, 3, 2, 2, 2, 1, 1, 3, 2, 5),  1,  1,   3), ## 26
                          LS_jj_me( LS_jj_qn(0, 3, 2, 2, 2, 3, 1, 3, 2, 5),  1,  2,   3), ## 28
                          LS_jj_me( LS_jj_qn(0, 1, 2, 2, 2, 1, 1, 3, 2, 5),  1, 14,  75), ## 30
                          LS_jj_me( LS_jj_qn(0, 1, 2, 2, 2, 2, 0, 4, 1, 4),  1,  9,  25), ## 31
                          LS_jj_me( LS_jj_qn(0, 1, 2, 2, 2, 3, 1, 3, 2, 5), -1,  7,  75), ## 32
                          LS_jj_me( LS_jj_qn(0, 1, 2, 2, 2, 3, 1, 3, 0, 3),  1,  9,  25), ## 33
                          LS_jj_me( LS_jj_qn(0, 1, 4, 4, 2, 1, 1, 3, 2, 5),  1, 12,  25), ## 34
                          LS_jj_me( LS_jj_qn(0, 1, 4, 4, 2, 2, 0, 4, 1, 4), -1,  7,  50), ## 35
                          LS_jj_me( LS_jj_qn(0, 1, 4, 4, 2, 3, 1, 3, 2, 5), -1,  6,  25), ## 36
                          LS_jj_me( LS_jj_qn(0, 1, 4, 4, 2, 3, 1, 3, 0, 3), -1,  7,  50), ## 37
                          LS_jj_me( LS_jj_qn(0, 1, 4, 2, 2, 2, 0, 4, 1, 4), -1,  1,   2), ## 39
                          LS_jj_me( LS_jj_qn(0, 1, 4, 2, 2, 3, 1, 3, 0, 3),  1,  1,   2), ## 41
                          LS_jj_me( LS_jj_qn(0, 3, 2, 2, 4, 1, 1, 3, 2, 5),  1, 21, 125), ## 42
                          LS_jj_me( LS_jj_qn(0, 3, 2, 2, 4, 2, 2, 0, 1, 4),  1,112, 375), ## 43
                          LS_jj_me( LS_jj_qn(0, 3, 2, 2, 4, 2, 0, 4, 3, 0), -1,  6, 125), ## 44
                          LS_jj_me( LS_jj_qn(0, 3, 2, 2, 4, 3, 1, 3, 2, 5),  1, 42, 125), ## 47
                          LS_jj_me( LS_jj_qn(0, 3, 2, 2, 4, 4, 2, 0, 1, 4),  1, 56, 375), ## 49
                          LS_jj_me( LS_jj_qn(0, 1, 2, 2, 4, 1, 1, 3, 2, 5), -1,  6, 125), ## 50
                          LS_jj_me( LS_jj_qn(0, 1, 2, 2, 4, 2, 2, 0, 1, 4), -1, 32, 375), ## 51
                          LS_jj_me( LS_jj_qn(0, 1, 2, 2, 4, 2, 0, 4, 1, 4),  1,  3,  35), ## 53
                          LS_jj_me( LS_jj_qn(0, 1, 2, 2, 4, 2, 0, 4, 1, 8),  1,324, 875), ## 54
                          LS_jj_me( LS_jj_qn(0, 1, 2, 2, 4, 3, 1, 3, 2, 5),  1,  3, 125), ## 55
                          LS_jj_me( LS_jj_qn(0, 1, 2, 2, 4, 3, 1, 3, 0, 3),  1, 27, 125), ## 56
                          LS_jj_me( LS_jj_qn(0, 1, 2, 2, 4, 4, 2, 0, 1, 4),  1, 64, 375), ## 57
                          LS_jj_me( LS_jj_qn(0, 1, 4, 4, 4, 1, 1, 3, 2, 5),  1,  4,  25), ## 58
                          LS_jj_me( LS_jj_qn(0, 1, 4, 4, 4, 2, 2, 0, 1, 4),  1,  1,  25), ## 59
                          LS_jj_me( LS_jj_qn(0, 1, 4, 4, 4, 2, 0, 4, 1, 4), -1,  1,  14), ## 61
                          LS_jj_me( LS_jj_qn(0, 1, 4, 4, 4, 2, 0, 4, 1, 8),  1, 96, 175), ## 62
                          LS_jj_me( LS_jj_qn(0, 1, 4, 4, 4, 3, 1, 3, 2, 5), -1,  2,  25), ## 63
                          LS_jj_me( LS_jj_qn(0, 1, 4, 4, 4, 3, 1, 3, 0, 3), -1,  1,  50), ## 64
                          LS_jj_me( LS_jj_qn(0, 1, 4, 4, 4, 4, 2, 0, 1, 4), -1,  2,  25), ## 65
                          LS_jj_me( LS_jj_qn(0, 1, 4, 2, 4, 1, 1, 3, 2, 5),  1, 16,  75), ## 66
                          LS_jj_me( LS_jj_qn(0, 1, 4, 2, 4, 2, 2, 0, 1, 4), -1,  3,  25), ## 67
                          LS_jj_me( LS_jj_qn(0, 1, 4, 2, 4, 2, 0, 4, 1, 4), -1,  3,  14), ## 69
                          LS_jj_me( LS_jj_qn(0, 1, 4, 2, 4, 2, 0, 4, 1, 8), -1,  8, 175), ## 70
                          LS_jj_me( LS_jj_qn(0, 1, 4, 2, 4, 3, 1, 3, 2, 5), -1,  8,  75), ## 71
                          LS_jj_me( LS_jj_qn(0, 1, 4, 2, 4, 3, 1, 3, 0, 3),  1,  3,  50), ## 72
                          LS_jj_me( LS_jj_qn(0, 1, 4, 2, 4, 4, 2, 0, 1, 4),  1,  6,  25), ## 73
                          LS_jj_me( LS_jj_qn(0, 3, 4, 0, 4, 1, 1, 3, 2, 5), -1,  2,  25), ## 74
                          LS_jj_me( LS_jj_qn(0, 3, 4, 0, 4, 2, 2, 0, 1, 4),  1,  8,  25), ## 75
                          LS_jj_me( LS_jj_qn(0, 3, 4, 0, 4, 2, 0, 4, 3, 0),  1,  7,  25), ## 76
                          LS_jj_me( LS_jj_qn(0, 3, 4, 0, 4, 3, 1, 3, 2, 5), -1,  4,  25), ## 79
                          LS_jj_me( LS_jj_qn(0, 3, 4, 0, 4, 4, 2, 0, 1, 4),  1,  4,  25), ## 81
                          LS_jj_me( LS_jj_qn(0, 1, 4, 0, 4, 1, 1, 3, 2, 5), -1,  4,  25), ## 82
                          LS_jj_me( LS_jj_qn(0, 1, 4, 0, 4, 2, 2, 0, 1, 4), -1,  1,  25), ## 83
                          LS_jj_me( LS_jj_qn(0, 1, 4, 0, 4, 2, 0, 4, 1, 4), -1,  2,   7), ## 85
                          LS_jj_me( LS_jj_qn(0, 1, 4, 0, 4, 2, 0, 4, 1, 8),  1,  6, 175), ## 86
                          LS_jj_me( LS_jj_qn(0, 1, 4, 0, 4, 3, 1, 3, 2, 5),  1,  2,  25), ## 87
                          LS_jj_me( LS_jj_qn(0, 1, 4, 0, 4, 3, 1, 3, 0, 3), -1,  8,  25), ## 88
                          LS_jj_me( LS_jj_qn(0, 1, 4, 0, 4, 4, 2, 0, 1, 4),  1,  2,  25), ## 89
                          LS_jj_me( LS_jj_qn(0, 3, 6, 2, 4, 1, 1, 3, 2, 5),  1, 32, 375), ## 90
                          LS_jj_me( LS_jj_qn(0, 3, 6, 2, 4, 2, 2, 0, 1, 4), -1,  6, 125), ## 91
                          LS_jj_me( LS_jj_qn(0, 3, 6, 2, 4, 2, 0, 4, 3, 0),  1, 84, 125), ## 92
                          LS_jj_me( LS_jj_qn(0, 3, 6, 2, 4, 3, 1, 3, 2, 5),  1, 64, 375), ## 95
                          LS_jj_me( LS_jj_qn(0, 3, 6, 2, 4, 4, 2, 0, 1, 4), -1,  3, 125), ## 97
                          LS_jj_me( LS_jj_qn(0, 1, 6, 2, 4, 1, 1, 3, 2, 5),  1, 32, 375), ## 98
                          LS_jj_me( LS_jj_qn(0, 1, 6, 2, 4, 2, 2, 0, 1, 4), -1,  6, 125), ## 99
                          LS_jj_me( LS_jj_qn(0, 1, 6, 2, 4, 2, 0, 4, 1, 4),  1, 12,  35), ## 101
                          LS_jj_me( LS_jj_qn(0, 1, 6, 2, 4, 2, 0, 4, 1, 8),  1,  1, 875), ## 102
                          LS_jj_me( LS_jj_qn(0, 1, 6, 2, 4, 3, 1, 3, 2, 5), -1, 16, 375), ## 103
                          LS_jj_me( LS_jj_qn(0, 1, 6, 2, 4, 3, 1, 3, 0, 3), -1, 48, 125), ## 104
                          LS_jj_me( LS_jj_qn(0, 1, 6, 2, 4, 4, 2, 0, 1, 4),  1, 12, 125), ## 105
                          LS_jj_me( LS_jj_qn(0, 1, 4, 4, 6, 1, 1, 3, 2, 5),  1,  2,  25), ## 106
                          LS_jj_me( LS_jj_qn(0, 1, 4, 4, 6, 2, 0, 4, 1, 4), -1,  2, 175), ## 107
                          LS_jj_me( LS_jj_qn(0, 1, 4, 4, 6, 2, 0, 4, 1, 8),  1,  3,   7), ## 108
                          LS_jj_me( LS_jj_qn(0, 1, 4, 4, 6, 3, 1, 3, 2, 5), -1,  1,  25), ## 109
                          LS_jj_me( LS_jj_qn(0, 1, 4, 4, 6, 3, 1, 3, 0, 3),  1,  2, 175), ## 110
                          LS_jj_me( LS_jj_qn(0, 1, 4, 4, 6, 3, 1, 3, 0, 9), -1,  3,   7), ## 111
                          LS_jj_me( LS_jj_qn(0, 1, 4, 2, 6, 1, 1, 3, 2, 5),  1,  4,  21), ## 112
                          LS_jj_me( LS_jj_qn(0, 1, 4, 2, 6, 2, 0, 4, 1, 4), -1, 12,  49), ## 113
                          LS_jj_me( LS_jj_qn(0, 1, 4, 2, 6, 2, 0, 4, 1, 8),  1,  2,  49), ## 114
                          LS_jj_me( LS_jj_qn(0, 1, 4, 2, 6, 3, 1, 3, 2, 5), -1,  2,  21), ## 115
                          LS_jj_me( LS_jj_qn(0, 1, 4, 2, 6, 3, 1, 3, 0, 3),  1,  3,  49), ## 116
                          LS_jj_me( LS_jj_qn(0, 1, 4, 2, 6, 3, 1, 3, 0, 9),  1, 18,  49), ## 117
                          LS_jj_me( LS_jj_qn(0, 3, 6, 2, 6, 1, 1, 3, 2, 5),  1,  1,   3), ## 118
                          LS_jj_me( LS_jj_qn(0, 3, 6, 2, 6, 3, 1, 3, 2, 5),  1,  2,   3), ## 121
                          LS_jj_me( LS_jj_qn(0, 1, 6, 2, 6, 1, 1, 3, 2, 5),  1,  1, 300), ## 124
                          LS_jj_me( LS_jj_qn(0, 1, 6, 2, 6, 2, 0, 4, 1, 4),  1, 48, 175), ## 125
                          LS_jj_me( LS_jj_qn(0, 1, 6, 2, 6, 2, 0, 4, 1, 8),  1,  2,   7), ## 126
                          LS_jj_me( LS_jj_qn(0, 1, 6, 2, 6, 3, 1, 3, 2, 5), -1,  1, 600), ## 127
                          LS_jj_me( LS_jj_qn(0, 1, 6, 2, 6, 3, 1, 3, 0, 3), -1, 48, 175), ## 128
                          LS_jj_me( LS_jj_qn(0, 1, 6, 2, 6, 3, 1, 3, 0, 9),  1,  9,  56), ## 129
                          LS_jj_me( LS_jj_qn(0, 1, 6, 0, 6, 1, 1, 3, 2, 5), -1,  1,   5), ## 130
                          LS_jj_me( LS_jj_qn(0, 1, 6, 0, 6, 2, 0, 4, 1, 4),  1,  1,  35), ## 131
                          LS_jj_me( LS_jj_qn(0, 1, 6, 0, 6, 2, 0, 4, 1, 8),  1,  6,  35), ## 132
                          LS_jj_me( LS_jj_qn(0, 1, 6, 0, 6, 3, 1, 3, 2, 5),  1,  1,  10), ## 133
                          LS_jj_me( LS_jj_qn(0, 1, 6, 0, 6, 3, 1, 3, 0, 3),  1, 16,  35), ## 134
                          LS_jj_me( LS_jj_qn(0, 1, 6, 0, 6, 3, 1, 3, 0, 9),  1,  3,  70), ## 135
                          LS_jj_me( LS_jj_qn(0, 1, 8, 2, 6, 1, 1, 3, 2, 5),  1, 27, 140), ## 136
                          LS_jj_me( LS_jj_qn(0, 1, 8, 2, 6, 2, 0, 4, 1, 4),  1,108, 245), ## 137
                          LS_jj_me( LS_jj_qn(0, 1, 8, 2, 6, 2, 0, 4, 1, 8), -1, 18, 245), ## 138
                          LS_jj_me( LS_jj_qn(0, 1, 8, 2, 6, 3, 1, 3, 2, 5), -1, 27, 280), ## 139
                          LS_jj_me( LS_jj_qn(0, 1, 8, 2, 6, 3, 1, 3, 0, 3),  1, 48, 245), ## 140
                          LS_jj_me( LS_jj_qn(0, 1, 8, 2, 6, 3, 1, 3, 0, 9), -1,  1,1960), ## 141
                          LS_jj_me( LS_jj_qn(0, 1, 4, 4, 8, 1, 1, 3, 2, 5),  1,  6, 125), ## 142
                          LS_jj_me( LS_jj_qn(0, 1, 4, 4, 8, 2, 2, 0, 1, 8),  1, 12, 125), ## 143
                          LS_jj_me( LS_jj_qn(0, 1, 4, 4, 8, 2, 0, 4, 1, 4),  1,  2, 175), ## 144
                          LS_jj_me( LS_jj_qn(0, 1, 4, 4, 8, 2, 0, 4, 1, 8),  1, 33, 175), ## 145
                          LS_jj_me( LS_jj_qn(0, 1, 4, 4, 8, 3, 1, 3, 2, 5), -1,  3, 125), ## 146
                          LS_jj_me( LS_jj_qn(0, 1, 4, 4, 8, 3, 1, 3, 0, 9), -1, 11,  25), ## 147
                          LS_jj_me( LS_jj_qn(0, 1, 4, 4, 8, 4, 2, 0, 1, 8), -1, 24, 125), ## 148
                          LS_jj_me( LS_jj_qn(0, 3, 6, 2, 8, 1, 1, 3, 2, 5),  1,  1,  15), ## 149
                          LS_jj_me( LS_jj_qn(0, 3, 6, 2, 8, 2, 2, 0, 1, 8),  1,  8,  15), ## 150
                          LS_jj_me( LS_jj_qn(0, 3, 6, 2, 8, 3, 1, 3, 2, 5),  1,  2,  15), ## 153
                          LS_jj_me( LS_jj_qn(0, 3, 6, 2, 8, 4, 2, 0, 1, 8),  1,  4,  15), ## 155
                          LS_jj_me( LS_jj_qn(0, 1, 6, 2, 8, 1, 1, 3, 2, 5),  1,289,1500), ## 156
                          LS_jj_me( LS_jj_qn(0, 1, 6, 2, 8, 2, 2, 0, 1, 8), -1, 32, 375), ## 157
                          LS_jj_me( LS_jj_qn(0, 1, 6, 2, 8, 2, 0, 4, 1, 4),  1,  4, 175), ## 158
                          LS_jj_me( LS_jj_qn(0, 1, 6, 2, 8, 2, 0, 4, 1, 8),  1, 66, 175), ## 159
                          LS_jj_me( LS_jj_qn(0, 1, 6, 2, 8, 3, 1, 3, 2, 5), -1,289,3000), ## 160
                          LS_jj_me( LS_jj_qn(0, 1, 6, 2, 8, 3, 1, 3, 0, 9),  1, 11, 200), ## 161
                          LS_jj_me( LS_jj_qn(0, 1, 6, 2, 8, 4, 2, 0, 1, 8),  1, 64, 375), ## 162
                          LS_jj_me( LS_jj_qn(0, 1, 8, 2, 8, 1, 1, 3, 2, 5),  1, 81, 500), ## 163
                          LS_jj_me( LS_jj_qn(0, 1, 8, 2, 8, 2, 2, 0, 1, 8),  1,  8, 125), ## 164
                          LS_jj_me( LS_jj_qn(0, 1, 8, 2, 8, 2, 0, 4, 1, 4),  1, 48, 175), ## 165
                          LS_jj_me( LS_jj_qn(0, 1, 8, 2, 8, 2, 0, 4, 1, 8), -1, 22, 175), ## 166
                          LS_jj_me( LS_jj_qn(0, 1, 8, 2, 8, 3, 1, 3, 2, 5), -1, 81,1000), ## 167
                          LS_jj_me( LS_jj_qn(0, 1, 8, 2, 8, 3, 1, 3, 0, 9),  1, 33, 200), ## 168
                          LS_jj_me( LS_jj_qn(0, 1, 8, 2, 8, 4, 2, 0, 1, 8), -1, 16, 125), ## 169
                          LS_jj_me( LS_jj_qn(0, 3, 8, 0, 8, 1, 1, 3, 2, 5), -1,  4,  15), ## 170
                          LS_jj_me( LS_jj_qn(0, 3, 8, 0, 8, 2, 2, 0, 1, 8),  1,  2,  15), ## 171
                          LS_jj_me( LS_jj_qn(0, 3, 8, 0, 8, 3, 1, 3, 2, 5), -1,  8,  15), ## 174
                          LS_jj_me( LS_jj_qn(0, 3, 8, 0, 8, 4, 2, 0, 1, 8),  1,  1,  15), ## 176
                          LS_jj_me( LS_jj_qn(0, 1, 8, 0, 8, 1, 1, 3, 2, 5), -1, 11, 375), ## 177
                          LS_jj_me( LS_jj_qn(0, 1, 8, 0, 8, 2, 2, 0, 1, 8), -1, 22, 375), ## 178
                          LS_jj_me( LS_jj_qn(0, 1, 8, 0, 8, 2, 0, 4, 1, 4),  1, 99, 175), ## 179
                          LS_jj_me( LS_jj_qn(0, 1, 8, 0, 8, 2, 0, 4, 1, 8), -1,  6, 175), ## 180
                          LS_jj_me( LS_jj_qn(0, 1, 8, 0, 8, 3, 1, 3, 2, 5),  1, 11, 750), ## 181
                          LS_jj_me( LS_jj_qn(0, 1, 8, 0, 8, 3, 1, 3, 0, 9), -1,  9,  50), ## 182
                          LS_jj_me( LS_jj_qn(0, 1, 8, 0, 8, 4, 2, 0, 1, 8),  1, 44, 375), ## 183
                          LS_jj_me( LS_jj_qn(0, 1,10, 2, 8, 1, 1, 3, 2, 5), -1, 88, 375), ## 184
                          LS_jj_me( LS_jj_qn(0, 1,10, 2, 8, 2, 2, 0, 1, 8),  1, 11, 375), ## 185
                          LS_jj_me( LS_jj_qn(0, 1,10, 2, 8, 2, 0, 4, 1, 4),  1, 22, 175), ## 186
                          LS_jj_me( LS_jj_qn(0, 1,10, 2, 8, 2, 0, 4, 1, 8),  1, 48, 175), ## 187
                          LS_jj_me( LS_jj_qn(0, 1,10, 2, 8, 3, 1, 3, 2, 5),  1, 44, 375), ## 188
                          LS_jj_me( LS_jj_qn(0, 1,10, 2, 8, 3, 1, 3, 0, 9),  1,  4,  25), ## 189
                          LS_jj_me( LS_jj_qn(0, 1,10, 2, 8, 4, 2, 0, 1, 8), -1, 22, 375), ## 190
                          LS_jj_me( LS_jj_qn(0, 1, 8, 2,10, 2, 0, 4, 1, 8), -1,  9,  25), ## 191
                          LS_jj_me( LS_jj_qn(0, 1, 8, 2,10, 3, 1, 3, 0, 9),  1, 16,  25), ## 192
                          LS_jj_me( LS_jj_qn(0, 1,10, 2,10, 2, 0, 4, 1, 8),  1, 16,  25), ## 193
                          LS_jj_me( LS_jj_qn(0, 1,10, 2,10, 3, 1, 3, 0, 9),  1,  9,  25), ## 194
                          LS_jj_me( LS_jj_qn(0, 1,10, 2,12, 2, 0, 4, 1, 8),  1,  2,   5), ## 195
                          LS_jj_me( LS_jj_qn(0, 1,10, 2,12, 3, 1, 3, 0, 9),  1,  3,   5), ## 196
                          LS_jj_me( LS_jj_qn(0, 1,12, 0,12, 2, 0, 4, 1, 8), -1,  3,   5), ## 197
                          LS_jj_me( LS_jj_qn(0, 1,12, 0,12, 3, 1, 3, 0, 9),  1,  2,   5) ]## 198
    ## d^7  shell;  73 ME in total                                                      
    const   LS_jj_d_7 = [ LS_jj_me( LS_jj_qn(0, 2, 2, 3, 1, 2, 0, 4, 2, 5), -1,  8, 15), ## 1
                          LS_jj_me( LS_jj_qn(0, 2, 2, 3, 1, 3, 1, 3, 1, 4), -1,  7, 15), ## 2
                          LS_jj_me( LS_jj_qn(0, 2, 2, 1, 1, 2, 0, 4, 2, 5),  1,  7, 15), ## 3
                          LS_jj_me( LS_jj_qn(0, 2, 2, 1, 1, 3, 1, 3, 1, 4), -1,  8, 15), ## 4
                          LS_jj_me( LS_jj_qn(0, 2, 2, 3, 3, 1, 1, 3, 3, 0),  1,  8,125), ## 5
                          LS_jj_me( LS_jj_qn(0, 2, 2, 3, 3, 2, 0, 4, 2, 5), -1, 28, 75), ## 6
                          LS_jj_me( LS_jj_qn(0, 2, 2, 3, 3, 3, 1, 3, 3, 0), -1,  8,375), ## 7
                          LS_jj_me( LS_jj_qn(0, 2, 2, 3, 3, 3, 1, 3, 1, 4), -1, 28, 75), ## 8
                          LS_jj_me( LS_jj_qn(0, 2, 2, 3, 3, 4, 2, 0, 0, 3), -1, 21,125), ## 9
                          LS_jj_me( LS_jj_qn(0, 2, 2, 1, 3, 1, 1, 3, 3, 0), -1,  7, 50), ## 10
                          LS_jj_me( LS_jj_qn(0, 2, 2, 1, 3, 2, 0, 4, 2, 5),  1,  4, 15), ## 11
                          LS_jj_me( LS_jj_qn(0, 2, 2, 1, 3, 3, 1, 3, 3, 0),  1,  7,150), ## 12
                          LS_jj_me( LS_jj_qn(0, 2, 2, 1, 3, 3, 1, 3, 1, 4), -1,  1, 15), ## 13
                          LS_jj_me( LS_jj_qn(0, 2, 2, 1, 3, 4, 2, 0, 0, 3), -1, 12, 25), ## 14
                          LS_jj_me( LS_jj_qn(0, 4, 4, 1, 3, 1, 1, 3, 3, 0),  1,  1,  4), ## 15
                          LS_jj_me( LS_jj_qn(0, 4, 4, 1, 3, 3, 1, 3, 3, 0),  1,  3,  4), ## 17
                          LS_jj_me( LS_jj_qn(0, 2, 4, 1, 3, 1, 1, 3, 3, 0),  1, 21,100), ## 20
                          LS_jj_me( LS_jj_qn(0, 2, 4, 1, 3, 3, 1, 3, 3, 0), -1,  7,100), ## 22
                          LS_jj_me( LS_jj_qn(0, 2, 4, 1, 3, 3, 1, 3, 1, 4),  1,  2,  5), ## 23
                          LS_jj_me( LS_jj_qn(0, 2, 4, 1, 3, 4, 2, 0, 0, 3), -1,  8, 25), ## 24
                          LS_jj_me( LS_jj_qn(0, 2, 6, 3, 3, 1, 1, 3, 3, 0),  1, 42,125), ## 25
                          LS_jj_me( LS_jj_qn(0, 2, 6, 3, 3, 2, 0, 4, 2, 5),  1,  9, 25), ## 26
                          LS_jj_me( LS_jj_qn(0, 2, 6, 3, 3, 3, 1, 3, 3, 0), -1, 14,125), ## 27
                          LS_jj_me( LS_jj_qn(0, 2, 6, 3, 3, 3, 1, 3, 1, 4), -1,  4, 25), ## 28
                          LS_jj_me( LS_jj_qn(0, 2, 6, 3, 3, 4, 2, 0, 0, 3),  1,  4,125), ## 29
                          LS_jj_me( LS_jj_qn(0, 2, 2, 3, 5, 2, 2, 0, 2, 5),  1, 24,125), ## 30
                          LS_jj_me( LS_jj_qn(0, 2, 2, 3, 5, 3, 1, 3, 1, 4), -1,  1, 25), ## 32
                          LS_jj_me( LS_jj_qn(0, 2, 2, 3, 5, 3, 1, 3, 1, 8), -1, 72,125), ## 33
                          LS_jj_me( LS_jj_qn(0, 2, 2, 3, 5, 4, 2, 0, 2, 5), -1, 24,125), ## 34
                          LS_jj_me( LS_jj_qn(0, 4, 4, 1, 5, 2, 2, 0, 2, 5),  1,  1,  2), ## 35
                          LS_jj_me( LS_jj_qn(0, 4, 4, 1, 5, 4, 2, 0, 2, 5),  1,  1,  2), ## 39
                          LS_jj_me( LS_jj_qn(0, 2, 4, 1, 5, 2, 2, 0, 2, 5), -1,  7,150), ## 40
                          LS_jj_me( LS_jj_qn(0, 2, 4, 1, 5, 2, 0, 4, 2, 5),  1,  4, 15), ## 41
                          LS_jj_me( LS_jj_qn(0, 2, 4, 1, 5, 3, 1, 3, 1, 4),  1, 16, 35), ## 42
                          LS_jj_me( LS_jj_qn(0, 2, 4, 1, 5, 3, 1, 3, 1, 8), -1, 32,175), ## 43
                          LS_jj_me( LS_jj_qn(0, 2, 4, 1, 5, 4, 2, 0, 2, 5),  1,  7,150), ## 44
                          LS_jj_me( LS_jj_qn(0, 2, 6, 3, 5, 2, 2, 0, 2, 5),  1, 28,375), ## 45
                          LS_jj_me( LS_jj_qn(0, 2, 6, 3, 5, 2, 0, 4, 2, 5),  1,  2,  3), ## 46
                          LS_jj_me( LS_jj_qn(0, 2, 6, 3, 5, 3, 1, 3, 1, 4), -1,  8,175), ## 47
                          LS_jj_me( LS_jj_qn(0, 2, 6, 3, 5, 3, 1, 3, 1, 8),  1,121,875), ## 48
                          LS_jj_me( LS_jj_qn(0, 2, 6, 3, 5, 4, 2, 0, 2, 5), -1, 28,375), ## 49
                          LS_jj_me( LS_jj_qn(0, 2, 6, 1, 5, 2, 2, 0, 2, 5),  1, 14, 75), ## 50
                          LS_jj_me( LS_jj_qn(0, 2, 6, 1, 5, 2, 0, 4, 2, 5), -1,  1, 15), ## 51
                          LS_jj_me( LS_jj_qn(0, 2, 6, 1, 5, 3, 1, 3, 1, 4),  1, 16, 35), ## 52
                          LS_jj_me( LS_jj_qn(0, 2, 6, 1, 5, 3, 1, 3, 1, 8),  1, 18,175), ## 53
                          LS_jj_me( LS_jj_qn(0, 2, 6, 1, 5, 4, 2, 0, 2, 5), -1, 14, 75), ## 54
                          LS_jj_me( LS_jj_qn(0, 2, 6, 3, 7, 2, 0, 4, 2, 5),  1,  1,  5), ## 55
                          LS_jj_me( LS_jj_qn(0, 2, 6, 3, 7, 3, 1, 3, 1, 4), -1,  3, 35), ## 56
                          LS_jj_me( LS_jj_qn(0, 2, 6, 3, 7, 3, 1, 3, 1, 8),  1,  5,  7), ## 57
                          LS_jj_me( LS_jj_qn(0, 2, 6, 1, 7, 2, 0, 4, 2, 5), -1,  3, 10), ## 58
                          LS_jj_me( LS_jj_qn(0, 2, 6, 1, 7, 3, 1, 3, 1, 4),  1,121,280), ## 59
                          LS_jj_me( LS_jj_qn(0, 2, 6, 1, 7, 3, 1, 3, 1, 8),  1, 15, 56), ## 60
                          LS_jj_me( LS_jj_qn(0, 2, 8, 1, 7, 2, 0, 4, 2, 5), -1,  1,  2), ## 61
                          LS_jj_me( LS_jj_qn(0, 2, 8, 1, 7, 3, 1, 3, 1, 4), -1, 27, 56), ## 62
                          LS_jj_me( LS_jj_qn(0, 2, 8, 1, 7, 3, 1, 3, 1, 8),  1,  1, 56), ## 63
                          LS_jj_me( LS_jj_qn(0, 2, 6, 3, 9, 2, 0, 4, 2, 5),  1,  2, 25), ## 64
                          LS_jj_me( LS_jj_qn(0, 2, 6, 3, 9, 3, 1, 3, 1, 8),  1, 11, 25), ## 65
                          LS_jj_me( LS_jj_qn(0, 2, 6, 3, 9, 4, 2, 0, 0, 9), -1, 12, 25), ## 66
                          LS_jj_me( LS_jj_qn(0, 2, 8, 1, 9, 2, 0, 4, 2, 5), -1, 49,125), ## 67
                          LS_jj_me( LS_jj_qn(0, 2, 8, 1, 9, 3, 1, 3, 1, 8), -1, 22,125), ## 68
                          LS_jj_me( LS_jj_qn(0, 2, 8, 1, 9, 4, 2, 0, 0, 9), -1, 54,125), ## 69
                          LS_jj_me( LS_jj_qn(0, 2,10, 1, 9, 2, 0, 4, 2, 5),  1, 66,125), ## 70
                          LS_jj_me( LS_jj_qn(0, 2,10, 1, 9, 3, 1, 3, 1, 8), -1, 48,125), ## 71
                          LS_jj_me( LS_jj_qn(0, 2,10, 1, 9, 4, 2, 0, 0, 9), -1, 11,125), ## 72
                          LS_jj_me( LS_jj_qn(0, 2,10, 1,11, 3, 1, 3, 1, 8), -1,  1,  1) ]## 73
    ## d^8  shell;  19 ME in total                                                     
    const   LS_jj_d_8 = [ LS_jj_me( LS_jj_qn(0, 5, 0, 0, 0, 2, 2, 0, 3, 0),  1, 2,  5), ## 1
                          LS_jj_me( LS_jj_qn(0, 5, 0, 0, 0, 4, 2, 0, 3, 0),  1, 3,  5), ## 2
                          LS_jj_me( LS_jj_qn(0, 3, 2, 2, 0, 2, 2, 0, 3, 0), -1, 3,  5), ## 3
                          LS_jj_me( LS_jj_qn(0, 3, 2, 2, 0, 4, 2, 0, 3, 0),  1, 2,  5), ## 4
                          LS_jj_me( LS_jj_qn(0, 3, 2, 2, 2, 3, 1, 3, 2, 5),  1, 1,  1), ## 5
                          LS_jj_me( LS_jj_qn(0, 3, 2, 2, 4, 2, 0, 4, 3, 0), -1, 6,125), ## 6
                          LS_jj_me( LS_jj_qn(0, 3, 2, 2, 4, 3, 1, 3, 2, 5),  1,63,125), ## 7
                          LS_jj_me( LS_jj_qn(0, 3, 2, 2, 4, 4, 2, 0, 1, 4),  1,56,125), ## 8
                          LS_jj_me( LS_jj_qn(0, 3, 4, 0, 4, 2, 0, 4, 3, 0),  1, 7, 25), ## 9
                          LS_jj_me( LS_jj_qn(0, 3, 4, 0, 4, 3, 1, 3, 2, 5), -1, 6, 25), ## 10
                          LS_jj_me( LS_jj_qn(0, 3, 4, 0, 4, 4, 2, 0, 1, 4),  1,12, 25), ## 11
                          LS_jj_me( LS_jj_qn(0, 3, 6, 2, 4, 2, 0, 4, 3, 0),  1,84,125), ## 12
                          LS_jj_me( LS_jj_qn(0, 3, 6, 2, 4, 3, 1, 3, 2, 5),  1,32,125), ## 13
                          LS_jj_me( LS_jj_qn(0, 3, 6, 2, 4, 4, 2, 0, 1, 4), -1, 9,125), ## 14
                          LS_jj_me( LS_jj_qn(0, 3, 6, 2, 6, 3, 1, 3, 2, 5),  1, 1,  1), ## 15
                          LS_jj_me( LS_jj_qn(0, 3, 6, 2, 8, 3, 1, 3, 2, 5),  1, 1,  5), ## 16
                          LS_jj_me( LS_jj_qn(0, 3, 6, 2, 8, 4, 2, 0, 1, 8),  1, 4,  5), ## 17
                          LS_jj_me( LS_jj_qn(0, 3, 8, 0, 8, 3, 1, 3, 2, 5), -1, 4,  5), ## 18
                          LS_jj_me( LS_jj_qn(0, 3, 8, 0, 8, 4, 2, 0, 1, 8),  1, 1,  5) ]## 19
    ## d^9  shell;  2 ME in total                                                     
    const   LS_jj_d_9 = [ LS_jj_me( LS_jj_qn(0, 4, 4, 1, 3, 3, 1, 3, 3, 0),  1, 1, 1),  ## 1
                          LS_jj_me( LS_jj_qn(0, 4, 4, 1, 5, 4, 2, 0, 2, 5),  1, 1, 1) ] ## 2
    ## d^10  shell;  1 ME in total  LS_jj_qn(                                         
    const   LS_jj_d_10= [ LS_jj_me( LS_jj_qn(0, 5, 0, 0, 0, 4, 2, 0, 3, 0),  1, 1, 1) ] ## 1
                          
                              
    # Include coupling matrix elements for open f-shell
    ## include("../src/module-LSjj-inc.jl")
    
end # module			
