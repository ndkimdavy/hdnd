include("../../interface/isolver.jl")
include("../../interface/iopacity.jl")
include("../constant.jl")
include("../reject.jl")
include("../riemann.jl")
include("rmodel1d.jl")

using Base.Threads

#===================================================================
                            Radiation 1D
                  Multigroup radiation global update
                   (Er, Fr)[n] -> (Er, Fr, Pr)[n+1]
===================================================================#

struct Radiation1DInput{T} <: ISolverInput
    dt::T
    nx::Int
    dx::T
    domain::Domain1D
    mode::String
    tolsim::T
    gE::Array{T,1}
    opacity::IOpacity{T}
    nspectral::Int
    nangular::Int
    rho::Array{T,3}
    Th::Array{T,3}
    Er::Array{T,4}
    Fr::Array{T,5}
    wall::String
    mclosure::String
    hpclosure::Array{T,1}
    fclosure::String
    hprmodel::Array{T,1}
    rsolver::String
    retry::Int
    reject::Reject
end

struct Radiation1DOutput{T} <: ISolverOutput
    # Radiation energy: (nxlocal, ny, nz, ng)
    Er::Array{T,4}

    # Radiation flux: (nxlocal, ny, nz, ng, ndim)
    Fr::Array{T,5}

    # Radiation pressure: (nxlocal, ny, nz, ng, ndim, ndim)
    Pr::Array{T,6}

    # Radiation regime: (nxlocal, ny, nz, ng)
    regime::Array{Union{Nothing,Symbol},4}
end

struct Radiation1DSolver{Tin<:ISolverInput,Tout<:ISolverOutput} <: ISolver{Tin,Tout} end

function abort(solver::Radiation1DSolver{Radiation1DInput{T},Radiation1DOutput{T}},
    msg::String,
    dt::T,
    nx::Int, dx::T, domain::Domain1D,
    jx::Int,
    l::Int,
    rhoj::T, Thj::T,
    ceffjl::T, Erjl::T, Frjl::T, Prjl::T,
    fred::T, fedd::T) where {T<:AbstractFloat}

    log("ERROR: Radiation1D")
    log("solver  = $(typeof(solver))")
    log("dt      = $(dt)")
    log("nx      = $(nx)")
    log("dx      = $(dx)")
    log("nxlocal = $(domain.nxlocal)")
    log("rank    = $(domain.rank)")
    log("jx      = $(jx)")
    log("l       = $(l)")
    log("rhoj    = $(rhoj)")
    log("Thj     = $(Thj)")
    log("ceffjl  = $(ceffjl)")
    log("Erjl    = $(Erjl)")
    log("Frjl    = $(Frjl)")
    log("Prjl    = $(Prjl)")
    log("fred    = $(fred)")
    log("fedd    = $(fedd)")

    error("$(msg)")
end

function kernel(solver::Radiation1DSolver{Radiation1DInput{T},Radiation1DOutput{T}}, input::Radiation1DInput{T})::Radiation1DOutput{T} where {T<:AbstractFloat}

    dt = input.dt
    nx = input.nx
    dx = input.dx
    domain = input.domain
    nxlocal = domain.nxlocal
    mode = input.mode
    tolsim = input.tolsim
    gE = input.gE
    opacity = input.opacity
    nspectral = input.nspectral
    nangular = input.nangular
    rho = input.rho
    Th = input.Th
    Er = input.Er
    Fr = input.Fr
    wall = input.wall
    mclosure = input.mclosure
    hpclosure = input.hpclosure
    fclosure = input.fclosure
    hprmodel = input.hprmodel
    rsolver = input.rsolver
    retry = input.retry
    reject = input.reject

    reset!(reject)

    ng = ngroup(opacity)
    dtdx = dt / dx

    ceff = zeros(T, nxlocal, 1, 1, ng)
    kR = zeros(T, nxlocal, 1, 1, ng)

    Ernew = copy(Er)
    Frnew = copy(Fr)
    Prnew = zeros(T, nxlocal, 1, 1, ng, 1, 1)
    regime = Array{Union{Nothing,Symbol},4}(nothing, nxlocal, 1, 1, ng)

    nthd = Threads.nthreads()
    blockJx = [(div((ithd-1)*nxlocal, nthd)+1):div(ithd*nxlocal, nthd) for ithd in 1:nthd]
    rejectJx = [Reject(false, "") for ithd in 1:nthd]

    Threads.@threads :static for ithd in 1:nthd
        _reject = rejectJx[ithd]

        for jxlocal in blockJx[ithd]
            rhoj = rho[jxlocal, 1, 1]
            Thj = Th[jxlocal, 1, 1]

            for l in 1:ng
                kRjl = kappaR(opacity, Thj, rhoj, l)

                if !isfinite(kRjl) || kRjl <= T(0)
                    reject!(_reject, "Radiation1D invalid | kRjl = $(kRjl)")
                    break
                end

                ceffjl = mode == "classic" ? T(c) : min(T(c), T(c) / (rhoj * kRjl * dx))

                if !isfinite(ceffjl) || ceffjl <= T(0)
                    reject!(_reject, "Radiation1D invalid | ceffjl = $(ceffjl)")
                    break
                end

                ceff[jxlocal, 1, 1, l] = ceffjl
                kR[jxlocal, 1, 1, l] = kRjl
            end

            if _reject.flag
                break
            end
        end
    end

    # Local radiation equation
    radiation = (Erjl, Frjl, rhoj, Thj, ceffjl, kRjl, jx, l, gradErjl, gradFrjl, _reject) -> begin
        # ndim = 1 : U = [Er_l, Frx_l]

        # Project ghost cell index into physical domain
        jx = clamp(jx, 1, nx)

        Frjl = corrector(solver, Frjl, -(T(1) - eps(T)) * ceffjl * Erjl, (T(1) - eps(T)) * ceffjl * Erjl, tolsim, true)

        if abs(Frjl) > (T(1) - eps(T)) * ceffjl * Erjl
            reject!(_reject, "Radiation1D invalid | Frjl = $(Frjl)")
            return zeros(T, 2), zeros(T, 2), T(0)
        end

        _input = RModel1DInput(
            l,
            mode,
            tolsim,
            gE,
            kRjl,
            ceffjl,
            nspectral,
            nangular,
            rhoj,
            Thj,
            Erjl,
            Frjl,
            gradErjl,
            gradFrjl,
            mclosure,
            hpclosure,
            fclosure,
            hprmodel,
            retry,
            _reject
        )

        _solver = RModel1DSolver{RModel1DInput{T},RModel1DOutput{T}}()
        output = kernel(_solver, _input)

        if _reject.flag
            return zeros(T, 2), zeros(T, 2), T(0)
        end

        Erjl = output.Erjl
        Frjl = output.Frjl
        Prjl = output.Prjl

        if !isfinite(Erjl) || Erjl <= T(0)
            reject!(_reject, "Radiation1D invalid | Erjl = $(Erjl)")
            return zeros(T, 2), zeros(T, 2), T(0)
        end

        if !isfinite(Frjl)
            reject!(_reject, "Radiation1D invalid | Frjl = $(Frjl)")
            return zeros(T, 2), zeros(T, 2), T(0)
        end

        if !isfinite(Prjl) || Prjl < T(0)
            reject!(_reject, "Radiation1D invalid | Prjl = $(Prjl)")
            return zeros(T, 2), zeros(T, 2), T(0)
        end

        # Conservative state
        U = T[
            Erjl,
            Frjl
        ]

        # Radiation flux in x-direction
        F = T[
            Frjl,
            ceffjl*ceffjl*Prjl
        ]

        return U, F, ceffjl
    end

    # Local Riemann flux
    riemann = (ErL, FrLx, rhoL, ThL, ceffL, kRL, jxL,
        ErR, FrRx, rhoR, ThR, ceffR, kRR, jxR, l, _reject) -> begin

        gradErjl = (ErR - ErL) / dx
        gradFrjl = (FrRx - FrLx) / dx

        UL, FL, ceffL = radiation(ErL, FrLx, rhoL, ThL, ceffL, kRL, jxL, l, gradErjl, gradFrjl, _reject)

        if _reject.flag
            return zeros(T, 2)
        end

        UR, FR, ceffR = radiation(ErR, FrRx, rhoR, ThR, ceffR, kRR, jxR, l, gradErjl, gradFrjl, _reject)

        if _reject.flag
            return zeros(T, 2)
        end

        if rsolver == "rusanov"
            smax = max(ceffL, ceffR)

            _input = RusanovInput(UL, UR, FL, FR, smax, _reject)
            _solver = RusanovSolver{T}()
            output = kernel(_solver, _input)

            return output.out

        elseif rsolver in ("hll", "hlle")
            smax = max(ceffL, ceffR)

            SL = -smax
            SR = smax

            if rsolver == "hll"
                _input = HLLInput(UL, UR, FL, FR, SL, SR, _reject)
                _solver = HLLSolver{T}()
                output = kernel(_solver, _input)

                return output.out

            else
                _input = HLLEInput(UL, UR, FL, FR, SL, SR, _reject)
                _solver = HLLESolver{T}()
                output = kernel(_solver, _input)

                return output.out
            end

        else
            abort(solver, "Radiation1D unknown rsolver", dt, nx, dx, domain, 0, l, T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0))
        end
    end

    #=
    Godunov update radiation1D
    Spatial direction: x
        - jxlocal  : current local cell index
        - jx       : current global cell index
        - jxL      : left global cell index
        - jxR      : right global cell index
        - jxGL     : left global ghost cell index,  jxGL = 0
        - jxGR     : right global ghost cell index, jxGR = nx + 1
        - sendleft : local left boundary state sent to left rank
        - sendright: local right boundary state sent to right rank
        - recvleft : ghost state received from left rank
        - recvright: ghost state received from right rank
        - wall     : boundary wall selector
                     "none"  -> no wall
                     "left"  -> wall at (ghost left | cell 1)
                     "right" -> wall at (cell nx | ghost right)

    (Er, Frx): left group state
    (Er, Frx): current group state
    (Er, Frx): right group state

    For each group l in 1D:
        [Er_l, Frx_l](t) -> [Er_l, Frx_l](t + dt)
    =#

    jxGL = 0
    jxGR = nx + 1

    # Flat MPI buffers: [rho, Th, Er1, Fr1, ..., Erng, Frng]
    sendleft = zeros(T, 2 + 2 * ng)
    sendright = zeros(T, 2 + 2 * ng)
    recvleft = zeros(T, 2 + 2 * ng)
    recvright = zeros(T, 2 + 2 * ng)

    sendleft[1] = rho[1, 1, 1]
    sendleft[2] = Th[1, 1, 1]

    sendright[1] = rho[nxlocal, 1, 1]
    sendright[2] = Th[nxlocal, 1, 1]

    for l in 1:ng
        sendleft[2*l+1] = Er[1, 1, 1, l]
        sendleft[2*l+2] = Fr[1, 1, 1, l, 1]

        sendright[2*l+1] = Er[nxlocal, 1, 1, l]
        sendright[2*l+2] = Fr[nxlocal, 1, 1, l, 1]
    end

    exchange!(domain, sendleft, sendright, recvleft, recvright)

    for _reject in rejectJx
        if _reject.flag
            reject!(reject, _reject.msg)
            return Radiation1DOutput(Ernew, Frnew, Prnew, regime)
        end
    end

    Threads.@threads :static for ithd in 1:nthd
        _reject = rejectJx[ithd]

        for jxlocal in blockJx[ithd]
            jx = domain.jxbegin + jxlocal - 1
            jxL = jx - 1
            jxR = jx + 1

            iswallL = wall == "left" && jxL == jxGL
            iswallR = wall == "right" && jxR == jxGR

            rhoC = rho[jxlocal, 1, 1]
            ThC = Th[jxlocal, 1, 1]

            if jxlocal == 1
                if hasleft(domain)
                    rhoL = recvleft[1]
                    ThL = recvleft[2]

                elseif iswallL
                    rhoL = rhoC
                    ThL = ThC

                else
                    rhoL = rhoC
                    ThL = ThC
                end

            else
                rhoL = rho[jxlocal - 1, 1, 1]
                ThL = Th[jxlocal - 1, 1, 1]
            end

            if jxlocal == nxlocal
                if hasright(domain)
                    rhoR = recvright[1]
                    ThR = recvright[2]

                elseif iswallR
                    rhoR = rhoC
                    ThR = ThC

                else
                    rhoR = rhoC
                    ThR = ThC
                end

            else
                rhoR = rho[jxlocal + 1, 1, 1]
                ThR = Th[jxlocal + 1, 1, 1]
            end

            for l in 1:ng
                ErC = Er[jxlocal, 1, 1, l]
                FrCx = Fr[jxlocal, 1, 1, l, 1]

                ceffC = ceff[jxlocal, 1, 1, l]
                kRC = kR[jxlocal, 1, 1, l]

                if jxlocal == 1
                    if hasleft(domain)
                        ErL = recvleft[2*l+1]
                        FrLx = recvleft[2*l+2]

                        kRL = kappaR(opacity, ThL, rhoL, l)

                        if !isfinite(kRL) || kRL <= T(0)
                            reject!(_reject, "Radiation1D invalid | kRL = $(kRL)")
                            break
                        end

                        ceffL = mode == "classic" ? T(c) : min(T(c), T(c) / (rhoL * kRL * dx))

                        if !isfinite(ceffL) || ceffL <= T(0)
                            reject!(_reject, "Radiation1D invalid | ceffL = $(ceffL)")
                            break
                        end

                    elseif iswallL
                        ErL = ErC
                        FrLx = FrCx
                        ceffL = ceffC
                        kRL = kRC

                    else
                        ErL = ErC
                        FrLx = FrCx
                        ceffL = ceffC
                        kRL = kRC
                    end

                else
                    ErL = Er[jxlocal - 1, 1, 1, l]
                    FrLx = Fr[jxlocal - 1, 1, 1, l, 1]
                    ceffL = ceff[jxlocal - 1, 1, 1, l]
                    kRL = kR[jxlocal - 1, 1, 1, l]
                end

                if jxlocal == nxlocal
                    if hasright(domain)
                        ErR = recvright[2*l+1]
                        FrRx = recvright[2*l+2]

                        kRR = kappaR(opacity, ThR, rhoR, l)

                        if !isfinite(kRR) || kRR <= T(0)
                            reject!(_reject, "Radiation1D invalid | kRR = $(kRR)")
                            break
                        end

                        ceffR = mode == "classic" ? T(c) : min(T(c), T(c) / (rhoR * kRR * dx))

                        if !isfinite(ceffR) || ceffR <= T(0)
                            reject!(_reject, "Radiation1D invalid | ceffR = $(ceffR)")
                            break
                        end

                    elseif iswallR
                        ErR = ErC
                        FrRx = FrCx
                        ceffR = ceffC
                        kRR = kRC

                    else
                        ErR = ErC
                        FrRx = FrCx
                        ceffR = ceffC
                        kRR = kRC
                    end

                else
                    ErR = Er[jxlocal + 1, 1, 1, l]
                    FrRx = Fr[jxlocal + 1, 1, 1, l, 1]
                    ceffR = ceff[jxlocal + 1, 1, 1, l]
                    kRR = kR[jxlocal + 1, 1, 1, l]
                end

                phiL = riemann(ErL, FrLx, rhoL, ThL, ceffL, kRL, jxL,
                    ErC, FrCx, rhoC, ThC, ceffC, kRC, jx, l, _reject)

                if _reject.flag
                    break
                end

                phiR = riemann(ErC, FrCx, rhoC, ThC, ceffC, kRC, jx,
                    ErR, FrRx, rhoR, ThR, ceffR, kRR, jxR, l, _reject)

                if _reject.flag
                    break
                end

                Erjl = ErC - dtdx * (phiR[1] - phiL[1])
                Frjl = FrCx - dtdx * (phiR[2] - phiL[2])

                if !isfinite(Erjl) || Erjl <= T(0)
                    reject!(_reject, "Radiation1D invalid | Erjl = $(Erjl)")
                    break
                end

                ceffjl = ceffC
                kRjl = kRC

                Frjl = corrector(solver, Frjl, -(T(1) - eps(T)) * ceffjl * Erjl, (T(1) - eps(T)) * ceffjl * Erjl, tolsim, true)

                if abs(Frjl) > (T(1) - eps(T)) * ceffjl * Erjl
                    reject!(_reject, "Radiation1D invalid | Frjl = $(Frjl)")
                    break
                end

                gradErjl = (ErR - ErL) / (T(2) * dx)
                gradFrjl = (FrRx - FrLx) / (T(2) * dx)

                _input = RModel1DInput(
                    l,
                    mode,
                    tolsim,
                    gE,
                    kRjl,
                    ceffjl,
                    nspectral,
                    nangular,
                    rhoC,
                    ThC,
                    Erjl,
                    Frjl,
                    gradErjl,
                    gradFrjl,
                    mclosure,
                    hpclosure,
                    fclosure,
                    hprmodel,
                    retry,
                    _reject
                )

                _solver = RModel1DSolver{RModel1DInput{T},RModel1DOutput{T}}()
                output = kernel(_solver, _input)

                if _reject.flag
                    break
                end

                Ernew[jxlocal, 1, 1, l] = output.Erjl
                Frnew[jxlocal, 1, 1, l, 1] = output.Frjl
                Prnew[jxlocal, 1, 1, l, 1, 1] = output.Prjl
                regime[jxlocal, 1, 1, l] = output.regime
            end

            if _reject.flag
                break
            end
        end
    end

    for _reject in rejectJx
        if _reject.flag
            reject!(reject, _reject.msg)
            return Radiation1DOutput(Ernew, Frnew, Prnew, regime)
        end
    end

    return Radiation1DOutput(Ernew, Frnew, Prnew, regime)
end
