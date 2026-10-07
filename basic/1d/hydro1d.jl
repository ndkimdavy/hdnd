include("../../interface/isolver.jl")
include("../reject.jl")
include("../riemann.jl")

using Base.Threads

#===================================================================
                              Hydro 1D
                        Euler global update
                (rho, mom, Eh)[n] -> (rho, mom, Eh)[n+1]
===================================================================#

struct Hydro1DInput{T} <: ISolverInput
    dt::T
    nx::Int
    dx::T
    domain::Domain1D
    gamma::T
    tolsim::T
    rho::Array{T,3}
    mom::Array{T,4}
    Eh::Array{T,3}
    wall::String
    hsolver::String
    retry::Int
    reject::Reject
end

struct Hydro1DOutput{T} <: ISolverOutput
    # Density: (nxlocal, ny, nz)
    rho::Array{T,3}

    # Momentum: (nxlocal, ny, nz, ndim)
    mom::Array{T,4}

    # Hydro energy: (nxlocal, ny, nz)
    Eh::Array{T,3}
end

struct Hydro1DSolver{Tin<:ISolverInput,Tout<:ISolverOutput} <: ISolver{Tin,Tout} end

function abort(solver::Hydro1DSolver{Hydro1DInput{T},Hydro1DOutput{T}},
    msg::String,
    dt::T,
    nx::Int, dx::T, domain::Domain1D,
    jx::Int,
    vx::T, cs::T, p::T, rhoj::T, momj::T, Ehj::T) where {T<:AbstractFloat}

    log("ERROR: Hydro1D")
    log("solver  = $(typeof(solver))")
    log("dt      = $(dt)")
    log("nx      = $(nx)")
    log("dx      = $(dx)")
    log("nxlocal = $(domain.nxlocal)")
    log("rank    = $(domain.rank)")
    log("jx      = $(jx)")
    log("vx      = $(vx)")
    log("cs      = $(cs)")
    log("p       = $(p)")
    log("rhoj    = $(rhoj)")
    log("momj    = $(momj)")
    log("Ehj     = $(Ehj)")

    error("$(msg)")
end

function kernel(solver::Hydro1DSolver{Hydro1DInput{T},Hydro1DOutput{T}}, input::Hydro1DInput{T})::Hydro1DOutput{T} where {T<:AbstractFloat}

    dt = input.dt
    nx = input.nx
    dx = input.dx
    domain = input.domain
    nxlocal = domain.nxlocal
    gamma = input.gamma
    tolsim = input.tolsim
    rho = input.rho
    mom = input.mom
    Eh = input.Eh
    wall = input.wall
    hsolver = input.hsolver
    retry = input.retry
    reject = input.reject

    reset!(reject)

    dtdx = dt / dx

    rhonew = copy(rho)
    momnew = copy(mom)
    Ehnew = copy(Eh)

    nthd = Threads.nthreads()
    blockJx = [div((ithd - 1) * nxlocal, nthd)+1:div(ithd * nxlocal, nthd) for ithd in 1:nthd]
    rejectJx = [Reject(false, "") for ithd in 1:nthd]

    # Local hydro Euler equation
    hydro = (rhoj, momj, Ehj, jx, _reject) -> begin
        # ndim = 1 : U = [rho, rho*vx, Eh]

        vx = momj / rhoj

        # Kinetic energy density
        Ek = T(0.5) * rhoj * vx^2

        # Pressure from ideal gas law
        p = (gamma - T(1)) * (Ehj - Ek)

        if !isfinite(p) || p <= T(0)
            reject!(_reject, "Hydro1D invalid | p = $(p)")
            return zeros(T, 3), zeros(T, 3), T(0), T(0), T(0)
        end

        # Conservative state
        U = T[
            rhoj,
            momj,
            Ehj
        ]

        # Euler flux in x-direction
        F = T[
            rhoj*vx,
            rhoj*vx*vx+p,
            (Ehj+p)*vx
        ]

        # Sound speed
        cs = sqrt(gamma * p / rhoj)

        if !isfinite(cs) || cs <= T(0)
            reject!(_reject, "Hydro1D invalid | cs = $(cs)")
            return zeros(T, 3), zeros(T, 3), T(0), T(0), T(0)
        end

        return U, F, vx, p, cs
    end

    # Local Riemann flux
    riemann = (rhoL, momL, EhL, jxL, rhoR, momR, EhR, jxR, _reject) -> begin
        UL, FL, vxL, pL, csL = hydro(rhoL, momL, EhL, jxL, _reject)

        if _reject.flag
            return zeros(T, 3)
        end

        UR, FR, vxR, pR, csR = hydro(rhoR, momR, EhR, jxR, _reject)

        if _reject.flag
            return zeros(T, 3)
        end

        if hsolver == "rusanov"
            smax = max(abs(vxL) + csL, abs(vxR) + csR)

            _input = RusanovInput(UL, UR, FL, FR, smax, _reject)
            _solver = RusanovSolver{T}()
            output = kernel(_solver, _input)

            return output.out

        elseif hsolver in ("hll", "hlle", "hllc")
            SL = min(vxL - csL, vxR - csR)
            SR = max(vxL + csL, vxR + csR)

            if hsolver == "hll"
                _input = HLLInput(UL, UR, FL, FR, SL, SR, _reject)
                _solver = HLLSolver{T}()
                output = kernel(_solver, _input)

                return output.out

            elseif hsolver == "hlle"
                _input = HLLEInput(UL, UR, FL, FR, SL, SR, _reject)
                _solver = HLLESolver{T}()
                output = kernel(_solver, _input)

                return output.out

            else
                SMden = rhoL * (SL - vxL) - rhoR * (SR - vxR)

                # HLLE fallback for singular HLLC denominator
                if abs(SMden) <= eps(T)
                    _input = HLLEInput(UL, UR, FL, FR, SL, SR, _reject)
                    _solver = HLLESolver{T}()
                    output = kernel(_solver, _input)

                    return output.out
                end

                SM = (pR - pL + rhoL * vxL * (SL - vxL) - rhoR * vxR * (SR - vxR)) / SMden

                # HLLE fallback for singular HLLC speeds
                if !isfinite(SM) || abs(SL - SM) <= eps(T) || abs(SR - SM) <= eps(T)
                    _input = HLLEInput(UL, UR, FL, FR, SL, SR, _reject)
                    _solver = HLLESolver{T}()
                    output = kernel(_solver, _input)

                    return output.out
                end

                EhSL = (SL - vxL) / (SL - SM) * (EhL + (SM - vxL) * (rhoL * SM + pL / (SL - vxL)))
                EhSR = (SR - vxR) / (SR - SM) * (EhR + (SM - vxR) * (rhoR * SM + pR / (SR - vxR)))

                USL = T[
                    rhoL*(SL-vxL)/(SL-SM),
                    rhoL*(SL-vxL)/(SL-SM)*SM,
                    EhSL
                ]

                USR = T[
                    rhoR*(SR-vxR)/(SR-SM),
                    rhoR*(SR-vxR)/(SR-SM)*SM,
                    EhSR
                ]

                # HLLE fallback for invalid HLLC states
                if any(x -> !isfinite(x), USL) || any(x -> !isfinite(x), USR)
                    _input = HLLEInput(UL, UR, FL, FR, SL, SR, _reject)
                    _solver = HLLESolver{T}()
                    output = kernel(_solver, _input)

                    return output.out
                end

                _input = HLLCInput(UL, UR, FL, FR, USL, USR, SL, SM, SR, _reject)
                _solver = HLLCSolver{T}()
                output = kernel(_solver, _input)

                return output.out
            end

        else
            abort(solver, "Hydro1D unknown hsolver", dt, nx, dx, domain, 0, T(0), T(0), T(0), T(0), T(0), T(0))
        end
    end

    #=
    Godunov update hydro1D
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

    (rho, mom, Eh): left state
    (rho, mom, Eh): current state
    (rho, mom, Eh): right state

    For each cell in 1D:
        [rho, mom, Eh](t) -> [rho, mom, Eh](t + dt)
    =#

    jxGL = 0
    jxGR = nx + 1

    sendleft = T[
        rho[1, 1, 1],
        mom[1, 1, 1, 1],
        Eh[1, 1, 1]
    ]

    sendright = T[
        rho[nxlocal, 1, 1],
        mom[nxlocal, 1, 1, 1],
        Eh[nxlocal, 1, 1]
    ]

    recvleft = zeros(T, 3)
    recvright = zeros(T, 3)

    exchange!(domain, sendleft, sendright, recvleft, recvright)

    Threads.@threads :static for ithd in 1:nthd
        _reject = rejectJx[ithd]

        for jxlocal in blockJx[ithd]
            jx = domain.jxbegin + jxlocal - 1
            jxL = jx - 1
            jxR = jx + 1

            rhoC = rho[jxlocal, 1, 1]
            momC = mom[jxlocal, 1, 1, 1]
            EhC = Eh[jxlocal, 1, 1]

            iswallL = wall == "left" && jxL == jxGL
            iswallR = wall == "right" && jxR == jxGR

            if jxlocal == 1
                if hasleft(domain)
                    rhoL = recvleft[1]
                    momL = recvleft[2]
                    EhL = recvleft[3]

                elseif iswallL
                    rhoL = rhoC
                    momL = -momC
                    EhL = EhC

                else
                    rhoL = rhoC
                    momL = momC
                    EhL = EhC
                end

            else
                rhoL = rho[jxlocal - 1, 1, 1]
                momL = mom[jxlocal - 1, 1, 1, 1]
                EhL = Eh[jxlocal - 1, 1, 1]
            end

            if jxlocal == nxlocal
                if hasright(domain)
                    rhoR = recvright[1]
                    momR = recvright[2]
                    EhR = recvright[3]

                elseif iswallR
                    rhoR = rhoC
                    momR = -momC
                    EhR = EhC

                else
                    rhoR = rhoC
                    momR = momC
                    EhR = EhC
                end

            else
                rhoR = rho[jxlocal + 1, 1, 1]
                momR = mom[jxlocal + 1, 1, 1, 1]
                EhR = Eh[jxlocal + 1, 1, 1]
            end

            phiL = riemann(rhoL, momL, EhL, jxL, rhoC, momC, EhC, jx, _reject)

            if _reject.flag
                break
            end

            phiR = riemann(rhoC, momC, EhC, jx, rhoR, momR, EhR, jxR, _reject)

            if _reject.flag
                break
            end

            rhonew[jxlocal, 1, 1] = rhoC - dtdx * (phiR[1] - phiL[1])
            momnew[jxlocal, 1, 1, 1] = momC - dtdx * (phiR[2] - phiL[2])
            Ehnew[jxlocal, 1, 1] = EhC - dtdx * (phiR[3] - phiL[3])

            if !isfinite(rhonew[jxlocal, 1, 1]) || rhonew[jxlocal, 1, 1] <= T(0)
                reject!(_reject, "Hydro1D invalid | rhonew[jxlocal, 1, 1] = $(rhonew[jxlocal, 1, 1])")
                break
            end

            if !isfinite(momnew[jxlocal, 1, 1, 1]) || !isfinite(Ehnew[jxlocal, 1, 1])
                reject!(_reject, "Hydro1D invalid | momnew[jxlocal, 1, 1, 1] = $(momnew[jxlocal, 1, 1, 1]) Ehnew[jxlocal, 1, 1] = $(Ehnew[jxlocal, 1, 1])")
                break
            end

            vxnew = momnew[jxlocal, 1, 1, 1] / rhonew[jxlocal, 1, 1]
            Eknew = T(0.5) * rhonew[jxlocal, 1, 1] * vxnew^2
            pnew = (gamma - T(1)) * (Ehnew[jxlocal, 1, 1] - Eknew)

            if !isfinite(pnew) || pnew <= T(0)
                reject!(_reject, "Hydro1D invalid | pnew = $(pnew)")
                break
            end
        end
    end

    for _reject in rejectJx
        if _reject.flag
            reject!(reject, _reject.msg)
            return Hydro1DOutput(rhonew, momnew, Ehnew)
        end
    end

    return Hydro1DOutput(rhonew, momnew, Ehnew)
end
