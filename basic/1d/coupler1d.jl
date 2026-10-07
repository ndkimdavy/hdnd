include("../../interface/isolver.jl")
include("../../interface/iopacity.jl")
include("../constant.jl")
include("../reject.jl")
include("domain1d.jl")
include("hydro1d.jl")
include("radiation1d.jl")
include("source1d.jl")
include("monitor1d.jl")

#===========================================================================================================================
                                                    Coupler 1D
                                            Hydro + Radiation + Coupling
        (rho, mom, Eh, Er, Fr)[n] -> (rho, mom, Eh)* -> (Er, Fr, Pr)* -> Source1D -> (rho, mom, Eh, Er, Fr, Pr)[n+1]
===========================================================================================================================#

struct Coupler1DInput{T} <: ISolverInput
    dx::T
    domain::Domain1D
    tsim::T
    courant::T
    mode::String
    gamma::T
    mu::T
    tolsim::T
    gE::Array{T,1}
    opacity::IOpacity{T}
    nspectral::Int
    nangular::Int
    rho0::T
    vx0::T
    Th0::T
    Er0::Array{T,1}
    Frx0::Array{T,1}
    hrho0::Array{T,1}
    hvx0::Array{T,1}
    hTh0::Array{T,1}
    hEr0::Array{T,1}
    hFrx0::Array{T,1}
    wall::String
    gmax::T
    mcool::String
    hpcool::Array{T,1}
    mclosure::String
    hpclosure::Array{T,1}
    fclosure::String
    hprmodel::Array{T,1}
    hsolver::String
    rsolver::String
    retry::Int
    reject::Reject
end

struct Coupler1DOutput{T} <: ISolverOutput end

struct Coupler1DSolver{Tin<:ISolverInput,Tout<:ISolverOutput} <: ISolver{Tin,Tout} end

function abort(solver::Coupler1DSolver{Coupler1DInput{T},Coupler1DOutput{T}},
    msg::String,
    it::Int, t::T, dt::T,
    dx::T, domain::Domain1D,
    jx::Int,
    l::Int,
    vx::T, cs::T, p::T, Thj::T, rhoj::T, momj::T, Ehj::T, Q0::T, Qx::T,
    ceffjl::T, Erjl::T, Frjl::T, Prjl::T, Trjl::T, cS0::T, Sx::T) where {T<:AbstractFloat}

    log("ERROR: Coupler1D")
    log("solver  = $(typeof(solver))")
    log("it      = $(it)")
    log("t       = $(t)")
    log("dt      = $(dt)")
    log("nx      = $(domain.nx)")
    log("dx      = $(dx)")
    log("nxlocal = $(domain.nxlocal)")
    log("rank    = $(domain.rank)")
    log("jx      = $(jx)")
    log("l       = $(l)")
    log("vx      = $(vx)")
    log("cs      = $(cs)")
    log("p       = $(p)")
    log("Thj     = $(Thj)")
    log("rhoj    = $(rhoj)")
    log("momj    = $(momj)")
    log("Ehj     = $(Ehj)")
    log("Q0      = $(Q0)")
    log("Qx      = $(Qx)")
    log("ceffjl  = $(ceffjl)")
    log("Erjl    = $(Erjl)")
    log("Frjl    = $(Frjl)")
    log("Prjl    = $(Prjl)")
    log("Trjl    = $(Trjl)")
    log("cS0     = $(cS0)")
    log("Sx      = $(Sx)")

    error("$(msg)")
end

function dispatch(solver::Coupler1DSolver{Coupler1DInput{T},Coupler1DOutput{T}}, input::Coupler1DInput{T})::Coupler1DOutput{T} where {T<:AbstractFloat}

    dx = input.dx
    domain = input.domain
    nx = domain.nx
    nxlocal = domain.nxlocal
    tsim = input.tsim
    courant = input.courant
    mode = input.mode
    gamma = input.gamma
    mu = input.mu
    tolsim = input.tolsim
    gE = input.gE
    opacity = input.opacity
    nspectral = input.nspectral
    nangular = input.nangular
    rho0 = input.rho0
    vx0 = input.vx0
    Th0 = input.Th0
    Er0 = input.Er0
    Frx0 = input.Frx0
    hrho0 = input.hrho0
    hvx0 = input.hvx0
    hTh0 = input.hTh0
    hEr0 = input.hEr0
    hFrx0 = input.hFrx0
    wall = input.wall
    gmax = input.gmax
    mcool = input.mcool
    hpcool = input.hpcool
    mclosure = input.mclosure
    hpclosure = input.hpclosure
    fclosure = input.fclosure
    hprmodel = input.hprmodel
    hsolver = input.hsolver
    rsolver = input.rsolver
    retry = input.retry
    reject = input.reject

    reset!(reject)

    ng = ngroup(opacity)

    monitor = Monitor1D()

    if !isfinite(dx) || dx <= T(0)
        abort(solver, "Coupler1DInput invalid dx", 0, T(0), T(0), dx, domain, 0, 0, T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0))
    end

    if !isfinite(tsim) || tsim <= T(0)
        abort(solver, "Coupler1DInput invalid tsim", 0, T(0), T(0), dx, domain, 0, 0, T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0))
    end

    if !isfinite(courant) || courant <= T(0)
        abort(solver, "Coupler1DInput invalid courant", 0, T(0), T(0), dx, domain, 0, 0, T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0))
    end

    if mode != "classic" && mode != "meso"
        abort(solver, "Coupler1DInput invalid mode", 0, T(0), T(0), dx, domain, 0, 0, T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0))
    end

    if !isfinite(gamma) || gamma <= T(1)
        abort(solver, "Coupler1DInput invalid gamma", 0, T(0), T(0), dx, domain, 0, 0, T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0))
    end

    if !isfinite(mu) || mu <= T(0)
        abort(solver, "Coupler1DInput invalid mu", 0, T(0), T(0), dx, domain, 0, 0, T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0))
    end

    if !isfinite(tolsim) || tolsim <= T(0)
        abort(solver, "Coupler1DInput invalid tolsim", 0, T(0), T(0), dx, domain, 0, 0, T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0))
    end

    if ng < 1
        abort(solver, "Coupler1DInput invalid ngroup", 0, T(0), T(0), dx, domain, 0, 0, T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0))
    end

    if length(gE) != ng + 1
        abort(solver, "Coupler1DInput inconsistent gE size", 0, T(0), T(0), dx, domain, 0, 0, T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0))
    end

    if nspectral < 1
        abort(solver, "Coupler1DInput invalid nspectral", 0, T(0), T(0), dx, domain, 0, 0, T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0))
    end

    if nangular < 1
        abort(solver, "Coupler1DInput invalid nangular", 0, T(0), T(0), dx, domain, 0, 0, T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0))
    end

    if !isfinite(rho0) || rho0 <= T(0)
        abort(solver, "Coupler1DInput invalid rho0", 0, T(0), T(0), dx, domain, 0, 0, vx0, T(0), T(0), Th0, rho0, T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0))
    end

    if !isfinite(vx0) || abs(vx0) >= T(c)
        abort(solver, "Coupler1DInput invalid vx0", 0, T(0), T(0), dx, domain, 0, 0, vx0, T(0), T(0), Th0, rho0, T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0))
    end

    if !isfinite(Th0) || Th0 <= T(0)
        abort(solver, "Coupler1DInput invalid Th0", 0, T(0), T(0), dx, domain, 0, 0, vx0, T(0), T(0), Th0, rho0, T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0))
    end

    if length(Er0) != ng
        abort(solver, "Coupler1DInput inconsistent Er0 size", 0, T(0), T(0), dx, domain, 0, 0, T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0))
    end

    if length(Frx0) != ng
        abort(solver, "Coupler1DInput inconsistent Frx0 size", 0, T(0), T(0), dx, domain, 0, 0, T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0))
    end

    if length(hrho0) != 2
        abort(solver, "Coupler1DInput inconsistent hrho0 size", 0, T(0), T(0), dx, domain, 0, 0, T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0))
    end

    if length(hvx0) != 2
        abort(solver, "Coupler1DInput inconsistent hvx0 size", 0, T(0), T(0), dx, domain, 0, 0, T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0))
    end

    if length(hTh0) != 2
        abort(solver, "Coupler1DInput inconsistent hTh0 size", 0, T(0), T(0), dx, domain, 0, 0, T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0))
    end

    if length(hEr0) != 2
        abort(solver, "Coupler1DInput inconsistent hEr0 size", 0, T(0), T(0), dx, domain, 0, 0, T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0))
    end

    if length(hFrx0) != 2
        abort(solver, "Coupler1DInput inconsistent hFrx0 size", 0, T(0), T(0), dx, domain, 0, 0, T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0))
    end

    if isempty(strip(wall))
        abort(solver, "Coupler1DInput invalid wall", 0, T(0), T(0), dx, domain, 0, 0, T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0))
    end

    if !isfinite(gmax) || gmax < T(1)
        abort(solver, "Coupler1DInput invalid gmax", 0, T(0), T(0), dx, domain, 0, 0, T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0))
    end

    if mcool != "none" && mcool != "fplanck" && mcool != "fpower"
        abort(solver, "Coupler1DInput invalid mcool", 0, T(0), T(0), dx, domain, 0, 0, T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0))
    end

    if length(hpcool) != 6 * ng
        abort(solver, "Coupler1DInput inconsistent hpcool size", 0, T(0), T(0), dx, domain, 0, 0, T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0))
    end

    if any(x -> !isfinite(x), hpcool)
        abort(solver, "Coupler1DInput invalid hpcool", 0, T(0), T(0), dx, domain, 0, 0, T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0))
    end

    if isempty(strip(mclosure))
        abort(solver, "Coupler1DInput invalid mclosure", 0, T(0), T(0), dx, domain, 0, 0, T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0))
    end

    if any(x -> !isfinite(x), hpclosure)
        abort(solver, "Coupler1DInput invalid hpclosure", 0, T(0), T(0), dx, domain, 0, 0, T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0))
    end

    if isempty(strip(fclosure))
        abort(solver, "Coupler1DInput invalid fclosure", 0, T(0), T(0), dx, domain, 0, 0, T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0))
    end

    if length(hprmodel) != 2
        abort(solver, "Coupler1DInput inconsistent hprmodel size", 0, T(0), T(0), dx, domain, 0, 0, T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0))
    end

    if !isfinite(hprmodel[1]) || hprmodel[1] <= T(0)
        abort(solver, "Coupler1DInput invalid hprmodel", 0, T(0), T(0), dx, domain, 0, 0, T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0))
    end

    if !isfinite(hprmodel[2]) || hprmodel[2] <= hprmodel[1]
        abort(solver, "Coupler1DInput invalid hprmodel", 0, T(0), T(0), dx, domain, 0, 0, T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0))
    end

    if isempty(strip(hsolver))
        abort(solver, "Coupler1DInput invalid hsolver", 0, T(0), T(0), dx, domain, 0, 0, T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0))
    end

    if isempty(strip(rsolver))
        abort(solver, "Coupler1DInput invalid rsolver", 0, T(0), T(0), dx, domain, 0, 0, T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0))
    end

    if retry < 0
        abort(solver, "Coupler1DInput invalid retry", 0, T(0), T(0), dx, domain, 0, 0, T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0))
    end

    rho = zeros(T, nxlocal, 1, 1)
    mom = zeros(T, nxlocal, 1, 1, 1)
    Eh = zeros(T, nxlocal, 1, 1)
    Th = zeros(T, nxlocal, 1, 1)

    Er = zeros(T, nxlocal, 1, 1, ng)
    Fr = zeros(T, nxlocal, 1, 1, ng, 1)
    Pr = zeros(T, nxlocal, 1, 1, ng, 1, 1)
    regime = Array{Union{Nothing,Symbol},4}(nothing, nxlocal, 1, 1, ng)

    # Initial hydro state
    for jxlocal in 1:nxlocal
        jx = domain.jxbegin + jxlocal - 1
        x = (T(jx) - T(0.5)) * dx

        rhoj = x < hrho0[1] * nx * dx ? rho0 * (T(1) + hrho0[2]) : rho0 * (T(1) - hrho0[2])
        vxj = x < hvx0[1] * nx * dx ? vx0 * (T(1) + hvx0[2]) : vx0 * (T(1) - hvx0[2])
        Thj = x < hTh0[1] * nx * dx ? Th0 * (T(1) + hTh0[2]) : Th0 * (T(1) - hTh0[2])

        if !isfinite(rhoj) || rhoj <= T(0)
            abort(solver, "Coupler1D invalid initial density", 0, T(0), T(0), dx, domain, jx, 0, vxj, T(0), T(0), Thj, rhoj, T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0))
        end

        if !isfinite(vxj) || abs(vxj) >= T(c)
            abort(solver, "Coupler1D invalid initial velocity", 0, T(0), T(0), dx, domain, jx, 0, vxj, T(0), T(0), Thj, rhoj, T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0))
        end

        if !isfinite(Thj) || Thj <= T(0)
            abort(solver, "Coupler1D invalid initial Th", 0, T(0), T(0), dx, domain, jx, 0, vxj, T(0), T(0), Thj, rhoj, T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0))
        end

        pj = rhoj * T(kB) * Thj / (mu * T(mp))
        Ehj = pj / (gamma - T(1)) + T(0.5) * rhoj * vxj^2

        if !isfinite(pj) || pj <= T(0)
            abort(solver, "Coupler1D invalid initial pressure", 0, T(0), T(0), dx, domain, jx, 0, vxj, T(0), pj, Thj, rhoj, rhoj * vxj, Ehj, T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0))
        end

        if !isfinite(Ehj) || Ehj <= T(0)
            abort(solver, "Coupler1D invalid initial energy", 0, T(0), T(0), dx, domain, jx, 0, vxj, T(0), pj, Thj, rhoj, rhoj * vxj, Ehj, T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0))
        end

        rho[jxlocal, 1, 1] = rhoj
        mom[jxlocal, 1, 1, 1] = rhoj * vxj
        Eh[jxlocal, 1, 1] = Ehj
        Th[jxlocal, 1, 1] = Thj
    end

    # Initial radiation state
    for l in 1:ng
        for jxlocal in 1:nxlocal
            jx = domain.jxbegin + jxlocal - 1
            x = (T(jx) - T(0.5)) * dx

            rhoj = rho[jxlocal, 1, 1]
            momj = mom[jxlocal, 1, 1, 1]
            Ehj = Eh[jxlocal, 1, 1]
            Thj = Th[jxlocal, 1, 1]

            Er0l = x < hEr0[1] * nx * dx ? Er0[l] * (T(1) + hEr0[2]) : Er0[l] * (T(1) - hEr0[2])
            Frx0l = x < hFrx0[1] * nx * dx ? Frx0[l] * (T(1) + hFrx0[2]) : Frx0[l] * (T(1) - hFrx0[2])

            kRjl = kappaR(opacity, Thj, rhoj, l)

            if !isfinite(kRjl) || kRjl <= T(0)
                abort(solver, "Coupler1D invalid initial kappaR", 0, T(0), T(0), dx, domain, jx, l, T(0), T(0), T(0), Thj, rhoj, momj, Ehj, T(0), T(0), T(0), Er0l, Frx0l, T(0), T(0), T(0), T(0))
            end

            ceffjl = mode == "classic" ? T(c) : min(T(c), T(c) / (rhoj * kRjl * dx))

            if !isfinite(ceffjl) || ceffjl <= T(0)
                abort(solver, "Coupler1D invalid initial ceff", 0, T(0), T(0), dx, domain, jx, l, T(0), T(0), T(0), Thj, rhoj, momj, Ehj, T(0), T(0), ceffjl, Er0l, Frx0l, T(0), T(0), T(0), T(0))
            end

            if !isfinite(Er0l) || Er0l <= T(0)
                abort(solver, "Coupler1D invalid initial Er", 0, T(0), T(0), dx, domain, jx, l, T(0), T(0), T(0), Thj, rhoj, momj, Ehj, T(0), T(0), ceffjl, Er0l, Frx0l, T(0), T(0), T(0), T(0))
            end

            if !isfinite(Frx0l)
                abort(solver, "Coupler1D invalid initial Fr", 0, T(0), T(0), dx, domain, jx, l, T(0), T(0), T(0), Thj, rhoj, momj, Ehj, T(0), T(0), ceffjl, Er0l, Frx0l, T(0), T(0), T(0), T(0))
            end

            if abs(Frx0l) > ceffjl * Er0l
                abort(solver, "Coupler1DInput invalid initial reduced flux", 0, T(0), T(0), dx, domain, jx, l, T(0), T(0), T(0), Thj, rhoj, momj, Ehj, T(0), T(0), ceffjl, Er0l, Frx0l, T(0), T(0), T(0), T(0))
            end

            Er[jxlocal, 1, 1, l] = Er0l
            Fr[jxlocal, 1, 1, l, 1] = Frx0l
        end
    end

    t = T(0)
    it = 0
    t0 = clock_ns()

    _fprimitive = let

        _solver = Source1DSolver{Source1DInput{T},Source1DOutput{T}}()

        _input = Source1DInput(
            T(1),
            nx,
            dx,
            domain,
            mode,
            gamma,
            mu,
            tolsim,
            gE,
            opacity,
            nspectral,
            rho,
            mom,
            Eh,
            Er,
            Fr,
            Pr,
            gmax,
            mcool,
            hpcool,
            regime,
            reject
        )

        (rhoj, momj, Ehj) -> fprimitive(_solver, _input, rhoj, momj, Ehj)
    end

    while t < tsim
        it += 1

        dth = typemax(T)
        dtr = typemax(T)

        # Explicit time-step bounds
        for jxlocal in 1:nxlocal
            jx = domain.jxbegin + jxlocal - 1
            rhoj = rho[jxlocal, 1, 1]
            momj = mom[jxlocal, 1, 1, 1]
            Ehj = Eh[jxlocal, 1, 1]

            vx, cs, p, Thj, _, _ = _fprimitive(rhoj, momj, Ehj)

            # Hydro time-step bound
            lambda = abs(vx) + cs
            dth = min(dth, courant * dx / lambda)

            # Radiation time-step bound
            for l in 1:ng
                kRjl = kappaR(opacity, Thj, rhoj, l)

                if !isfinite(kRjl) || kRjl <= T(0)
                    abort(solver, "Coupler1D invalid local kappaR", it, t, T(0), dx, domain, jx, l, vx, cs, p, Thj, rhoj, momj, Ehj, T(0), T(0), T(0), Er[jxlocal, 1, 1, l], Fr[jxlocal, 1, 1, l, 1], T(0), T(0), T(0), T(0))
                end

                ceffjl = mode == "classic" ? T(c) : min(T(c), T(c) / (rhoj * kRjl * dx))

                if !isfinite(ceffjl) || ceffjl <= T(0)
                    abort(solver, "Coupler1D invalid local ceff", it, t, T(0), dx, domain, jx, l, vx, cs, p, Thj, rhoj, momj, Ehj, T(0), T(0), ceffjl, Er[jxlocal, 1, 1, l], Fr[jxlocal, 1, 1, l, 1], T(0), T(0), T(0), T(0))
                end

                dtr = min(dtr, courant * dx / ceffjl)
            end
        end

        # Local time step
        dtlocal = min(dth, dtr)

        # Global time step
        dt = globalmin(domain, dtlocal)

        if !isfinite(dt) || dt <= T(0)
            reject!(reject, "Coupler1D invalid | dt = $(dt)")
            abort(solver, "Coupler1D invalid dt", it, t, dt, dx, domain, 0, 0, T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0))
        end

        # Adjust final time step
        if t + dt > tsim
            dt = tsim - t
        end

        accepted = false

        for itry in 0:retry
            _reject = reject
            reset!(_reject)
            _dt = dt / T(2)^itry

            if !isfinite(_dt) || _dt <= T(0)
                reject!(_reject, "Coupler1D invalid | _dt = $(_dt)")
                continue
            end

            _rho = copy(rho)
            _mom = copy(mom)
            _Eh = copy(Eh)
            _Th = copy(Th)

            _Er = copy(Er)
            _Fr = copy(Fr)
            _Pr = copy(Pr)

            # Hydro transport
            _input = Hydro1DInput(
                _dt,
                nx,
                dx,
                domain,
                gamma,
                tolsim,
                _rho,
                _mom,
                _Eh,
                wall,
                hsolver,
                retry,
                _reject
            )

            _solver = Hydro1DSolver{Hydro1DInput{T},Hydro1DOutput{T}}()
            output = kernel(_solver, _input)

            if globalany(domain, _reject.flag)
                continue
            end

            # Intermediate Hydro state
            _rho = output.rho
            _mom = output.mom
            _Eh = output.Eh

            if any(x -> !isfinite(x), _rho) ||
               any(x -> !isfinite(x), _mom) ||
               any(x -> !isfinite(x), _Eh)
                reject!(_reject, "Coupler1D invalid | _rho = $(_rho) _mom = $(_mom) _Eh = $(_Eh)")
            end

            if globalany(domain, _reject.flag)
                continue
            end

            # Hydro temperature reconstruction
            for jxlocal in 1:nxlocal
                jx = domain.jxbegin + jxlocal - 1
                rhoj = _rho[jxlocal, 1, 1]
                momj = _mom[jxlocal, 1, 1, 1]
                Ehj = _Eh[jxlocal, 1, 1]

                if !isfinite(rhoj) || rhoj <= T(0)
                    reject!(_reject, "Coupler1D invalid | rhoj = $(rhoj)")
                    break
                end

                if !isfinite(momj) || !isfinite(Ehj)
                    reject!(_reject, "Coupler1D invalid | momj = $(momj) Ehj = $(Ehj)")
                    break
                end

                vxj = momj / rhoj
                Ekj = T(0.5) * rhoj * vxj^2
                Ethj = Ehj - Ekj

                Ethj = corrector(_solver, Ethj, tolsim, Ehj, tolsim, true)

                if !isfinite(Ethj) || Ethj < tolsim || Ethj > Ehj
                    reject!(_reject, "Coupler1D invalid | Ethj = $(Ethj)")
                    break
                end

                Ehj = Ekj + Ethj
                _Eh[jxlocal, 1, 1] = Ehj

                vx, cs, p, Thj, _, _ = _fprimitive(rhoj, momj, Ehj)

                if !isfinite(vx) || abs(vx) >= T(c)
                    reject!(_reject, "Coupler1D invalid | vx = $(vx)")
                    break
                end

                if !isfinite(cs) || cs <= T(0)
                    reject!(_reject, "Coupler1D invalid | cs = $(cs)")
                    break
                end

                if !isfinite(p) || p <= T(0)
                    reject!(_reject, "Coupler1D invalid | p = $(p)")
                    break
                end

                if !isfinite(Thj) || Thj <= T(0)
                    reject!(_reject, "Coupler1D invalid | Thj = $(Thj)")
                    break
                end

                _Th[jxlocal, 1, 1] = Thj
            end

            if globalany(domain, _reject.flag)
                continue
            end

            # Radiation transport
            _input = Radiation1DInput(
                _dt,
                nx,
                dx,
                domain,
                mode,
                tolsim,
                gE,
                opacity,
                nspectral,
                nangular,
                _rho,
                _Th,
                _Er,
                _Fr,
                wall,
                mclosure,
                hpclosure,
                fclosure,
                hprmodel,
                rsolver,
                retry,
                _reject
            )

            _solver = Radiation1DSolver{Radiation1DInput{T},Radiation1DOutput{T}}()
            output = kernel(_solver, _input)

            if globalany(domain, _reject.flag)
                continue
            end

            # Intermediate Radiation state
            _Er = output.Er
            _Fr = output.Fr
            _Pr = output.Pr
            _regime = output.regime

            if any(x -> !isfinite(x), _Er) ||
               any(x -> !isfinite(x), _Fr) ||
               any(x -> !isfinite(x), _Pr)
                reject!(_reject, "Coupler1D invalid | _Er = $(_Er) _Fr = $(_Fr) _Pr = $(_Pr)")
            end

            if any(x -> x !== :thick && x !== :m1 && x !== :thin, _regime)
                reject!(_reject, "Coupler1D invalid | _regime = $(_regime)")
            end

            if globalany(domain, _reject.flag)
                continue
            end

            # Source coupling
            _input = Source1DInput(
                _dt,
                nx,
                dx,
                domain,
                mode,
                gamma,
                mu,
                tolsim,
                gE,
                opacity,
                nspectral,
                _rho,
                _mom,
                _Eh,
                _Er,
                _Fr,
                _Pr,
                gmax,
                mcool,
                hpcool,
                _regime,
                _reject
            )

            _solver = Source1DSolver{Source1DInput{T},Source1DOutput{T}}()
            output = kernel(_solver, _input)

            if globalany(domain, _reject.flag)
                continue
            end

            cS0 = output.cS0
            Sx = output.Sx

            if any(x -> !isfinite(x), cS0) ||
               any(x -> !isfinite(x), Sx)
                reject!(_reject, "Coupler1D invalid | cS0 = $(cS0) Sx = $(Sx)")
            end

            if globalany(domain, _reject.flag)
                continue
            end

            for jxlocal in 1:nxlocal
                jx = domain.jxbegin + jxlocal - 1
                x = (T(jx) - T(0.5)) * dx

                rhoj = _rho[jxlocal, 1, 1]
                momj = _mom[jxlocal, 1, 1, 1]
                Ehj = _Eh[jxlocal, 1, 1]

                if !isfinite(rhoj) || rhoj <= T(0)
                    reject!(_reject, "Coupler1D invalid | rhoj = $(rhoj)")
                    break
                end

                if !isfinite(momj) || !isfinite(Ehj)
                    reject!(_reject, "Coupler1D invalid | momj = $(momj) Ehj = $(Ehj)")
                    break
                end

                vxj = momj / rhoj
                Ekj = T(0.5) * rhoj * vxj^2
                Ethj = Ehj - Ekj

                Ethj = corrector(_solver, Ethj, tolsim, Ehj, tolsim, true)

                if !isfinite(Ethj) || Ethj < tolsim || Ethj > Ehj
                    reject!(_reject, "Coupler1D invalid | Ethj = $(Ethj)")
                    break
                end

                Ehj = Ekj + Ethj
                _Eh[jxlocal, 1, 1] = Ehj

                vx, cs, p, Thj, _, _ = _fprimitive(rhoj, momj, Ehj)

                if !isfinite(vx) || abs(vx) >= T(c)
                    reject!(_reject, "Coupler1D invalid | vx = $(vx)")
                    break
                end

                if !isfinite(cs) || cs <= T(0)
                    reject!(_reject, "Coupler1D invalid | cs = $(cs)")
                    break
                end

                if !isfinite(p) || p <= T(0)
                    reject!(_reject, "Coupler1D invalid | p = $(p)")
                    break
                end

                if !isfinite(Thj) || Thj <= T(0)
                    reject!(_reject, "Coupler1D invalid | Thj = $(Thj)")
                    break
                end

                cS0j = cS0[jxlocal, 1, 1, :]
                Sxj = Sx[jxlocal, 1, 1, :, 1]

                dEh = T(0)
                dmom = T(0)
                Q0 = T(0)
                Qx = T(0)

                for l in 1:ng
                    kRjl = kappaR(opacity, Thj, rhoj, l)

                    if !isfinite(kRjl) || kRjl <= T(0)
                        reject!(_reject, "Coupler1D invalid | kRjl = $(kRjl)")
                        break
                    end

                    mfreejl = T(1) / (rhoj * kRjl)
                    ceffjl = mode == "classic" ? T(c) : min(T(c), T(c) / (rhoj * kRjl * dx))

                    if !isfinite(ceffjl) || ceffjl <= T(0)
                        reject!(_reject, "Coupler1D invalid | ceffjl = $(ceffjl)")
                        break
                    end

                    if mode == "meso"
                        cS0jl, Sxjl = tfmeso(_solver, _input, cS0j, Sxj, l, ceffjl)
                    else
                        cS0jl, Sxjl = cS0j[l], Sxj[l]
                    end

                    if !isfinite(cS0jl) || !isfinite(Sxjl)
                        reject!(_reject, "Coupler1D invalid | cS0jl = $(cS0jl) Sxjl = $(Sxjl)")
                        break
                    end

                    Erjl = _Er[jxlocal, 1, 1, l]
                    Frjl = _Fr[jxlocal, 1, 1, l, 1]
                    Prjl = _Pr[jxlocal, 1, 1, l, 1, 1]

                    if !isfinite(Erjl) || Erjl <= T(0)
                        reject!(_reject, "Coupler1D invalid | Erjl = $(Erjl)")
                        break
                    end

                    if !isfinite(Frjl)
                        reject!(_reject, "Coupler1D invalid | Frjl = $(Frjl)")
                        break
                    end

                    if !isfinite(Prjl) || Prjl < T(0)
                        reject!(_reject, "Coupler1D invalid | Prjl = $(Prjl)")
                        break
                    end

                    Ernew = Erjl - _dt * cS0jl
                    Frnew = Frjl - _dt * ceffjl^2 * Sxjl

                    if !isfinite(Ernew) || Ernew <= T(0)
                        reject!(_reject, "Coupler1D invalid | Ernew = $(Ernew)")
                        break
                    end

                    Frnew = corrector(_solver, Frnew, -(T(1) - eps(T)) * ceffjl * Ernew, (T(1) - eps(T)) * ceffjl * Ernew, tolsim, true)

                    if abs(Frnew) > (T(1) - eps(T)) * ceffjl * Ernew
                        reject!(_reject, "Coupler1D invalid | Frnew = $(Frnew)")
                        break
                    end

                    Sxjl = (Frjl - Frnew) / (_dt * ceffjl^2)

                    if !isfinite(Sxjl)
                        reject!(_reject, "Coupler1D invalid | Sxjl = $(Sxjl)")
                        break
                    end

                    _Er[jxlocal, 1, 1, l] = Ernew
                    _Fr[jxlocal, 1, 1, l, 1] = Frnew

                    dEh += _dt * cS0jl
                    dmom += _dt * Sxjl
                    Q0 += cS0jl
                    Qx += Sxjl

                    Trjl = (Ernew / T(aR))^(T(1) / T(4))

                    if !isfinite(Trjl) || Trjl <= T(0)
                        reject!(_reject, "Coupler1D invalid | Trjl = $(Trjl)")
                        break
                    end

                    monitor!(
                        monitor,
                        it,
                        t,
                        _dt,
                        jxlocal,
                        nxlocal,
                        x,
                        domain.rank,
                        l,
                        vx,
                        cs,
                        p,
                        Thj,
                        rhoj,
                        momj,
                        Ehj,
                        Q0,
                        Qx,
                        mfreejl,
                        ceffjl,
                        Ernew,
                        Frnew,
                        Prjl,
                        Trjl,
                        cS0jl,
                        Sxjl,
                        _regime[jxlocal, 1, 1, l]
                    )
                end

                if _reject.flag
                    break
                end

                _mom[jxlocal, 1, 1, 1] = momj + dmom
                _Eh[jxlocal, 1, 1] = Ehj + dEh

                if !isfinite(_mom[jxlocal, 1, 1, 1]) || !isfinite(_Eh[jxlocal, 1, 1])
                    reject!(_reject, "Coupler1D invalid | _mom[jxlocal, 1, 1, 1] = $(_mom[jxlocal, 1, 1, 1]) _Eh[jxlocal, 1, 1] = $(_Eh[jxlocal, 1, 1])")
                    break
                end

                vxnew = _mom[jxlocal, 1, 1, 1] / _rho[jxlocal, 1, 1]
                Eknew = T(0.5) * _rho[jxlocal, 1, 1] * vxnew^2
                Ethnew = _Eh[jxlocal, 1, 1] - Eknew

                Ethnew = corrector(_solver, Ethnew, tolsim, _Eh[jxlocal, 1, 1], tolsim, true)

                if !isfinite(Ethnew) || Ethnew < tolsim || Ethnew > _Eh[jxlocal, 1, 1]
                    reject!(_reject, "Coupler1D invalid | Ethnew = $(Ethnew)")
                    break
                end

                _Eh[jxlocal, 1, 1] = Eknew + Ethnew
            end

            if globalany(domain, _reject.flag)
                continue
            end

            if any(x -> !isfinite(x), _rho) ||
               any(x -> !isfinite(x), _mom) ||
               any(x -> !isfinite(x), _Eh)
                reject!(_reject, "Coupler1D invalid | _rho = $(_rho) _mom = $(_mom) _Eh = $(_Eh)")
            end

            if globalany(domain, _reject.flag)
                continue
            end

            if any(x -> !isfinite(x), _Er) ||
               any(x -> !isfinite(x), _Fr) ||
               any(x -> !isfinite(x), _Pr)
                reject!(_reject, "Coupler1D invalid | _Er = $(_Er) _Fr = $(_Fr) _Pr = $(_Pr)")
            end

            if globalany(domain, _reject.flag)
                continue
            end

            # Hydro temperature reconstruction
            for jxlocal in 1:nxlocal
                jx = domain.jxbegin + jxlocal - 1
                rhoj = _rho[jxlocal, 1, 1]
                momj = _mom[jxlocal, 1, 1, 1]
                Ehj = _Eh[jxlocal, 1, 1]

                if !isfinite(rhoj) || rhoj <= T(0)
                    reject!(_reject, "Coupler1D invalid | rhoj = $(rhoj)")
                    break
                end

                if !isfinite(momj) || !isfinite(Ehj)
                    reject!(_reject, "Coupler1D invalid | momj = $(momj) Ehj = $(Ehj)")
                    break
                end

                vxj = momj / rhoj
                Ekj = T(0.5) * rhoj * vxj^2
                Ethj = Ehj - Ekj

                Ethj = corrector(_solver, Ethj, tolsim, Ehj, tolsim, true)

                if !isfinite(Ethj) || Ethj < tolsim || Ethj > Ehj
                    reject!(_reject, "Coupler1D invalid | Ethj = $(Ethj)")
                    break
                end

                Ehj = Ekj + Ethj
                _Eh[jxlocal, 1, 1] = Ehj

                _, _, _, Thj, _, _ = _fprimitive(rhoj, momj, Ehj)

                if !isfinite(Thj) || Thj <= T(0)
                    reject!(_reject, "Coupler1D invalid | Thj = $(Thj)")
                    break
                end

                _Th[jxlocal, 1, 1] = Thj
            end

            if globalany(domain, _reject.flag)
                continue
            end

            # Final Hydro state
            rho = _rho
            mom = _mom
            Eh = _Eh
            Th = _Th

            # Final Radiation state
            Er = _Er
            Fr = _Fr
            Pr = _Pr

            t += _dt
            accepted = true

            break
        end

        if !accepted
            reject!(reject, "Coupler1D invalid | accepted = $(accepted)")
            abort(solver, "Coupler1D retry failed", it, t, dt, dx, domain, 0, 0, T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0))
        end

        if !isfinite(t)
            reject!(reject, "Coupler1D invalid | t = $(t)")
            abort(solver, "Coupler1D invalid time", it, t, dt, dx, domain, 0, 0, T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0))
        end

        if domain.rank == 0
            elapsed = (clock_ns() - t0) / 1.0e9
            ratio = Float64(t / tsim)
            progress("HDND", ratio, elapsed, width=30, height=1)
        end
    end

    return Coupler1DOutput{T}()
end
