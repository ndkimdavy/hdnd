include("../../interface/isolver.jl")
include("../../interface/iopacity.jl")
include("../constant.jl")
include("../reject.jl")

using Base.Threads

#===================================================================
                            Source 1D
          (rho, mom, Eh, Er, Fr, Pr)[*] -> (S0, Sx)[lab]
===================================================================#

struct Source1DInput{T} <: ISolverInput
    dt::T
    nx::Int
    dx::T
    domain::Domain1D
    mode::String
    gamma::T
    mu::T
    tolsim::T
    gE::Array{T,1}
    opacity::IOpacity{T}
    nspectral::Int
    rho::Array{T,3}
    mom::Array{T,4}
    Eh::Array{T,3}
    Er::Array{T,4}
    Fr::Array{T,5}
    Pr::Array{T,6}
    gmax::T
    mcool::String
    hpcool::Array{T,1}
    regime::Array{Union{Nothing,Symbol},4}
    reject::Reject
end

struct Source1DOutput{T} <: ISolverOutput
    cS0::Array{T,4}
    Sx::Array{T,5}
end

struct Source1DSolver{Tin<:ISolverInput,Tout<:ISolverOutput} <: ISolver{Tin,Tout} end

function abort(solver::Source1DSolver{Source1DInput{T},Source1DOutput{T}},
    msg::String,
    dt::T,
    nx::Int, dx::T, domain::Domain1D,
    jx::Int,
    l::Int,
    vx::T, cs::T, p::T, Thj::T, rhoj::T, momj::T, Ehj::T, Q0::T, Qx::T,
    ceffjl::T, Erjl::T, Frjl::T, Prjl::T, Trjl::T, cS0::T, Sx::T) where {T<:AbstractFloat}

    log("ERROR: Source1D")
    log("solver  = $(typeof(solver))")
    log("dt      = $(dt)")
    log("nx      = $(nx)")
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

function fprimitive(solver::Source1DSolver{Source1DInput{T},Source1DOutput{T}}, input::Source1DInput{T}, rhoj::T, momj::T, Ehj::T)::Tuple{T,T,T,T,T,T} where {T<:AbstractFloat}

    gamma = input.gamma
    mu = input.mu

    vx = momj / rhoj
    Ekj = T(0.5) * rhoj * vx^2
    Ethj = Ehj - Ekj
    p = (gamma - T(1)) * Ethj
    Thj = (mu * T(mp) / T(kB)) * (p / rhoj)
    cs = sqrt(gamma * p / rhoj)

    return vx, cs, p, Thj, Ekj, Ethj
end

function fBl(solver::Source1DSolver{Source1DInput{T},Source1DOutput{T}}, input::Source1DInput{T}, Thj::T, l::Int)::T where {T<:AbstractFloat}

    gE = input.gE
    nspectral = input.nspectral

    E1 = gE[l] * T(keV_to_J)
    E2 = gE[l+1] * T(keV_to_J)

    dE = (E2 - E1) / T(nspectral)
    Bjl = T(0)

    for m in 1:nspectral
        Ephm = E1 + (T(m) - T(0.5)) * dE
        nu = Ephm / T(h)
        dnu = dE / T(h)
        argexp = Ephm / (T(kB) * Thj)

        if isinf(argexp)
            Bjlm = T(0)
        else
            den = exp(argexp) - T(1)
            Bjlm = isinf(den) ? T(0) : (T(2) * T(h) * nu^3 / T(c)^2) / den
        end

        Bjl += Bjlm * dnu
    end

    return Bjl
end

function fplanck(solver::Source1DSolver{Source1DInput{T},Source1DOutput{T}}, input::Source1DInput{T}, rhoj::T, Thj::T, kPjl::T, l::Int)::T where {T<:AbstractFloat}

    Bjl = fBl(solver, input, Thj, l)
    Eeq = T(4) * T(pi) * Bjl / T(c)
    Lambda = rhoj * kPjl * Eeq

    return Lambda
end

function fpower(solver::Source1DSolver{Source1DInput{T},Source1DOutput{T}}, input::Source1DInput{T}, rhoj::T, p::T, x::T, l::Int)::T where {T<:AbstractFloat}

    hpcool = input.hpcool

    i = 6 * (l - 1)

    C0 = hpcool[i+1]
    a = hpcool[i+2]
    b = hpcool[i+3]
    c = hpcool[i+4]
    d = hpcool[i+5]
    e = hpcool[i+6]

    y = T(1)
    z = T(1)

    Lambda = C0 * rhoj^a * p^b * x^c * y^d * z^e

    return Lambda
end

function tfmeso(solver::Source1DSolver{Source1DInput{T},Source1DOutput{T}}, input::Source1DInput{T}, cS0j::Array{T,1}, Sxj::Array{T,1}, l::Int, ceffjl::T)::Tuple{T,T} where {T<:AbstractFloat}

    gmax = input.gmax

    scale = min(T(c) / ceffjl, gmax)

    return scale * cS0j[l], scale * Sxj[l]
end

function fexplicit(solver::Source1DSolver{Source1DInput{T},Source1DOutput{T}}, input::Source1DInput{T}, jxlocal::Int, reject::Reject)::Tuple{Array{T,1},Array{T,1}} where {T<:AbstractFloat}

    dx = input.dx
    domain = input.domain
    opacity = input.opacity
    rho = input.rho
    mom = input.mom
    Eh = input.Eh
    Er = input.Er
    Fr = input.Fr
    Pr = input.Pr
    mcool = input.mcool
    regime = input.regime

    ng = ngroup(opacity)

    cS0 = zeros(T, ng)
    Sx = zeros(T, ng)

    jx = domain.jxbegin + jxlocal - 1
    x = (T(jx) - T(0.5)) * dx

    rhoj = rho[jxlocal, 1, 1]
    momj = mom[jxlocal, 1, 1, 1]
    Ehj = Eh[jxlocal, 1, 1]

    vx, _, p, Thj, _, _ = fprimitive(solver, input, rhoj, momj, Ehj)

    b = vx / T(c)
    g = T(1) / sqrt(T(1) - b^2)

    # Lorentz matrix: laboratory -> comobile
    L = T[
        g -b*g
        -b*g g
    ]

    # Lorentz matrix: comobile -> laboratory
    iL = T[
        g b*g
        b*g g
    ]

    for l in 1:ng
        regimejl = regime[jxlocal, 1, 1, l]

        # Group source in comobile frame
        if regimejl === :thin && mcool == "fpower"
            Lambda = fpower(solver, input, rhoj, p, x, l)

            if !isfinite(Lambda) || Lambda < T(0)
                reject!(reject, "Source1D invalid | Lambda = $(Lambda)")
                return cS0, Sx
            end

            S0com = -Lambda
            Sxcom = T(0)

        elseif regimejl === :thin && mcool == "fplanck"
            kPjl = kappaP(opacity, Thj, rhoj, l)

            if !isfinite(kPjl) || kPjl <= T(0)
                reject!(reject, "Source1D invalid | kPjl = $(kPjl)")
                return cS0, Sx
            end

            Lambda = fplanck(solver, input, rhoj, Thj, kPjl, l)

            if !isfinite(Lambda) || Lambda < T(0)
                reject!(reject, "Source1D invalid | Lambda = $(Lambda)")
                return cS0, Sx
            end

            S0com = -Lambda
            Sxcom = T(0)

        elseif regimejl === :thick || regimejl === :m1 || regimejl === :thin
            Erjl = Er[jxlocal, 1, 1, l]
            Frjl = Fr[jxlocal, 1, 1, l, 1]
            Prjl = Pr[jxlocal, 1, 1, l, 1, 1]

            kRjl = kappaR(opacity, Thj, rhoj, l)
            kPjl = kappaP(opacity, Thj, rhoj, l)

            if !isfinite(kRjl) || kRjl <= T(0)
                reject!(reject, "Source1D invalid | kRjl = $(kRjl)")
                return cS0, Sx
            end

            if !isfinite(kPjl) || kPjl <= T(0)
                reject!(reject, "Source1D invalid | kPjl = $(kPjl)")
                return cS0, Sx
            end

            # Radiation matrix in laboratory frame
            R = T[
                Erjl Frjl/T(c)
                Frjl/T(c) Prjl
            ]

            # Laboratory -> comobile
            Rc = L * R * transpose(L)

            Erc = Rc[1, 1]
            Frc = T(c) * Rc[1, 2]
            Prc = Rc[2, 2]

            if !isfinite(Erc) || Erc <= T(0)
                reject!(reject, "Source1D invalid | Erc = $(Erc)")
                return cS0, Sx
            end

            if !isfinite(Frc)
                reject!(reject, "Source1D invalid | Frc = $(Frc)")
                return cS0, Sx
            end

            if !isfinite(Prc) || Prc < T(0)
                reject!(reject, "Source1D invalid | Prc = $(Prc)")
                return cS0, Sx
            end

            Bjl = fBl(solver, input, Thj, l)

            if !isfinite(Bjl) || Bjl < T(0)
                reject!(reject, "Source1D invalid | Bjl = $(Bjl)")
                return cS0, Sx
            end

            # Group equilibrium energy
            Eeq = T(4) * T(pi) * Bjl / T(c)

            if !isfinite(Eeq) || Eeq < T(0)
                reject!(reject, "Source1D invalid | Eeq = $(Eeq)")
                return cS0, Sx
            end

            S0com = rhoj * kPjl * (Erc - Eeq)
            Sxcom = rhoj * kRjl * Frc / T(c)

        else
            reject!(reject, "Source1D invalid | regimejl = $(regimejl)")
            return cS0, Sx
        end

        if !isfinite(S0com) || !isfinite(Sxcom)
            reject!(reject, "Source1D invalid | S0com = $(S0com) Sxcom = $(Sxcom)")
            return cS0, Sx
        end

        Scom = T[
            S0com
            Sxcom
        ]

        # Comobile -> laboratory
        Slab = iL * Scom

        if !all(isfinite, Slab)
            reject!(reject, "Source1D invalid | Slab = $(Slab)")
            return cS0, Sx
        end

        # Group source in laboratory frame
        cS0[l] = T(c) * Slab[1]
        Sx[l] = Slab[2]
    end

    return cS0, Sx
end

function kernel(solver::Source1DSolver{Source1DInput{T},Source1DOutput{T}}, input::Source1DInput{T})::Source1DOutput{T} where {T<:AbstractFloat}

    domain = input.domain
    nxlocal = domain.nxlocal
    opacity = input.opacity
    reject = input.reject

    reset!(reject)

    ng = ngroup(opacity)

    cS0 = zeros(T, nxlocal, 1, 1, ng)
    Sx = zeros(T, nxlocal, 1, 1, ng, 1)

    nthd = Threads.nthreads()
    blockJx = [(div((ithd-1)*nxlocal, nthd)+1):div(ithd*nxlocal, nthd) for ithd in 1:nthd]
    rejectJx = [Reject(false, "") for ithd in 1:nthd]

    Threads.@threads :static for ithd in 1:nthd
        _reject = rejectJx[ithd]

        for jxlocal in blockJx[ithd]
            cS0j, Sxj = fexplicit(solver, input, jxlocal, _reject)

            if _reject.flag
                break
            end

            if !all(isfinite, cS0j) || !all(isfinite, Sxj)
                reject!(_reject, "Source1D invalid | cS0j = $(cS0j) Sxj = $(Sxj)")
                break
            end

            cS0[jxlocal, 1, 1, :] .= cS0j
            Sx[jxlocal, 1, 1, :, 1] .= Sxj
        end
    end

    for _reject in rejectJx
        if _reject.flag
            reject!(reject, _reject.msg)
            return Source1DOutput(cS0, Sx)
        end
    end

    return Source1DOutput(cS0, Sx)
end
