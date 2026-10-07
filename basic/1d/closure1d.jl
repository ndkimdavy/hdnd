include("../../interface/iclosure.jl")
include("../constant.jl")
include("../mem.jl")
include("../reject.jl")
include("../randomsearch.jl")

using JSON
using ONNXRunTime

#===================================================================
                             Closure 1D
                        (Er, Fr, ceff) -> fedd
===================================================================#

struct Closure1DInput{T} <: IClosureInput
    l::Int
    tolsim::T
    gE::Array{T,1}
    ceffjl::T
    nspectral::Int
    nangular::Int
    Erjl::T
    Frjl::T
    mclosure::String
    hpclosure::Array{T,1}
    fclosure::String
    retry::Int
    reject::Reject
end

struct Closure1DOutput{T} <: IClosureOutput
    fedd::T
end

struct Closure1D{Tin<:IClosureInput,Tout<:IClosureOutput} <: IClosure{Tin,Tout} end

function abort(closure::Closure1D{Closure1DInput{T},Closure1DOutput{T}},
    msg::String,
    l::Int,
    Erjl::T, Frjl::T, ceffjl::T,
    fred::T, fedd::T) where {T<:AbstractFloat}

    log("ERROR: Closure1D")
    log("closure = $(typeof(closure))")
    log("l       = $(l)")
    log("Erjl    = $(Erjl)")
    log("Frjl    = $(Frjl)")
    log("ceffjl  = $(ceffjl)")
    log("fred    = $(fred)")
    log("fedd    = $(fedd)")

    error("$(msg)")
end

function ftf(closure::Closure1D{Closure1DInput{T},Closure1DOutput{T}}, input::Closure1DInput{T}, Alpha::T, Beta::T)::Tuple{T,T} where {T<:AbstractFloat}

    alpha = Alpha
    beta = Beta

    return alpha, beta
end

function fmoments(closure::Closure1D{Closure1DInput{T},Closure1DOutput{T}}, input::Closure1DInput{T}, Alpha::T, Beta::T)::Tuple{T,T,T} where {T<:AbstractFloat}

    l = input.l
    gE = input.gE
    ceffjl = input.ceffjl
    nspectral = input.nspectral
    nangular = input.nangular

    alpha, beta = ftf(closure, input, Alpha, Beta)

    if !isfinite(alpha) || alpha <= T(0)
        return typemax(T), typemax(T), typemax(T)
    end

    if !isfinite(beta) || beta <= -T(1) || beta >= T(1)
        return typemax(T), typemax(T), typemax(T)
    end

    if !isfinite(ceffjl) || ceffjl <= T(0)
        return typemax(T), typemax(T), typemax(T)
    end

    E1 = gE[l] * T(keV_to_J)
    E2 = gE[l+1] * T(keV_to_J)

    if E2 <= E1
        return typemax(T), typemax(T), typemax(T)
    end

    # Spectral quadrature
    dE = (E2 - E1) / T(nspectral)

    # Angular quadrature
    domega = T(4) * T(pi) / T(nangular)

    Ec = T(0)
    Fc = T(0)
    Pc = T(0)

    # Angular directions
    for k in 1:nangular
        omega = -T(1) + (T(k) - T(0.5)) * T(2) / T(nangular)

        for m in 1:nspectral
            Eph = E1 + (T(m) - T(0.5)) * dE
            nu = Eph / T(h)
            dnu = dE / T(h)
            argexp = (Eph / T(kB)) * alpha * (T(1) + beta * omega)

            if isnan(argexp) || argexp <= T(0)
                return typemax(T), typemax(T), typemax(T)
            end

            if isinf(argexp)
                Im = T(0)
            else
                den = exp(argexp) - T(1)

                if isinf(den)
                    Im = T(0)
                elseif !isfinite(den) || den <= T(0)
                    return typemax(T), typemax(T), typemax(T)
                else
                    Im = (T(2) * T(h) * nu^3 / T(c)^2) / den
                end
            end

            if !isfinite(Im) || Im < T(0)
                return typemax(T), typemax(T), typemax(T)
            end

            Ec += Im * dnu * domega / ceffjl
            Fc += omega * Im * dnu * domega
            Pc += omega * omega * Im * dnu * domega / ceffjl
        end
    end

    if !isfinite(Ec) || !isfinite(Fc) || !isfinite(Pc)
        return typemax(T), typemax(T), typemax(T)
    end

    if Ec <= T(0) || Pc <= T(0)
        return typemax(T), typemax(T), typemax(T)
    end

    return Ec, Fc, Pc
end

function fphi(closure::Closure1D{Closure1DInput{T},Closure1DOutput{T}}, input::Closure1DInput{T}, Alpha::T, Beta::T, fred::T)::T where {T<:AbstractFloat}

    tolsim = input.tolsim
    ceffjl = input.ceffjl
    Erjl = input.Erjl

    Ec, Fc, Pc = fmoments(closure, input, Alpha, Beta)

    if Ec == typemax(T) || Fc == typemax(T) || Pc == typemax(T)
        return typemax(T)
    end

    frc = abs(Fc) / (ceffjl * Ec)

    frc = corrector(closure, frc, T(0), T(1), tolsim, false)

    if frc <= eps(T) || frc > T(1)
        return typemax(T)
    end

    if fred <= eps(T)
        return typemax(T)
    end

    L1 = Base.log(Ec / Erjl)
    L2 = Base.log(frc / fred)

    J = L1^2 + L2^2

    if !isfinite(J)
        return typemax(T)
    end

    return J
end

function fm1gray(closure::Closure1D{Closure1DInput{T},Closure1DOutput{T}}, input::Closure1DInput{T}, fred::T)::Tuple{Bool,T} where {T<:AbstractFloat}

    tolsim = input.tolsim

    rd = T(4) - T(3) * fred^2

    if rd < T(0)
        return false, T(0)
    end

    fedd = (T(3) + T(4) * fred^2) / (T(5) + T(2) * sqrt(rd))

    if !isfinite(fedd)
        return false, T(0)
    end

    fedd = corrector(closure, fedd, T(1) / T(3), T(1), tolsim, false)

    if fedd < T(1) / T(3) || fedd > T(1)
        return false, T(0)
    end

    return true, fedd
end

function fmgpolynomial(closure::Closure1D{Closure1DInput{T},Closure1DOutput{T}}, input::Closure1DInput{T}, fred::T, fa::T, ta::T, sqfa::T)::Tuple{Bool,T} where {T<:AbstractFloat}

    tolsim = input.tolsim

    fedd = T(0)

    if (T(0) < fa <= T(0.9) && T(0) <= ta <= T(1.0e-5) * sqfa) || (T(0.9) <= fa <= T(1) && T(0) <= ta <= T(1.0e-3) / sqfa)
        f2 = fred^2
        f3 = fred^3
        f4 = fred^4
        f5 = fred^5
        f6 = fred^6
        f7 = fred^7
        f8 = fred^8
        f9 = fred^9
        f10 = fred^10
        f11 = fred^11

        fedd = T(1) / T(3) + T(2) / T(3) * f2 +
               T(0.3058945350580) * (f3 - f2) +
               T(-3.5960491961545) * (f4 - f2) +
               T(24.1130825450229) * (f5 - f2) +
               T(-92.2832025267787) * (f6 - f2) +
               T(220.1853600902954) * (f7 - f2) +
               T(-329.3727903123292) * (f8 - f2) +
               T(299.9957746394199) * (f9 - f2) +
               T(-151.3384873479037) * (f10 - f2) +
               T(32.2668624654447) * (f11 - f2)

    elseif (T(0) < fa <= T(0.9) && T(10.0^0.6) / sqfa <= ta) || (T(0.9) <= fa <= T(1) && T(10.0^0.6) / sqfa <= ta)
        f2 = fred^2
        f3 = fred^3
        f4 = fred^4
        f5 = fred^5
        f6 = fred^6
        f7 = fred^7
        f8 = fred^8
        f9 = fred^9
        f10 = fred^10
        f11 = fred^11
        f12 = fred^12

        fedd = T(1) / T(3) + f2 - T(1) / T(3) * f3 +
               T(6.5813886320063) * (f4 - T(2) * f3 + f2) +
               T(-44.5593808930324) * (f5 - T(3) * f3 + T(2) * f2) +
               T(179.9251066153469) * (f6 - T(4) * f3 + T(3) * f2) +
               T(-463.2145547920471) * (f7 - T(5) * f3 + T(4) * f2) +
               T(776.0741150675088) * (f8 - T(6) * f3 + T(5) * f2) +
               T(-841.0048410765069) * (f9 - T(7) * f3 + T(6) * f2) +
               T(566.6515584347650) * (f10 - T(8) * f3 + T(7) * f2) +
               T(-215.1853394583736) * (f11 - T(9) * f3 + T(8) * f2) +
               T(35.1130440991836) * (f12 - T(10) * f3 + T(9) * f2)

    elseif (T(0) < fa <= T(1.0e-3) && T(10.0^1.1) * sqfa <= ta <= T(10.0^-1.5) / sqfa) || (T(1.0e-3) <= fa <= T(10.0^-2.2) && T(10.0^1.1) * sqfa <= ta <= T(1))
        return fm1gray(closure, input, fred)

    else
        return false, T(0)
    end

    if !isfinite(fedd)
        return false, T(0)
    end

    fedd = corrector(closure, fedd, T(1) / T(3), T(1), tolsim, false)

    if fedd < T(1) / T(3) || fedd > T(1)
        return false, T(0)
    end

    return true, fedd
end

function fmgrandomsearch(closure::Closure1D{Closure1DInput{T},Closure1DOutput{T}}, input::Closure1DInput{T}, fred::T)::Tuple{Bool,T} where {T<:AbstractFloat}

    l = input.l
    tolsim = input.tolsim
    ceffjl = input.ceffjl
    Erjl = input.Erjl
    Frjl = input.Frjl
    hpclosure = input.hpclosure
    retry = input.retry
    reject = input.reject

    if length(hpclosure) != 10
        abort(closure, "Closure1DInput inconsistent hpclosure size for randomsearch", l, Erjl, Frjl, ceffjl, fred, T(0))
    end

    Alphamin = hpclosure[1]
    Alphamax = hpclosure[2]
    Betamin = hpclosure[3]
    Betamax = hpclosure[4]
    shrink = hpclosure[5]
    epsilon = hpclosure[6]
    nmax = Int(round(hpclosure[7]))
    stallmax = Int(round(hpclosure[8]))
    nsamp = Int(round(hpclosure[9]))
    seed0 = Int(round(hpclosure[10]))

    if Alphamax <= Alphamin
        abort(closure, "Closure1D invalid Alpha bounds for randomsearch", l, Erjl, Frjl, ceffjl, fred, T(0))
    end

    if Betamax <= Betamin
        abort(closure, "Closure1D invalid Beta bounds for randomsearch", l, Erjl, Frjl, ceffjl, fred, T(0))
    end

    if shrink <= T(0) || shrink >= T(1)
        abort(closure, "Closure1D invalid shrink for randomsearch", l, Erjl, Frjl, ceffjl, fred, T(0))
    end

    if epsilon <= T(0)
        abort(closure, "Closure1D invalid epsilon for randomsearch", l, Erjl, Frjl, ceffjl, fred, T(0))
    end

    if nmax < 1 || stallmax < 1 || nsamp < 1
        abort(closure, "Closure1D invalid integer parameter for randomsearch", l, Erjl, Frjl, ceffjl, fred, T(0))
    end

    dom = T[
        Alphamin Alphamax
        Betamin Betamax
    ]

    for itry in 0:retry
        seed = seed0 + itry
        objective = (x) -> fphi(closure, input, x[1], x[2], fred)

        _input = RandomSearchInput(
            dom,
            shrink,
            epsilon,
            nmax,
            stallmax,
            nsamp,
            seed,
            objective,
            reject
        )

        _solver = RandomSearchSolver{typeof(_input),RandomSearchOutput{T}}()
        output = kernel(_solver, _input)

        if reject.flag
            reset!(reject)

            if itry < retry
                continue
            end

            return false, T(0)
        end

        Alpha = output.out[1]
        Beta = output.out[2]

        Ec, Fc, Pc = fmoments(closure, input, Alpha, Beta)

        if Ec == typemax(T) || Fc == typemax(T) || Pc == typemax(T)
            if itry < retry
                continue
            end

            return false, T(0)
        end

        fedd = Pc / Ec

        if !isfinite(fedd)
            if itry < retry
                continue
            end

            return false, T(0)
        end

        fedd = corrector(closure, fedd, T(1) / T(3), T(1), tolsim, false)

        if fedd < T(1) / T(3) || fedd > T(1)
            if itry < retry
                continue
            end

            return false, T(0)
        end

        return true, fedd
    end

    return false, T(0)
end

function fmglinesearch(closure::Closure1D{Closure1DInput{T},Closure1DOutput{T}}, input::Closure1DInput{T}, fred::T)::Tuple{Bool,T} where {T<:AbstractFloat}

    l = input.l
    ceffjl = input.ceffjl
    Erjl = input.Erjl
    Frjl = input.Frjl

    #=
    NOTE: Linesearch is temporarily disabled for the multigroup M1 closure.
    The current objective is not suitable for finite-difference linesearch.
    =#

    abort(closure, "Closure1D fmglinesearch function not implemented", l, Erjl, Frjl, ceffjl, fred, T(0))
end

function fmgai(closure::Closure1D{Closure1DInput{T},Closure1DOutput{T}}, input::Closure1DInput{T}, fred::T, ta::T, fa::T)::Tuple{Bool,T} where {T<:AbstractFloat}

    l = input.l
    tolsim = input.tolsim
    ceffjl = input.ceffjl
    Erjl = input.Erjl
    Frjl = input.Frjl
    fclosure = input.fclosure

    if isempty(strip(fclosure)) || fclosure == "none"
        abort(closure, "Closure1D invalid fclosure for ai", l, Erjl, Frjl, ceffjl, fred, T(0))
    end

    key = (:closure_ai, T, abspath(fclosure))

    ai = nothing

    lock(CMEM.lock)
    try
        if !haskey(CMEM.data, key)
            info = JSON.parsefile(fclosure)

            fonnx = joinpath(dirname(fclosure), info["onnx"])
            session = ONNXRunTime.load_inference(fonnx, execution_provider=:cpu)

            Xnames = info["input"]
            Yname = String(info["output"])
            normalization = info["normalization"]

            CMEM.data[key] = (
                session=session,
                Xnames=Xnames,
                Yname=Yname,
                Xonnx=String(ONNXRunTime.input_names(session)[1]),
                Yonnx=Yname,
                Xmin=T[normalization["min"][name] for name in Xnames],
                Xmax=T[normalization["max"][name] for name in Xnames],
            )
        end

        ai = CMEM.data[key]
    finally
        unlock(CMEM.lock)
    end

    r = ceffjl / T(c)

    X = T[
        fred,
        log10(ta),
        log10(fa),
        log10(r)
    ]

    Xnorm = Float32.((X .- ai.Xmin) ./ (ai.Xmax .- ai.Xmin))
    Xbatch = reshape(Xnorm, 1, :)

    Y = ai.session(Dict(ai.Xonnx => Xbatch), [ai.Yonnx])
    fedd = T(Y[ai.Yonnx][1])

    if !isfinite(fedd)
        return false, T(0)
    end

    fedd = corrector(closure, fedd, T(1) / T(3), T(1), tolsim, true)

    if fedd < T(1) / T(3) || fedd > T(1)
        return false, T(0)
    end

    return true, fedd
end

function closure!(closure::Closure1D{Closure1DInput{T},Closure1DOutput{T}}, input::Closure1DInput{T})::Closure1DOutput{T} where {T<:AbstractFloat}

    l = input.l
    tolsim = input.tolsim
    gE = input.gE
    ceffjl = input.ceffjl
    Erjl = input.Erjl
    Frjl = input.Frjl
    mclosure = input.mclosure
    reject = input.reject

    reset!(reject)

    fred = abs(Frjl) / (ceffjl * Erjl)

    fred = corrector(closure, fred, T(0), T(1), tolsim, false)

    if fred < T(0) || fred > T(1)
        reject!(reject, "Closure1D invalid | fred = $(fred)")
        return Closure1DOutput(T(0))
    end

    if fred <= tolsim
        return Closure1DOutput(T(1) / T(3))
    end

    E1 = gE[l] * T(keV_to_J)
    E2 = gE[l+1] * T(keV_to_J)

    nu1 = E1 / T(h)
    nu2 = E2 / T(h)

    fa = nu1 / nu2

    if !isfinite(fa) || fa <= T(0) || fa > T(1)
        abort(closure, "Closure1D invalid group affinity", l, Erjl, Frjl, ceffjl, fred, T(0))
    end

    ta = (T(kB) / T(h)) * (Erjl / T(aR))^T(0.25) / sqrt(nu1 * nu2)

    if !isfinite(ta) || ta <= T(0)
        abort(closure, "Closure1D invalid adimensional temperature", l, Erjl, Frjl, ceffjl, fred, T(0))
    end

    sqfa = sqrt(fa)

    found, fedd = fmgpolynomial(closure, input, fred, fa, ta, sqfa)

    if found
        return Closure1DOutput(fedd)
    end

    if mclosure == "randomsearch"
        found, fedd = fmgrandomsearch(closure, input, fred)

    elseif mclosure == "linesearch"
        found, fedd = fmglinesearch(closure, input, fred)

    elseif mclosure == "ai"
        found, fedd = fmgai(closure, input, fred, ta, fa)

    else
        abort(closure, "Closure1DInput invalid mclosure", l, Erjl, Frjl, ceffjl, fred, T(0))
    end

    if found
        return Closure1DOutput(fedd)
    end

    #=
    NOTE: M1-gray fallback used when the multigroup closure fails.
    =#
    
    found, fedd = fm1gray(closure, input, fred)

    if found
        return Closure1DOutput(fedd)
    end

    if !found
        reject!(reject, "Closure1D invalid | found = $(found)")
        return Closure1DOutput(T(0))
    end

    return Closure1DOutput(T(0))
end
