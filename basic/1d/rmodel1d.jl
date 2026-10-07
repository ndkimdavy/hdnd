include("../../interface/isolver.jl")
include("../reject.jl")
include("closure1d.jl")

#===================================================================
                            RModel 1D
               (Er, Fr, kR, ceff) -> (Er, Fr, Pr, fedd)
===================================================================#

struct RModel1DInput{T} <: ISolverInput
    l::Int
    mode::String
    tolsim::T
    gE::Array{T,1}
    kRjl::T
    ceffjl::T
    nspectral::Int
    nangular::Int
    rhoj::T
    Thj::T
    Erjl::T
    Frjl::T
    gradErjl::T
    gradFrjl::T
    mclosure::String
    hpclosure::Array{T,1}
    fclosure::String
    hprmodel::Array{T,1}
    retry::Int
    reject::Reject
end

struct RModel1DOutput{T} <: ISolverOutput
    Erjl::T
    Frjl::T
    Prjl::T
    fedd::T
    regime::Union{Nothing,Symbol}
end

struct RModel1DSolver{Tin<:ISolverInput,Tout<:ISolverOutput} <: ISolver{Tin,Tout} end

function abort(solver::RModel1DSolver{RModel1DInput{T},RModel1DOutput{T}},
    msg::String,
    l::Int,
    tolsim::T,
    nspectral::Int, nangular::Int,
    rhoj::T, Thj::T,
    Erjl::T, Frjl::T,
    gradErjl::T, gradFrjl::T,
    kRjl::T, ceffjl::T,
    fred::T, fedd::T, mfreejl::T) where {T<:AbstractFloat}

    log("ERROR: RModel1D")
    log("solver    = $(typeof(solver))")
    log("l         = $(l)")
    log("tolsim    = $(tolsim)")
    log("nspectral = $(nspectral)")
    log("nangular  = $(nangular)")
    log("rhoj      = $(rhoj)")
    log("Thj       = $(Thj)")
    log("Erjl      = $(Erjl)")
    log("Frjl      = $(Frjl)")
    log("gradErjl  = $(gradErjl)")
    log("gradFrjl  = $(gradFrjl)")
    log("kRjl      = $(kRjl)")
    log("ceffjl    = $(ceffjl)")
    log("fred      = $(fred)")
    log("fedd      = $(fedd)")
    log("mfreejl   = $(mfreejl)")

    error("$(msg)")
end

function fthin(solver::RModel1DSolver{RModel1DInput{T},RModel1DOutput{T}}, input::RModel1DInput{T}, ceffjl::T)::RModel1DOutput{T} where {T<:AbstractFloat}

    tolsim = input.tolsim
    Erjl = input.Erjl
    Frjl = input.Frjl
    reject = input.reject

    fred = abs(Frjl) / (ceffjl * Erjl)

    if fred < T(0) || fred > T(1)
        reject!(reject, "RModel1D invalid | fred = $(fred)")
        return RModel1DOutput(T(0), T(0), T(0), T(0), nothing)
    end

    if fred <= tolsim
        fedd = T(1) / T(3)
    elseif fred >= T(1) - tolsim
        fedd = T(1)
    else
        output = fm1closure(solver, input, ceffjl)

        if reject.flag
            return output
        end

        return RModel1DOutput(output.Erjl, output.Frjl, output.Prjl, output.fedd, :thin)
    end

    Prjl = fedd * Erjl

    if !isfinite(Prjl) || Prjl < T(0)
        reject!(reject, "RModel1D invalid | Prjl = $(Prjl)")
        return RModel1DOutput(T(0), T(0), T(0), T(0), nothing)
    end

    return RModel1DOutput(Erjl, Frjl, Prjl, fedd, :thin)
end

function fthick(solver::RModel1DSolver{RModel1DInput{T},RModel1DOutput{T}}, input::RModel1DInput{T}, ceffjl::T, kRjl::T)::RModel1DOutput{T} where {T<:AbstractFloat}

    tolsim = input.tolsim
    rhoj = input.rhoj
    Erjl = input.Erjl
    gradErjl = input.gradErjl
    reject = input.reject

    fedd = T(1) / T(3)
    Prjl = fedd * Erjl

    if !isfinite(Prjl) || Prjl < T(0)
        reject!(reject, "RModel1D invalid | Prjl = $(Prjl)")
        return RModel1DOutput(T(0), T(0), T(0), T(0), nothing)
    end

    # Radiative diffusive heat flux equation
    FDjl = -ceffjl * gradErjl / (T(3) * rhoj * kRjl)

    if !isfinite(FDjl)
        reject!(reject, "RModel1D invalid | FDjl = $(FDjl)")
        return RModel1DOutput(T(0), T(0), T(0), T(0), nothing)
    end

    if abs(FDjl) <= tolsim
        Frjl = T(0)
    else
        scale = ceffjl * Erjl / abs(FDjl)

        if !isfinite(scale)
            reject!(reject, "RModel1D invalid | scale = $(scale)")
            return RModel1DOutput(T(0), T(0), T(0), T(0), nothing)
        end

        scale = corrector(solver, scale, T(0), T(1), tolsim, true)

        if scale < T(0) || scale > T(1)
            reject!(reject, "RModel1D invalid | scale = $(scale)")
            return RModel1DOutput(T(0), T(0), T(0), T(0), nothing)
        end

        Frjl = scale * FDjl
    end

    if !isfinite(Frjl)
        reject!(reject, "RModel1D invalid | Frjl = $(Frjl)")
        return RModel1DOutput(T(0), T(0), T(0), T(0), nothing)
    end

    return RModel1DOutput(Erjl, Frjl, Prjl, fedd, :thick)
end

function fm1(solver::RModel1DSolver{RModel1DInput{T},RModel1DOutput{T}}, input::RModel1DInput{T}, ceffjl::T, fedd::T)::RModel1DOutput{T} where {T<:AbstractFloat}

    Erjl = input.Erjl
    Frjl = input.Frjl
    reject = input.reject

    if !isfinite(fedd) || fedd < T(1) / T(3) || fedd > T(1)
        reject!(reject, "RModel1D invalid | fedd = $(fedd)")
        return RModel1DOutput(T(0), T(0), T(0), T(0), nothing)
    end

    Prjl = fedd * Erjl

    if !isfinite(Prjl) || Prjl < T(0)
        reject!(reject, "RModel1D invalid | Prjl = $(Prjl)")
        return RModel1DOutput(T(0), T(0), T(0), T(0), nothing)
    end

    return RModel1DOutput(Erjl, Frjl, Prjl, fedd, :m1)
end

function fm1closure(solver::RModel1DSolver{RModel1DInput{T},RModel1DOutput{T}}, input::RModel1DInput{T}, ceffjl::T)::RModel1DOutput{T} where {T<:AbstractFloat}

    l = input.l
    tolsim = input.tolsim
    gE = input.gE
    nspectral = input.nspectral
    nangular = input.nangular
    Erjl = input.Erjl
    Frjl = input.Frjl
    mclosure = input.mclosure
    hpclosure = input.hpclosure
    fclosure = input.fclosure
    retry = input.retry
    reject = input.reject

    closure = Closure1D{Closure1DInput{T},Closure1DOutput{T}}()

    _input = Closure1DInput(
        l,
        tolsim,
        gE,
        ceffjl,
        nspectral,
        nangular,
        Erjl,
        Frjl,
        mclosure,
        hpclosure,
        fclosure,
        retry,
        reject
    )

    output = closure!(closure, _input)
    fedd = output.fedd

    if reject.flag
        return RModel1DOutput(T(0), T(0), T(0), T(0), nothing)
    end

    return fm1(solver, input, ceffjl, fedd)
end

function kernel(solver::RModel1DSolver{RModel1DInput{T},RModel1DOutput{T}}, input::RModel1DInput{T})::RModel1DOutput{T} where {T<:AbstractFloat}

    mode = input.mode
    hprmodel = input.hprmodel
    kRjl = input.kRjl
    ceffjl = input.ceffjl
    rhoj = input.rhoj
    reject = input.reject

    reset!(reject)

    if mode == "classic"
        return fm1closure(solver, input, ceffjl)
    end

    mfreethick = hprmodel[1]
    mfreethin = hprmodel[2]

    mfreejl = T(1) / (rhoj * kRjl)

    if !isfinite(mfreejl) || mfreejl <= T(0)
        reject!(reject, "RModel1D invalid | mfreejl = $(mfreejl)")
        return RModel1DOutput(T(0), T(0), T(0), T(0), nothing)
    end

    if mfreejl <= mfreethick
        return fthick(solver, input, ceffjl, kRjl)

    elseif mfreejl >= mfreethin
        return fthin(solver, input, ceffjl)

    else
        return fm1closure(solver, input, ceffjl)
    end
end
