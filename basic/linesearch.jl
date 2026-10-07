include("../interface/isolver.jl")
include("reject.jl")

using Random
using LinearAlgebra
using Base.Threads
using Optim

#===================================================================
                        Line Search solver
===================================================================#

const USE_LS_BRENT = lowercase(get(ENV, "USE_LS_BRENT", "false")) in ("1", "true")

struct LineSearchInput{T,F} <: ISolverInput
    dom::Array{T,2}
    amin::T
    amax::T
    epsilon::T
    nmax::Int
    cgmode::String
    delta::T
    nsamp::Int
    seed::Int
    objective::F
    reject::Reject
end

struct LineSearchOutput{T} <: ISolverOutput
    out::Array{T,1}
end

struct LineSearchSolver{Tin<:ISolverInput,Tout<:ISolverOutput} <: ISolver{Tin,Tout} end

function abort(solver::LineSearchSolver{LineSearchInput{T,F},LineSearchOutput{T}}, msg::String, dom::Array{T,2}, amin::T, amax::T, epsilon::T, nmax::Int, cgmode::String, delta::T, nsamp::Int, seed::Int) where {T<:AbstractFloat,F<:Function}
    log("ERROR: LineSearch")
    log("solver  = $(typeof(solver))")
    log("dom     = $(dom)")
    log("amin    = $(amin)")
    log("amax    = $(amax)")
    log("epsilon = $(epsilon)")
    log("nmax    = $(nmax)")
    log("cgmode  = $(cgmode)")
    log("delta   = $(delta)")
    log("nsamp   = $(nsamp)")
    log("seed    = $(seed)")

    error("$(msg)")
end

function kernel(solver::LineSearchSolver{LineSearchInput{T,F},LineSearchOutput{T}}, input::LineSearchInput{T,F})::LineSearchOutput{T} where {T<:AbstractFloat,F<:Function}

    dom = input.dom
    amin = input.amin
    amax = input.amax
    epsilon = input.epsilon
    nmax = input.nmax
    cgmode = input.cgmode
    delta = input.delta
    nsamp = input.nsamp
    seed = input.seed
    objective = input.objective
    reject = input.reject

    reset!(reject)

    ndim = size(dom, 1)

    # Golden-section line search
    fls1 = (x, p, amin, amax, epsilon, nmax, objective) -> begin
        out = similar(x)

        _amin = amin
        _amax = amax

        r = T((sqrt(5.0) - 1.0) / 2.0)

        a1 = _amax - r * (_amax - _amin)
        a2 = _amin + r * (_amax - _amin)

        out .= x .+ a1 .* p
        phi1 = objective(out)

        out .= x .+ a2 .* p
        phi2 = objective(out)

        for n in 1:nmax
            if _amax - _amin <= epsilon
                break
            end

            if phi1 < phi2
                _amax = a2
                a2 = a1
                phi2 = phi1

                a1 = _amax - r * (_amax - _amin)
                out .= x .+ a1 .* p
                phi1 = objective(out)
            else
                _amin = a1
                a1 = a2
                phi1 = phi2

                a2 = _amin + r * (_amax - _amin)
                out .= x .+ a2 .* p
                phi2 = objective(out)
            end
        end

        a = (_amin + _amax) / T(2)
        out .= x .+ a .* p

        return out
    end

    # Brent line search
    fls2 = (x, p, amin, amax, objective) -> begin
        out = similar(x)
        tmp = similar(x)

        phi = a -> begin
            tmp .= x .+ a .* p
            return objective(tmp)
        end

        result = optimize(phi, amin, amax, Brent())
        a = T(Optim.minimizer(result))

        out .= x .+ a .* p

        return out
    end

    nthd = Threads.nthreads()

    best_y = fill(typemax(T), nthd)
    best_x = [zeros(T, ndim) for _ in 1:nthd]
    best_accepted = fill(false, nthd)

    Threads.@threads for ithd in 1:nthd
        for isamp in ithd:nthd:nsamp
            rng = MersenneTwister(seed + isamp)

            x = zeros(T, ndim)
            out = zeros(T, ndim)

            grad = zeros(T, ndim)
            grad0 = zeros(T, ndim)

            p = zeros(T, ndim)
            p0 = zeros(T, ndim)

            dgrad = zeros(T, ndim)

            accepted = false

            x .= dom[:, 1] .+ rand(rng, T, ndim) .* (dom[:, 2] .- dom[:, 1])

            for n in 1:nmax
                for j in 1:ndim
                    xj = x[j]
                    x[j] = xj - delta
                    phi1 = objective(x)
                    x[j] = xj + delta
                    phi2 = objective(x)
                    x[j] = xj

                    grad[j] = (phi2 - phi1) / (T(2) * delta)
                end

                if n == 1 || cgmode == "none"
                    beta = zero(T)

                elseif cgmode == "FR"
                    dgrad .= grad .- grad0
                    beta = dot(grad, grad) / max(dot(grad0, grad0), eps(T))

                elseif cgmode == "PR"
                    dgrad .= grad .- grad0
                    beta = dot(grad, dgrad) / max(dot(grad0, grad0), eps(T))

                elseif cgmode == "PR+"
                    dgrad .= grad .- grad0
                    beta = max(zero(T), dot(grad, dgrad) / max(dot(grad0, grad0), eps(T)))

                elseif cgmode == "HS"
                    dgrad .= grad .- grad0
                    beta = dot(grad, dgrad) / max(dot(p0, dgrad), eps(T))

                elseif cgmode == "DY"
                    dgrad .= grad .- grad0
                    beta = dot(grad, grad) / max(dot(p0, dgrad), eps(T))

                else
                    abort(solver, "LineSearchInput invalid cgmode", dom, amin, amax, epsilon, nmax, cgmode, delta, nsamp, seed)
                end

                p .= -grad .+ beta .* p0

                if dot(grad, p) >= zero(T)
                    p .= -grad
                end

                if USE_LS_BRENT
                    out .= fls2(x, p, amin, amax, objective)
                else
                    out .= fls1(x, p, amin, amax, epsilon, nmax, objective)
                end

                if norm(out - x) < epsilon
                    x .= out
                    accepted = true
                    break
                end

                grad0 .= grad
                p0 .= p
                x .= out
            end

            y = objective(x)

            if y < best_y[ithd]
                best_y[ithd] = y
                best_x[ithd] .= x
                best_accepted[ithd] = accepted
            end
        end
    end

    ibest = argmin(best_y)
    out = copy(best_x[ibest])

    if best_y[ibest] == typemax(T) || !best_accepted[ibest]
        reject!(reject, "LineSearch invalid | best_y[ibest] = $(best_y[ibest]) best_accepted[ibest] = $(best_accepted[ibest])")
    end

    return LineSearchOutput(out)
end

#===================================================================
                        Line Search test
===================================================================#

function test_line_search()
    println("\n============= Line Search test ===============")

    dom = [
        -10.0 10.0;
        -10.0 10.0
    ]

    amin = 0.0
    amax = 1.0
    epsilon = 1.0e-9
    nmax = 100
    delta = 1.0e-9
    nsamp = 1000
    seed = 42

    objective = (x) -> -3.5 * exp(-((x[1] - 3.0)^2 + (x[2] - 3.0)^2)) -
                       2.8 * exp(-((x[1] - 1.0)^2 + (x[2] + 3.0)^2)) -
                       2.4 * exp(-((x[1] + 2.0)^2 + (x[2] - 1.0)^2)) -
                       2.2 * exp(-((x[1] - 4.5)^2 + (x[2] + 4.0)^2)) -
                       2.0 * exp(-((x[1] + 4.0)^2 + (x[2] + 4.5)^2)) -
                       1.9 * exp(-((x[1] - 6.0)^2 + (x[2] - 1.5)^2)) -
                       1.7 * exp(-((x[1] + 5.5)^2 + (x[2] - 5.0)^2)) -
                       1.6 * exp(-((x[1] - 0.0)^2 + (x[2] - 6.0)^2))

    println("dom          = ", dom)
    println("amin         = ", amin)
    println("amax         = ", amax)
    println("epsilon      = ", epsilon)
    println("nmax         = ", nmax)
    println("delta        = ", delta)
    println("nsamp        = ", nsamp)
    println("seed         = ", seed)
    println("nthd         = ", Threads.nthreads())
    println("USE_LS_BRENT = ", USE_LS_BRENT)

    cgmodes = ["none", "FR", "PR", "PR+", "HS", "DY"]

    for cgmode in cgmodes
        println("\n------------- cgmode = ", cgmode, " -------------")

        reject = Reject(false, "")
        reset!(reject)

        _solver = LineSearchSolver{LineSearchInput{Float64,typeof(objective)},LineSearchOutput{Float64}}()
        input = LineSearchInput(dom, amin, amax, epsilon, nmax, cgmode, delta, nsamp, seed, objective, reject)

        kout = kernel(_solver, input)

        println("kernel(...) out       = ", kout.out)
        println("kernel(...) objective = ", objective(kout.out))
        println("kernel(...) reject    = ", reject.flag)

        reset!(reject)

        sout = solver(_solver, input)

        println("solver(...) out       = ", sout.out)
        println("solver(...) objective = ", objective(sout.out))
        println("solver(...) reject    = ", reject.flag)
    end
end

function test()
    test_line_search()
end

if abspath(PROGRAM_FILE) == @__FILE__
    test()
end
