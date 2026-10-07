include("../interface/isolver.jl")
include("reject.jl")

using Random
using LinearAlgebra
using Base.Threads

#===================================================================
                        Random Search solver
===================================================================#

struct RandomSearchInput{T,F} <: ISolverInput
    dom::Array{T,2}
    shrink::T
    epsilon::T
    nmax::Int
    stallmax::Int
    nsamp::Int
    seed::Int
    objective::F
    reject::Reject
end

struct RandomSearchOutput{T} <: ISolverOutput
    out::Array{T,1}
end

struct RandomSearchSolver{Tin<:ISolverInput,Tout<:ISolverOutput} <: ISolver{Tin,Tout} end

function abort(solver::RandomSearchSolver{RandomSearchInput{T,F},RandomSearchOutput{T}}, msg::String, dom::Array{T,2}, shrink::T, epsilon::T, nmax::Int, stallmax::Int, nsamp::Int, seed::Int) where {T<:AbstractFloat,F<:Function}
    log("ERROR: RandomSearch")
    log("solver   = $(typeof(solver))")
    log("dom      = $(dom)")
    log("shrink   = $(shrink)")
    log("epsilon  = $(epsilon)")
    log("nmax     = $(nmax)")
    log("stallmax = $(stallmax)")
    log("nsamp    = $(nsamp)")
    log("seed     = $(seed)")

    error("$(msg)")
end

function kernel(solver::RandomSearchSolver{RandomSearchInput{T,F},RandomSearchOutput{T}}, input::RandomSearchInput{T,F})::RandomSearchOutput{T} where {T<:AbstractFloat,F<:Function}

    dom = input.dom
    shrink = input.shrink
    epsilon = input.epsilon
    nmax = input.nmax
    stallmax = input.stallmax
    nsamp = input.nsamp
    seed = input.seed
    objective = input.objective
    reject = input.reject

    reset!(reject)

    ndim = size(dom, 1)
    nthd = Threads.nthreads()

    _dom = copy(dom)

    xbest = zeros(T, ndim)
    ybest = typemax(T)

    stall = 0
    accepted = false

    for n in 1:nmax
        yold = ybest

        best_y = fill(typemax(T), nthd)
        best_x = [zeros(T, ndim) for _ in 1:nthd]

        Threads.@threads for ithd in 1:nthd
            for isamp in ithd:nthd:nsamp
                rng = MersenneTwister(seed + n + isamp)

                x = zeros(T, ndim)

                x .= _dom[:, 1] .+ rand(rng, T, ndim) .* (_dom[:, 2] .- _dom[:, 1])

                y = objective(x)

                if y < best_y[ithd]
                    best_y[ithd] = y
                    best_x[ithd] .= x
                end
            end
        end

        ibest = argmin(best_y)

        if best_y[ibest] < ybest
            ybest = best_y[ibest]
            xbest .= best_x[ibest]
            accepted = true
        end

        if abs(yold - ybest) <= epsilon
            stall += 1
        else
            stall = 0
        end

        if stall >= stallmax
            break
        end

        for j in 1:ndim
            width = _dom[j, 2] - _dom[j, 1]
            hwidth = T(0.5) * shrink * width

            _dom[j, 1] = max(dom[j, 1], xbest[j] - hwidth)
            _dom[j, 2] = min(dom[j, 2], xbest[j] + hwidth)
        end

        if maximum(_dom[:, 2] .- _dom[:, 1]) <= epsilon
            break
        end
    end

    out = copy(xbest)

    if ybest == typemax(T) || !accepted
        reject!(reject, "RandomSearch invalid | ybest = $(ybest) accepted = $(accepted)")
    end

    return RandomSearchOutput(out)
end

#===================================================================
                        Random Search test
===================================================================#

function test_random_search()
    println("\n============ Random Search test ==============")

    dom = [
        -10.0 10.0;
        -10.0 10.0
    ]

    shrink = 0.5
    epsilon = 1.0e-9
    nmax = 100
    stallmax = 5
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

    println("dom         = ", dom)
    println("shrink      = ", shrink)
    println("epsilon     = ", epsilon)
    println("nmax        = ", nmax)
    println("stallmax    = ", stallmax)
    println("nsamp       = ", nsamp)
    println("seed        = ", seed)
    println("nthd        = ", Threads.nthreads())

    reject = Reject(false, "")
    reset!(reject)

    _solver = RandomSearchSolver{RandomSearchInput{Float64,typeof(objective)},RandomSearchOutput{Float64}}()
    input = RandomSearchInput(dom, shrink, epsilon, nmax, stallmax, nsamp, seed, objective, reject)

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

function test()
    test_random_search()
end

if abspath(PROGRAM_FILE) == @__FILE__
    test()
end
