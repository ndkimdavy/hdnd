include("../interface/isolver.jl")
include("reject.jl")

#===================================================================
                         Rusanov solver
          out = 1/2 * (FL + FR) - 1/2 * smax * (UR - UL)
====================================================================#

struct RusanovInput{T} <: ISolverInput
    UL::Array{T,1}
    UR::Array{T,1}
    FL::Array{T,1}
    FR::Array{T,1}
    smax::T
    reject::Reject
end

struct RusanovOutput{T} <: ISolverOutput
    out::Array{T,1}
end

struct RusanovSolver{T} <: ISolver{RusanovInput{T},RusanovOutput{T}} end

function abort(solver::RusanovSolver{T}, msg::String, UL::Array{T,1}, UR::Array{T,1}, FL::Array{T,1}, FR::Array{T,1}, smax::T) where {T<:AbstractFloat}
    log("ERROR: Rusanov")
    log("solver = $(typeof(solver))")
    log("UL     = $(UL)")
    log("UR     = $(UR)")
    log("FL     = $(FL)")
    log("FR     = $(FR)")
    log("smax   = $(smax)")

    error("$(msg)")
end

function kernel(solver::RusanovSolver{T}, input::RusanovInput{T})::RusanovOutput{T} where {T<:AbstractFloat}

    UL = input.UL
    UR = input.UR
    FL = input.FL
    FR = input.FR
    smax = input.smax
    reject = input.reject

    reset!(reject)

    half = T(0.5)
    out = half .* (FL .+ FR) .- half .* smax .* (UR .- UL)

    if any(x -> !isfinite(x), out)
        reject!(reject, "Rusanov invalid | out = $(out)")
    end

    return RusanovOutput(out)
end

#===================================================================
                            HLL solver
    if 0 <= SL
        out = FL
    elseif SL <= 0 <= SR
        out = (SR * FL - SL * FR + SL * SR * (UR - UL)) / (SR - SL)
    else
        out = FR
====================================================================#

struct HLLInput{T} <: ISolverInput
    UL::Array{T,1}
    UR::Array{T,1}
    FL::Array{T,1}
    FR::Array{T,1}
    SL::T
    SR::T
    reject::Reject
end

struct HLLOutput{T} <: ISolverOutput
    out::Array{T,1}
end

struct HLLSolver{T} <: ISolver{HLLInput{T},HLLOutput{T}} end

function abort(solver::HLLSolver{T}, msg::String, UL::Array{T,1}, UR::Array{T,1}, FL::Array{T,1}, FR::Array{T,1}, SL::T, SR::T) where {T<:AbstractFloat}
    log("ERROR: HLL")
    log("solver = $(typeof(solver))")
    log("UL     = $(UL)")
    log("UR     = $(UR)")
    log("FL     = $(FL)")
    log("FR     = $(FR)")
    log("SL     = $(SL)")
    log("SR     = $(SR)")

    error("$(msg)")
end

function kernel(solver::HLLSolver{T}, input::HLLInput{T})::HLLOutput{T} where {T<:AbstractFloat}

    UL = input.UL
    UR = input.UR
    FL = input.FL
    FR = input.FR
    SL = input.SL
    SR = input.SR
    reject = input.reject

    reset!(reject)

    zeroT = T(0)

    if zeroT <= SL
        out = FL
    elseif SL <= zeroT <= SR
        out = (SR .* FL .- SL .* FR .+ (SL * SR) .* (UR .- UL)) ./ (SR - SL)
    else
        out = FR
    end

    if any(x -> !isfinite(x), out)
        reject!(reject, "HLL invalid | out = $(out)")
    end

    return HLLOutput(out)
end

#===================================================================
                            HLLE solver
    Choice of wave speeds for Euler:
        SL = min(uL - cL, uR - cR)
        SR = max(uL + cL, uR + cR)
    if 0 <= SL
        out = FL
    elseif SL <= 0 <= SR
        out = (SR * FL - SL * FR + SL * SR * (UR - UL)) / (SR - SL)
    else
        out = FR
====================================================================#

struct HLLEInput{T} <: ISolverInput
    UL::Array{T,1}
    UR::Array{T,1}
    FL::Array{T,1}
    FR::Array{T,1}
    SL::T
    SR::T
    reject::Reject
end

struct HLLEOutput{T} <: ISolverOutput
    out::Array{T,1}
end

struct HLLESolver{T} <: ISolver{HLLEInput{T},HLLEOutput{T}} end

function abort(solver::HLLESolver{T}, msg::String, UL::Array{T,1}, UR::Array{T,1}, FL::Array{T,1}, FR::Array{T,1}, SL::T, SR::T) where {T<:AbstractFloat}
    log("ERROR: HLLE")
    log("solver = $(typeof(solver))")
    log("UL     = $(UL)")
    log("UR     = $(UR)")
    log("FL     = $(FL)")
    log("FR     = $(FR)")
    log("SL     = $(SL)")
    log("SR     = $(SR)")

    error("$(msg)")
end

function kernel(solver::HLLESolver{T}, input::HLLEInput{T})::HLLEOutput{T} where {T<:AbstractFloat}

    UL = input.UL
    UR = input.UR
    FL = input.FL
    FR = input.FR
    SL = input.SL
    SR = input.SR
    reject = input.reject

    reset!(reject)

    zeroT = T(0)

    if zeroT <= SL
        out = FL
    elseif SL <= zeroT <= SR
        out = (SR .* FL .- SL .* FR .+ (SL * SR) .* (UR .- UL)) ./ (SR - SL)
    else
        out = FR
    end

    if any(x -> !isfinite(x), out)
        reject!(reject, "HLLE invalid | out = $(out)")
    end

    return HLLEOutput(out)
end

#===================================================================
                            HLLC solver
    Choice of wave speeds for Euler:
        SL = min(uL - cL, uR - cR)
        SR = max(uL + cL, uR + cR)
    Contact wave speed:
        SM = (pR - pL + rhoL * uL * (SL - uL) - rhoR * uR * (SR - uR)) /
             (rhoL * (SL - uL) - rhoR * (SR - uR))
    if 0 <= SL
        out = FL
    elseif SL <= 0 <= SM
        out = FL + SL * (USL - UL)
    elseif SM <= 0 <= SR
        out = FR + SR * (USR - UR)
    else
        out = FR
====================================================================#

struct HLLCInput{T} <: ISolverInput
    UL::Array{T,1}
    UR::Array{T,1}
    FL::Array{T,1}
    FR::Array{T,1}
    USL::Array{T,1}
    USR::Array{T,1}
    SL::T
    SM::T
    SR::T
    reject::Reject
end

struct HLLCOutput{T} <: ISolverOutput
    out::Array{T,1}
end

struct HLLCSolver{T} <: ISolver{HLLCInput{T},HLLCOutput{T}} end

function abort(solver::HLLCSolver{T}, msg::String, UL::Array{T,1}, UR::Array{T,1}, FL::Array{T,1}, FR::Array{T,1}, USL::Array{T,1}, USR::Array{T,1}, SL::T, SM::T, SR::T) where {T<:AbstractFloat}
    log("ERROR: HLLC")
    log("solver = $(typeof(solver))")
    log("UL     = $(UL)")
    log("UR     = $(UR)")
    log("FL     = $(FL)")
    log("FR     = $(FR)")
    log("USL    = $(USL)")
    log("USR    = $(USR)")
    log("SL     = $(SL)")
    log("SM     = $(SM)")
    log("SR     = $(SR)")

    error("$(msg)")
end

function kernel(solver::HLLCSolver{T}, input::HLLCInput{T})::HLLCOutput{T} where {T<:AbstractFloat}

    UL = input.UL
    UR = input.UR
    FL = input.FL
    FR = input.FR
    USL = input.USL
    USR = input.USR
    SL = input.SL
    SM = input.SM
    SR = input.SR
    reject = input.reject

    reset!(reject)

    zeroT = T(0)

    if zeroT <= SL
        out = FL
    elseif SL <= zeroT <= SM
        out = FL .+ SL .* (USL .- UL)
    elseif SM <= zeroT <= SR
        out = FR .+ SR .* (USR .- UR)
    else
        out = FR
    end

    if any(x -> !isfinite(x), out)
        reject!(reject, "HLLC invalid | out = $(out)")
    end

    return HLLCOutput(out)
end

#===================================================================
                                TESTS
====================================================================#

function test_rusanov()
    println("\n================ Rusanov test ================")

    UL = [1.0, 2.0, 3.0]
    UR = [1.5, 2.5, 3.5]
    FL = [10.0, 20.0, 30.0]
    FR = [40.0, 50.0, 60.0]
    smax = 2.0

    println("UL   = ", UL)
    println("UR   = ", UR)
    println("FL   = ", FL)
    println("FR   = ", FR)
    println("smax = ", smax)

    reject = Reject(false, "")
    reset!(reject)

    _solver = RusanovSolver{Float64}()
    input = RusanovInput(UL, UR, FL, FR, smax, reject)

    println("kernel(...) out    = ", kernel(_solver, input).out)
    println("kernel(...) reject = ", reject.flag)

    reset!(reject)

    println("solver(...) out    = ", solver(_solver, input).out)
    println("solver(...) reject = ", reject.flag)
end

function test_hll()
    println("\n================== HLL test ==================")

    UL = [1.0, 2.0, 3.0]
    UR = [1.5, 2.5, 3.5]
    FL = [10.0, 20.0, 30.0]
    FR = [40.0, 50.0, 60.0]

    uL = 1.0
    cL = 2.0
    uR = -0.5
    cR = 1.5

    SL = min(uL - cL, uR - cR)
    SR = max(uL + cL, uR + cR)

    println("UL = ", UL)
    println("UR = ", UR)
    println("FL = ", FL)
    println("FR = ", FR)
    println("uL = ", uL, ", cL = ", cL)
    println("uR = ", uR, ", cR = ", cR)
    println("SL = ", SL)
    println("SR = ", SR)

    reject = Reject(false, "")
    reset!(reject)

    _solver = HLLSolver{Float64}()
    input = HLLInput(UL, UR, FL, FR, SL, SR, reject)

    println("kernel(...) out    = ", kernel(_solver, input).out)
    println("kernel(...) reject = ", reject.flag)

    reset!(reject)

    println("solver(...) out    = ", solver(_solver, input).out)
    println("solver(...) reject = ", reject.flag)
end

function test_hlle()
    println("\n================= HLLE test ==================")

    UL = [1.0, 2.0, 3.0]
    UR = [1.5, 2.5, 3.5]
    FL = [10.0, 20.0, 30.0]
    FR = [40.0, 50.0, 60.0]

    uL = 1.0
    cL = 2.0
    uR = -0.5
    cR = 1.5

    SL = min(uL - cL, uR - cR)
    SR = max(uL + cL, uR + cR)

    println("UL = ", UL)
    println("UR = ", UR)
    println("FL = ", FL)
    println("FR = ", FR)
    println("uL = ", uL, ", cL = ", cL)
    println("uR = ", uR, ", cR = ", cR)
    println("SL = ", SL)
    println("SR = ", SR)

    reject = Reject(false, "")
    reset!(reject)

    _solver = HLLESolver{Float64}()
    input = HLLEInput(UL, UR, FL, FR, SL, SR, reject)

    println("kernel(...) out    = ", kernel(_solver, input).out)
    println("kernel(...) reject = ", reject.flag)

    reset!(reject)

    println("solver(...) out    = ", solver(_solver, input).out)
    println("solver(...) reject = ", reject.flag)
end

function test_hllc()
    println("\n================= HLLC test ==================")

    UL = [1.0, 2.0, 3.0]
    UR = [1.5, 2.5, 3.5]
    FL = [10.0, 20.0, 30.0]
    FR = [40.0, 50.0, 60.0]

    rhoL = 1.0
    uL = 1.0
    pL = 1.0
    cL = 2.0

    rhoR = 0.8
    uR = -0.5
    pR = 0.7
    cR = 1.5

    SL = min(uL - cL, uR - cR)
    SR = max(uL + cL, uR + cR)

    SM = (pR - pL + rhoL * uL * (SL - uL) - rhoR * uR * (SR - uR)) /
         (rhoL * (SL - uL) - rhoR * (SR - uR))

    USL = [1.1, 2.1, 3.1]
    USR = [1.4, 2.4, 3.4]

    println("UL   = ", UL)
    println("UR   = ", UR)
    println("FL   = ", FL)
    println("FR   = ", FR)
    println("rhoL = ", rhoL, ", uL = ", uL, ", pL = ", pL, ", cL = ", cL)
    println("rhoR = ", rhoR, ", uR = ", uR, ", pR = ", pR, ", cR = ", cR)
    println("SL   = ", SL)
    println("SM   = ", SM)
    println("SR   = ", SR)
    println("USL  = ", USL)
    println("USR  = ", USR)

    reject = Reject(false, "")
    reset!(reject)

    _solver = HLLCSolver{Float64}()
    input = HLLCInput(UL, UR, FL, FR, USL, USR, SL, SM, SR, reject)

    println("kernel(...) out    = ", kernel(_solver, input).out)
    println("kernel(...) reject = ", reject.flag)

    reset!(reject)

    println("solver(...) out    = ", solver(_solver, input).out)
    println("solver(...) reject = ", reject.flag)
end

function test()
    test_rusanov()
    test_hll()
    test_hlle()
    test_hllc()
end

if abspath(PROGRAM_FILE) == @__FILE__
    test()
end
