include("../basic/constant.jl")
include("../basic/log.jl")
include("../basic/1d/closure1d.jl")

using Printf
using Base.Threads

#===================================================================
              Generate HDND Closure1D natural dataset
===================================================================#

#===================================================================
                               Inputs
===================================================================#

# Numerical type
T = Float64

# Domain limits
domfred = Float64[0.0, 0.999]
domta = Float64[1.0e-12, 1.0e12]
domfa = Float64[1.0e-12, 0.999]
domr = Float64[1.0e-12, 1.0]

# Sampling
nfred = 128
nta = 128
nfa = 64
nr = 64

# Reference group
eref = 1.0

# Tolerance
tolsim = 1.0e-2

# Numerical quadrature
nspectral = 64
nangular = 64

# Closure
mclosure = "randomsearch"
hpclosure = Float64[1.0e-12, 2.5e-10, -0.99, 0.99, 0.5, 1.0e-9, 10, 5, 100, 42]
fclosure = "none"

# Solver control
retry = 0

# Output
fcsv = "genclosure1d_128x128x64x64.csv"
sep = ","
pd = "%24.15e"

#===================================================================
                              Functions
===================================================================#

function generate()
    vlogr = collect(range(log10(domr[1]), log10(domr[2]), length=nr))
    vfred = collect(range(domfred[1], domfred[2], length=nfred))
    vlogta = collect(range(log10(domta[1]), log10(domta[2]), length=nta))
    vlogfa = collect(range(log10(domfa[1]), log10(domfa[2]), length=nfa))

    PD = Printf.Format(pd)

    npoint = nr * nfred * nta * nfa
    nprint = max(1, div(npoint, 1000))

    nthd = Threads.nthreads()
    blockFred = [(div((ithd-1)*nfred, nthd)+1):div(ithd*nfred, nthd) for ithd in 1:nthd]
    rejectFred = [Reject(false, "") for ithd in 1:nthd]
    bufFred = [IOBuffer() for ithd in 1:nthd]

    npointFred = zeros(Int, nthd)
    nacceptFred = zeros(Int, nthd)
    nrejectFred = zeros(Int, nthd)

    lockProgress = ReentrantLock()

    t0 = clock_ns()

    Threads.@threads :static for ithd in 1:nthd
        _reject = rejectFred[ithd]

        solver = Closure1D{Closure1DInput{T},Closure1DOutput{T}}()

        io = bufFred[ithd]

        for ifred in blockFred[ithd]
            fred = vfred[ifred]

            for logta in vlogta
                ta = T(10)^T(logta)

                for logfa in vlogfa
                    fa = T(10)^T(logfa)

                    for logr in vlogr
                        r = T(10)^T(logr)
                        ceffjl = r * T(c)

                        gE = T[
                            fa*T(eref),
                            T(eref)
                        ]

                        l = 1

                        E1 = gE[l] * T(keV_to_J)
                        E2 = gE[l+1] * T(keV_to_J)

                        nu1 = E1 / T(h)
                        nu2 = E2 / T(h)

                        Erjl = T(aR) * (ta * sqrt(nu1 * nu2) * T(h) / T(kB))^4
                        Frjl = fred * ceffjl * Erjl

                        reset!(_reject)

                        input = Closure1DInput(
                            l,
                            T(tolsim),
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
                            _reject
                        )

                        output = closure!(solver, input)

                        if _reject.flag
                            nrejectFred[ithd] += 1

                        else
                            fedd = output.fedd

                            if !isfinite(fedd)
                                nrejectFred[ithd] += 1

                            elseif fedd < T(1) / T(3) - T(tolsim) || fedd > T(1) + T(tolsim)
                                nrejectFred[ithd] += 1

                            else
                                nacceptFred[ithd] += 1

                                println(io, join((
                                        Printf.format(PD, fred),
                                        Printf.format(PD, logta),
                                        Printf.format(PD, logfa),
                                        Printf.format(PD, logr),
                                        Printf.format(PD, fedd)
                                    ), sep))
                            end
                        end

                        npointFred[ithd] += 1

                        if npointFred[ithd] % nprint == 0
                            lock(lockProgress)

                            elapsed = Float64(clock_ns() - t0) * 1.0e-9
                            progress("GENCLOSURE1D", Float64(sum(npointFred)) / Float64(npoint), elapsed)

                            unlock(lockProgress)
                        end
                    end
                end
            end
        end

        lock(lockProgress)

        elapsed = Float64(clock_ns() - t0) * 1.0e-9
        progress("GENCLOSURE1D", Float64(sum(npointFred)) / Float64(npoint), elapsed)

        unlock(lockProgress)
    end

    rm(fcsv; force=true)

    open(fcsv, "w") do io
        println(io, join(("fred", "logta", "logfa", "logr", "fedd"), sep))

        for ithd in 1:nthd
            write(io, take!(bufFred[ithd]))
        end
    end

    naccept = sum(nacceptFred)
    nreject = sum(nrejectFred)

    log("")
    log("file    = $(fcsv)")
    log("points  = $(npoint)")
    log("accept  = $(naccept)")
    log("reject  = $(nreject)")

    return nothing
end

generate()
