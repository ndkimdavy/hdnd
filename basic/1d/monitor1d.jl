include("../../interface/imonitor.jl")
include("../mem.jl")

#===================================================================
                            Monitor 1D
===================================================================#

struct Monitor1D <: IMonitor end

function monitor!(monitor::Monitor1D,
    it::Int,
    t::T,
    dt::T,
    jxlocal::Int,
    nxlocal::Int,
    x::T,
    rank::Int,
    l::Int,
    vx::T,
    cs::T,
    p::T,
    Thj::T,
    rhoj::T,
    momj::T,
    Ehj::T,
    Q0::T,
    Qx::T,
    mfreejl::T,
    ceffjl::T,
    Erjl::T,
    Frjl::T,
    Prjl::T,
    Trjl::T,
    cS0::T,
    Sx::T,
    regime::Symbol) where {T<:AbstractFloat}

    fedd = Prjl / Erjl
    fred = abs(Frjl) / (ceffjl * Erjl)
    iregime = regime === :thick ? 1.0 : regime === :m1 ? 2.0 : regime === :thin ? 3.0 : error("Monitor1D invalid regime")

    lock(DBMEM.lock)

    try
        update! = (key::Symbol, value) -> begin
            if !haskey(DBMEM.data, key)
                DBMEM.data[key] = fill(NaN, nxlocal)
            end

            col = DBMEM.data[key]
            v = Float64(value)

            if !isequal(col[jxlocal], v)
                col[jxlocal] = v
            end
        end

        update!(:it, it)
        update!(:t, t)
        update!(:dt, dt)
        update!(:x, x)
        update!(:rank, rank)
        update!(:vx, vx)
        update!(:cs, cs)
        update!(:p, p)
        update!(:Th, Thj)
        update!(:rho, rhoj)
        update!(:mom, momj)
        update!(:Eh, Ehj)
        update!(:Q0, Q0)
        update!(:Qx, Qx)
        update!(Symbol("mfree", l), mfreejl)
        update!(Symbol("ceff", l), ceffjl)
        update!(Symbol("Er", l), Erjl)
        update!(Symbol("Frx", l), Frjl)
        update!(Symbol("Pr", l), Prjl)
        update!(Symbol("Tr", l), Trjl)
        update!(Symbol("cS0", l), cS0)
        update!(Symbol("Sx", l), Sx)
        update!(Symbol("fedd", l), fedd)
        update!(Symbol("fred", l), fred)
        update!(Symbol("regime", l), iregime)

    finally
        unlock(DBMEM.lock)
    end

    return
end
