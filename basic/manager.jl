include("../interface/iio.jl")
include("reject.jl")
include("oplosalamos.jl")
include("1d/coupler1d.jl")

#===================================================================
                            Manager
===================================================================#

struct ManagerConfig{T}
    ndim::Int
    nx::Int
    ny::Int
    nz::Int
    dx::T
    dy::T
    dz::T
    nprocx::Int
    nprocy::Int
    nprocz::Int
    tsim::T
    courant::T
    mode::String
    gamma::T
    mu::T
    tolsim::T
    gE::Array{T,1}
    srcop::String
    fop::String
    ratiokR::Array{T,1}
    ratiokP::Array{T,1}
    nspectral::Int
    nangular::Int
    rho0::T
    vx0::T
    vy0::T
    vz0::T
    Th0::T
    Er0::Array{T,1}
    Frx0::Array{T,1}
    Fry0::Array{T,1}
    Frz0::Array{T,1}
    hrho0::Array{T,1}
    hvx0::Array{T,1}
    hvy0::Array{T,1}
    hvz0::Array{T,1}
    hTh0::Array{T,1}
    hEr0::Array{T,1}
    hFrx0::Array{T,1}
    hFry0::Array{T,1}
    hFrz0::Array{T,1}
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
    outpath::String
    outfile::String
    ioformat::String
    fps::T
end

struct Manager{T}
    config::ManagerConfig{T}
    opacity::IOpacity{T}
end

function init(config::ManagerConfig{T}) where {T<:AbstractFloat}

    if config.srcop == "losalamos"
        opacity = OpLosAlamos{T}(config.fop, config.ratiokR, config.ratiokP)
        load!(opacity)
    else
        error("Manager invalid opacity source: $(config.srcop)")
    end

    return Manager{T}(
        config,
        opacity
    )
end

function run(manager::Manager{T}, io::Union{IIO,Nothing}) where {T<:AbstractFloat}

    config = manager.config
    opacity = manager.opacity

    reject = Reject(false, "")
    reset!(reject)

    if config.ndim == 1
        domain = Domain1D(config.nx, 0, 0, 0, config.nprocx, 0, 0, 0, 0, nothing)
        domain = init(domain)

        if domain.rank == 0
            log("\n================ App config ================")
            log("ndim      = $(config.ndim)")
            log("nx        = $(config.nx)")
            log("ny        = $(config.ny)")
            log("nz        = $(config.nz)")
            log("dx        = $(config.dx)")
            log("dy        = $(config.dy)")
            log("dz        = $(config.dz)")
            log("nprocx    = $(config.nprocx)")
            log("nprocy    = $(config.nprocy)")
            log("nprocz    = $(config.nprocz)")
            log("tsim      = $(config.tsim)")
            log("courant   = $(config.courant)")
            log("mode      = $(config.mode)")
            log("gamma     = $(config.gamma)")
            log("mu        = $(config.mu)")
            log("tolsim    = $(config.tolsim)")
            log("gE        = $(config.gE)")
            log("srcop     = $(config.srcop)")
            log("fop       = $(config.fop)")
            log("ratiokR   = $(config.ratiokR)")
            log("ratiokP   = $(config.ratiokP)")
            log("nspectral = $(config.nspectral)")
            log("nangular  = $(config.nangular)")
            log("rho0      = $(config.rho0)")
            log("vx0       = $(config.vx0)")
            log("vy0       = $(config.vy0)")
            log("vz0       = $(config.vz0)")
            log("Th0       = $(config.Th0)")
            log("Er0       = $(config.Er0)")
            log("Frx0      = $(config.Frx0)")
            log("Fry0      = $(config.Fry0)")
            log("Frz0      = $(config.Frz0)")
            log("hrho0     = $(config.hrho0)")
            log("hvx0      = $(config.hvx0)")
            log("hvy0      = $(config.hvy0)")
            log("hvz0      = $(config.hvz0)")
            log("hTh0      = $(config.hTh0)")
            log("hEr0      = $(config.hEr0)")
            log("hFrx0     = $(config.hFrx0)")
            log("hFry0     = $(config.hFry0)")
            log("hFrz0     = $(config.hFrz0)")
            log("wall      = $(config.wall)")
            log("gmax      = $(config.gmax)")
            log("mcool     = $(config.mcool)")
            log("hpcool    = $(config.hpcool)")
            log("mclosure  = $(config.mclosure)")
            log("hpclosure = $(config.hpclosure)")
            log("fclosure  = $(config.fclosure)")
            log("hprmodel  = $(config.hprmodel)")
            log("hsolver   = $(config.hsolver)")
            log("rsolver   = $(config.rsolver)")
            log("retry     = $(config.retry)")
            log("outpath   = $(config.outpath)")
            log("outfile   = $(config.outfile)")
            log("ioformat  = $(config.ioformat)")
            log("fps       = $(config.fps)")
        end

        input = Coupler1DInput(
            config.dx,
            domain,
            config.tsim,
            config.courant,
            config.mode,
            config.gamma,
            config.mu,
            config.tolsim,
            config.gE,
            opacity,
            config.nspectral,
            config.nangular,
            config.rho0,
            config.vx0,
            config.Th0,
            config.Er0,
            config.Frx0,
            config.hrho0,
            config.hvx0,
            config.hTh0,
            config.hEr0,
            config.hFrx0,
            config.wall,
            config.gmax,
            config.mcool,
            config.hpcool,
            config.mclosure,
            config.hpclosure,
            config.fclosure,
            config.hprmodel,
            config.hsolver,
            config.rsolver,
            config.retry,
            reject
        )

        solver = Coupler1DSolver{Coupler1DInput{T},Coupler1DOutput{T}}()

    else
        error("Manager currently supports ndim = 1 only")
    end

    if domain.rank == 0 && lowercase(strip(config.outpath)) != "none"
        rm(config.outpath; force=true, recursive=true)
        mkpath(config.outpath)
    end

    dbmem = DBMem{T}(OrderedDict{Symbol,Array{T,1}}(), nothing)

    task = Threads.@spawn dispatch(solver, input)

    lastit = -1

    try
        while !istaskdone(task)
            lastit = sync(manager, io, lastit, domain, dbmem)
            sleep(T(1) / config.fps)
        end

        wait(task)

        if reject.flag
            error("Manager Coupler$(config.ndim)D rejected")
        end

        lastit = sync(manager, io, lastit, domain, dbmem, true)

    finally
        if io !== nothing
            release(io)
        end
    end

    return nothing
end

function sync(manager::Manager{T}, io::Union{IIO,Nothing}, lastit::Int, domain::Domain1D, dbmem::DBMem{T}, final::Bool=false) where {T<:AbstractFloat}

    dbmemlocal = DBMem{T}(OrderedDict{Symbol,Array{T,1}}(), nothing)

    lock(DBMEM.lock)

    try
        for (key, data) in DBMEM.data
            dbmemlocal.data[key] = copy(data)
        end

    finally
        unlock(DBMEM.lock)
    end

    copydb! = (data, jxbegin::Int, jxend::Int) -> begin
        for (key, values) in data
            if !haskey(dbmem.data, key)
                dbmem.data[key] = fill(T(NaN), domain.nx)
            end

            dbmem.data[key][jxbegin:jxend] .= values
        end
    end

    if domain.rank != 0
        if !haskey(dbmemlocal.data, :it)
            if final
                send(domain, (domain.jxbegin, domain.jxend, dbmemlocal.data), 0, 3)
            end

            return lastit
        end

        it = Int(dbmemlocal.data[:it][1])

        if !final && it == lastit
            return lastit
        end

        send(domain, (domain.jxbegin, domain.jxend, dbmemlocal.data), 0, final ? 3 : 2)

        return it
    end

    if haskey(dbmemlocal.data, :it)
        copydb!(dbmemlocal.data, domain.jxbegin, domain.jxend)
    end

    if final
        for irank in 1:(domain.nprocx-1)
            msg = recv(domain, irank, 3)

            jxbegin, jxend, data = msg
            copydb!(data, jxbegin, jxend)
        end

    else
        for irank in 1:(domain.nprocx-1)
            while probe(domain, irank, 2)
                msg = recv(domain, irank, 2)

                jxbegin, jxend, data = msg
                copydb!(data, jxbegin, jxend)
            end
        end
    end

    if !haskey(dbmemlocal.data, :it)
        return lastit
    end

    it = Int(dbmemlocal.data[:it][1])

    if !final && it == lastit
        return lastit
    end

    if !haskey(dbmem.data, :rank) || any(x -> !isfinite(x), dbmem.data[:rank])
        return lastit
    end

    if io !== nothing
        write(io, dbmem)
    end

    return it
end
