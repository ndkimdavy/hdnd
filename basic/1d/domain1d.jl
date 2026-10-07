include("../../interface/idomain.jl")
include("../log.jl")

using MPI

#===================================================================
                            Domain 1D
===================================================================#

struct Domain1D <: IDomain
    nx::Int
    nxlocal::Int

    jxbegin::Int
    jxend::Int

    nprocx::Int
    coordx::Int

    rank::Int
    rankleft::Int
    rankright::Int

    comm::Union{MPI.Comm,Nothing}
end

function abort(domain::Domain1D, msg::String)

    log("ERROR: Domain1D")
    log("domain    = $(typeof(domain))")
    log("nx        = $(domain.nx)")
    log("nxlocal   = $(domain.nxlocal)")
    log("jxbegin   = $(domain.jxbegin)")
    log("jxend     = $(domain.jxend)")
    log("nprocx    = $(domain.nprocx)")
    log("coordx    = $(domain.coordx)")
    log("rank      = $(domain.rank)")
    log("rankleft  = $(domain.rankleft)")
    log("rankright = $(domain.rankright)")

    error("$(msg)")
end

function init(domain::Domain1D)::Domain1D

    nx = domain.nx
    nprocx = domain.nprocx

    if nx <= 0
        abort(domain, "Domain1D invalid nx=$(nx)")
    end

    if nprocx <= 0
        abort(domain, "Domain1D invalid nprocx=$(nprocx)")
    end

    if nprocx > nx
        abort(domain, "Domain1D invalid nprocx=$(nprocx), nx=$(nx)")
    end

    if !MPI.Initialized()
        provided = MPI.Init(threadlevel=:multiple)
    else
        provided = MPI.Query_thread()
    end

    if provided != MPI.THREAD_MULTIPLE
        abort(domain, "Domain1D requires MPI.THREAD_MULTIPLE")
    end

    comm = MPI.COMM_WORLD
    nproc = MPI.Comm_size(comm)

    if nprocx != nproc
        abort(domain, "Domain1D invalid nprocx=$(nprocx), nproc=$(nproc)")
    end

    commcart = MPI.Cart_create(
        comm,
        (nprocx,);
        periodic=(false,),
        reorder=false
    )

    rank = MPI.Comm_rank(commcart)
    coordx = Int(MPI.Cart_coords(commcart)[1])

    rankleft, rankright = MPI.Cart_shift(commcart, 0, 1)

    q, r = divrem(nx, nprocx)

    nxlocal = q + (coordx < r ? 1 : 0)

    jxbegin = coordx * q + min(coordx, r) + 1
    jxend = jxbegin + nxlocal - 1

    return Domain1D(
        nx,
        nxlocal,
        jxbegin,
        jxend,
        nprocx,
        coordx,
        rank,
        Int(rankleft),
        Int(rankright),
        commcart
    )
end

function hasleft(domain::Domain1D)::Bool
    return domain.rankleft != MPI.PROC_NULL
end

function hasright(domain::Domain1D)::Bool
    return domain.rankright != MPI.PROC_NULL
end

function exchange!(domain::Domain1D, sendleft::Array{T,1}, sendright::Array{T,1}, recvleft::Array{T,1}, recvright::Array{T,1}) where {T<:AbstractFloat}

    MPI.Sendrecv!(
        sendright,
        recvleft,
        domain.comm;
        dest=domain.rankright,
        source=domain.rankleft,
        sendtag=0,
        recvtag=0
    )

    MPI.Sendrecv!(
        sendleft,
        recvright,
        domain.comm;
        dest=domain.rankleft,
        source=domain.rankright,
        sendtag=1,
        recvtag=1
    )

    return nothing
end

function send(domain::Domain1D, msg::Tuple, dest::Int, tag::Int)
    MPI.send(msg, domain.comm; dest=dest, tag=tag)

    return nothing
end

function recv(domain::Domain1D, source::Int, tag::Int)::Tuple
    return MPI.recv(domain.comm; source=source, tag=tag)
end

function probe(domain::Domain1D, source::Int, tag::Int)::Bool
    return MPI.Iprobe(domain.comm; source=source, tag=tag)
end

function globalmin(domain::Domain1D, value::T)::T where {T<:AbstractFloat}
    return MPI.Allreduce(value, min, domain.comm)
end

function globalmax(domain::Domain1D, value::T)::T where {T<:AbstractFloat}
    return MPI.Allreduce(value, max, domain.comm)
end

function globalany(domain::Domain1D, value::Bool)::Bool
    return MPI.Allreduce(value ? 1 : 0, max, domain.comm) != 0
end
