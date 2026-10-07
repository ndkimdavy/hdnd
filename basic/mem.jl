using OrderedCollections

#===================================================================
                              DBMem
                          Shared Data Base
===================================================================#

struct DBMem{T}
    data::OrderedDict{Symbol,Array{T,1}}
    lock::Union{ReentrantLock,Nothing}
end

#===================================================================
                              CMem
                        Shared Cache Memory
===================================================================#

struct CMem
    data::OrderedDict{Any,Any}
    lock::Union{ReentrantLock,Nothing}
end

#===================================================================
                            Global DB
===================================================================#

const DBMEM = DBMem{Float64}(
    OrderedDict{Symbol,Array{Float64,1}}(),
    ReentrantLock()
)

#===================================================================
                           Global Cache
===================================================================#

const CMEM = CMem(
    OrderedDict{Any,Any}(),
    ReentrantLock()
)