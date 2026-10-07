abstract type IOpacity{T<:AbstractFloat} end

function abort(opacity::IOpacity{T}, msg::String, args...) where {T<:AbstractFloat}
    error("opacity=$(typeof(opacity)), msg=$(msg), args=$(args)")
end

function load!(opacity::IOpacity{T})::IOpacity{T} where {T<:AbstractFloat}
    error("load! function not implemented for opacity=$(typeof(opacity))")
end

function ngroup(opacity::IOpacity{T})::Int where {T<:AbstractFloat}
    error("ngroup function not implemented for opacity=$(typeof(opacity))")
end

function kappaR(opacity::IOpacity{T}, Th::T, rho::T, ig::Int)::T where {T<:AbstractFloat}
    error("kappaR function not implemented for opacity=$(typeof(opacity)), Th=$(T), rho=$(T), ig=$(Int)")
end

function kappaP(opacity::IOpacity{T}, Th::T, rho::T, ig::Int)::T where {T<:AbstractFloat}
    error("kappaP function not implemented for opacity=$(typeof(opacity)), Th=$(T), rho=$(T), ig=$(Int)")
end
