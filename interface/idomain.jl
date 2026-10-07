abstract type IDomain end

function abort(domain::IDomain, msg::String, args...)
    error("domain=$(typeof(domain)), msg=$(msg), args=$(args)")
end

function init(domain::IDomain, args...)::IDomain
    error("init function not implemented for domain=$(typeof(domain)), args=$(typeof(args))")
end

function hasleft(domain::IDomain)::Bool
    error("hasleft function not implemented for domain=$(typeof(domain))")
end

function hasright(domain::IDomain)::Bool
    error("hasright function not implemented for domain=$(typeof(domain))")
end

function hasbottom(domain::IDomain)::Bool
    error("hasbottom function not implemented for domain=$(typeof(domain))")
end

function hastop(domain::IDomain)::Bool
    error("hastop function not implemented for domain=$(typeof(domain))")
end

function hasback(domain::IDomain)::Bool
    error("hasback function not implemented for domain=$(typeof(domain))")
end

function hasfront(domain::IDomain)::Bool
    error("hasfront function not implemented for domain=$(typeof(domain))")
end

function exchange!(domain::IDomain, args...)
    error("exchange! function not implemented for domain=$(typeof(domain)), args=$(typeof(args))")
end

function send(domain::IDomain, msg::Tuple, dest::Int, tag::Int)
    error("send function not implemented for domain=$(typeof(domain)), msg=$(typeof(msg)), dest=$(Int), tag=$(Int)")
end

function recv(domain::IDomain, source::Int, tag::Int)::Tuple
    error("recv function not implemented for domain=$(typeof(domain)), source=$(Int), tag=$(Int)")
end

function probe(domain::IDomain, source::Int, tag::Int)::Bool
    error("probe function not implemented for domain=$(typeof(domain)), source=$(Int), tag=$(Int)")
end

function globalmin(domain::IDomain, value::T)::T where {T<:AbstractFloat}
    error("globalmin function not implemented for domain=$(typeof(domain)), value=$(T)")
end

function globalmax(domain::IDomain, value::T)::T where {T<:AbstractFloat}
    error("globalmax function not implemented for domain=$(typeof(domain)), value=$(T)")
end

function globalany(domain::IDomain, value::Bool)::Bool
    error("globalany function not implemented for domain=$(typeof(domain)), value=$(Bool)")
end
