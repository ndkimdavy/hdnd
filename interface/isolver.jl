abstract type ISolverInput end
abstract type ISolverOutput end
abstract type ISolver{Tin<:ISolverInput,Tout<:ISolverOutput} end

function abort(solver::ISolver{Tin,Tout}, msg::String, args...) where {Tin<:ISolverInput,Tout<:ISolverOutput}
    error("solver=$(typeof(solver)), msg=$(msg), args=$(args)")
end

function kernel(solver::ISolver{Tin,Tout}, input::Tin)::Tout where {Tin<:ISolverInput,Tout<:ISolverOutput}
    error("kernel function not implemented for solver=$(typeof(solver)), input=$(Tin), output=$(Tout)")
end

function solver(solver::ISolver{Tin,Tout}, input::Tin)::Tout where {Tin<:ISolverInput,Tout<:ISolverOutput}
    return kernel(solver, input)
end

function corrector(solver::ISolver{Tin,Tout}, x::T, xmin::T, xmax::T, tol::T, force::Bool=false) where {Tin<:ISolverInput,Tout<:ISolverOutput,T<:AbstractFloat}

    if !isfinite(x)
        return x
    end

    if force
        return clamp(x, xmin, xmax)
    end

    if abs(x - xmin) <= tol || abs(x - xmax) <= tol
        return clamp(x, xmin, xmax)
    end

    return x
end
