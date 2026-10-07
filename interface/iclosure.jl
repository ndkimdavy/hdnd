include("isolver.jl")

abstract type IClosureInput <: ISolverInput end
abstract type IClosureOutput <: ISolverOutput end
abstract type IClosure{Tin<:IClosureInput,Tout<:IClosureOutput} <: ISolver{Tin,Tout} end

function abort(closure::IClosure{Tin,Tout}, msg::String, args...) where {Tin<:IClosureInput,Tout<:IClosureOutput}
    error("closure=$(typeof(closure)), msg=$(msg), args=$(args)")
end

function closure!(closure::IClosure{Tin,Tout}, input::Tin)::Tout where {Tin<:IClosureInput,Tout<:IClosureOutput}
    error("closure! function not implemented for closure=$(typeof(closure)), input=$(Tin), output=$(Tout)")
end