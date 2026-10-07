include("../interface/iopacity.jl")
include("constant.jl")
include("log.jl")

#===================================================================
                            OpLosAlamos
===================================================================#

mutable struct OpLosAlamos{T<:AbstractFloat} <: IOpacity{T}
    fop::String

    Ti::Array{T,1}
    rhoj::Array{T,1}

    kR::Array{T,2}
    kP::Array{T,2}

    kRmg::Array{T,3}
    kPmg::Array{T,3}

    ratiokR::Array{T,1}
    ratiokP::Array{T,1}
end

function OpLosAlamos{T}(fop::String, ratiokR::Array{T,1}, ratiokP::Array{T,1}) where {T<:AbstractFloat}
    return OpLosAlamos{T}(
        fop,
        zeros(T, 0),
        zeros(T, 0),
        zeros(T, 0, 0),
        zeros(T, 0, 0),
        zeros(T, 0, 0, 0),
        zeros(T, 0, 0, 0),
        copy(ratiokR),
        copy(ratiokP)
    )
end

function OpLosAlamos{T}(fop::String) where {T<:AbstractFloat}
    return OpLosAlamos{T}(fop, zeros(T, 0), zeros(T, 0))
end

function findline(lines::Array{String,1}, key::String, iline::Int)::Int

    for i in iline:length(lines)
        if occursin(key, lines[i])
            return i
        end
    end

    return 0
end

function getline(::Type{T}, line::String)::Array{T,1} where {T<:AbstractFloat}

    val = Array{T,1}()
    regex = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][-+]?\d+)?"

    for m in eachmatch(regex, line)
        push!(val, parse(T, m.match))
    end

    return val
end

function interp2d(xi::Array{T,1}, yj::Array{T,1}, zij::AbstractArray{T,2}, x::T, y::T)::T where {T<:AbstractFloat}

    x = clamp(x, xi[1], xi[end])
    y = clamp(y, yj[1], yj[end])

    i = searchsortedlast(xi, x)
    j = searchsortedlast(yj, y)

    i = clamp(i, 1, length(xi) - 1)
    j = clamp(j, 1, length(yj) - 1)

    x1 = xi[i]
    x2 = xi[i+1]
    y1 = yj[j]
    y2 = yj[j+1]

    wx = (x - x1) / (x2 - x1)
    wy = (y - y1) / (y2 - y1)

    z11 = zij[i, j]
    z21 = zij[i+1, j]
    z12 = zij[i, j+1]
    z22 = zij[i+1, j+1]

    return (T(1) - wx) * (T(1) - wy) * z11 +
           wx * (T(1) - wy) * z21 +
           (T(1) - wx) * wy * z12 +
           wx * wy * z22
end

function abort(opacity::OpLosAlamos{T}, msg::String, iline::Int, key::String, val::Array{T,1}, iT::Int, irho::Int, ig::Int) where {T<:AbstractFloat}
    log("ERROR: OpLosAlamos")
    log("opacity = $(typeof(opacity))")
    log("fop     = $(opacity.fop)")
    log("iline   = $(iline)")
    log("key     = $(key)")
    log("val     = $(val)")
    log("iT      = $(iT)")
    log("irho    = $(irho)")
    log("ig      = $(ig)")

    error("$(msg)")
end

function load!(opacity::OpLosAlamos{T})::OpLosAlamos{T} where {T<:AbstractFloat}

    if !isfile(opacity.fop)
        abort(opacity, "OpLosAlamos file not found", 0, opacity.fop, zeros(T, 0), 0, 0, 0)
    end

    lines = readlines(opacity.fop)

    key = "Temperature grid used the following"
    iline = findline(lines, key, 1)

    if iline == 0
        abort(opacity, "OpLosAlamos line not found", iline, key, zeros(T, 0), 0, 0, 0)
    end

    val = getline(T, lines[iline])

    if isempty(val)
        abort(opacity, "OpLosAlamos invalid temperature header", iline, key, val, 0, 0, 0)
    end

    nTi = Int(val[end])

    Ti = zeros(T, 0)
    iline += 1

    while length(Ti) < nTi
        append!(Ti, getline(T, lines[iline]))
        iline += 1
    end

    key = "Density grid used the following"
    iline = findline(lines, key, iline)

    if iline == 0
        abort(opacity, "OpLosAlamos line not found", iline, key, zeros(T, 0), 0, 0, 0)
    end

    val = getline(T, lines[iline])

    if isempty(val)
        abort(opacity, "OpLosAlamos invalid density header", iline, key, val, 0, 0, 0)
    end

    nrhoj = Int(val[end])

    rhoj = zeros(T, 0)
    iline += 1

    while length(rhoj) < nrhoj
        append!(rhoj, getline(T, lines[iline]))
        iline += 1
    end

    key = "Photon grid used the following"
    iline = findline(lines, key, iline)

    if iline == 0
        abort(opacity, "OpLosAlamos line not found", iline, key, zeros(T, 0), 0, 0, 0)
    end

    val = getline(T, lines[iline])

    if isempty(val)
        abort(opacity, "OpLosAlamos invalid photon header", iline, key, val, 0, 0, 0)
    end

    ng = Int(val[end])

    gEi = zeros(T, 0)
    iline += 1

    while length(gEi) < ng
        append!(gEi, getline(T, lines[iline]))
        iline += 1
    end

    if isempty(opacity.ratiokR)
        opacity.ratiokR = ones(T, ng)
    end

    if isempty(opacity.ratiokP)
        opacity.ratiokP = ones(T, ng)
    end

    if length(opacity.ratiokR) != ng
        abort(opacity, "OpLosAlamos inconsistent ratiokR size", 0, "ratiokR", opacity.ratiokR, 0, 0, 0)
    end

    if length(opacity.ratiokP) != ng
        abort(opacity, "OpLosAlamos inconsistent ratiokP size", 0, "ratiokP", opacity.ratiokP, 0, 0, 0)
    end

    if any(x -> !isfinite(x) || x <= T(0), opacity.ratiokR)
        abort(opacity, "OpLosAlamos invalid ratiokR", 0, "ratiokR", opacity.ratiokR, 0, 0, 0)
    end

    if any(x -> !isfinite(x) || x <= T(0), opacity.ratiokP)
        abort(opacity, "OpLosAlamos invalid ratiokP", 0, "ratiokP", opacity.ratiokP, 0, 0, 0)
    end

    kR = zeros(T, nTi, nrhoj)
    kP = zeros(T, nTi, nrhoj)
    kRmg = zeros(T, nTi, nrhoj, ng)
    kPmg = zeros(T, nTi, nrhoj, ng)

    key = "Rosseland and Planck opacities and free electrons"

    for iT in 1:nTi
        iline = findline(lines, key, iline)

        if iline == 0
            abort(opacity, "OpLosAlamos line not found", iline, key, zeros(T, 0), iT, 0, 0)
        end

        iline += 2

        for irho in 1:nrhoj
            val = getline(T, lines[iline])

            if length(val) < 3
                abort(opacity, "OpLosAlamos invalid gray opacity line", iline, key, val, iT, irho, 0)
            end

            kR[iT, irho] = val[2] * T(cm2g_to_m2kg)
            kP[iT, irho] = val[3] * T(cm2g_to_m2kg)

            iline += 1
        end
    end

    key = "Multigroup opacities"
    iline = findline(lines, key, iline)

    if iline == 0
        abort(opacity, "OpLosAlamos line not found", iline, key, zeros(T, 0), 0, 0, 0)
    end

    iline += 1

    key = "for T, density ="

    for iT in 1:nTi
        for irho in 1:nrhoj
            iline = findline(lines, key, iline)

            if iline == 0
                abort(opacity, "OpLosAlamos line not found", iline, key, zeros(T, 0), iT, irho, 0)
            end

            iline += 1

            for ig in 1:ng
                val = getline(T, lines[iline])

                if length(val) < 3
                    abort(opacity, "OpLosAlamos invalid multigroup opacity line", iline, key, val, iT, irho, ig)
                end

                kRmg[iT, irho, ig] = val[2] * T(cm2g_to_m2kg)
                kPmg[iT, irho, ig] = val[3] * T(cm2g_to_m2kg)

                iline += 1
            end
        end
    end

    opacity.Ti = Ti .* T(keV_to_K)
    opacity.rhoj = rhoj .* T(gcm3_to_kgm3)
    opacity.kR = kR
    opacity.kP = kP
    opacity.kRmg = kRmg
    opacity.kPmg = kPmg

    return opacity
end

function ngroup(opacity::OpLosAlamos{T})::Int where {T<:AbstractFloat}
    return size(opacity.kRmg, 3)
end

function kappaR(opacity::OpLosAlamos{T}, Th::T, rho::T, ig::Int)::T where {T<:AbstractFloat}

    if ig == -1
        return interp2d(opacity.Ti, opacity.rhoj, opacity.kR, Th, rho)
    end

    if ig < 1 || ig > ngroup(opacity)
        abort(opacity, "OpLosAlamos invalid group index", 0, "", zeros(T, 0), 0, 0, ig)
    end

    return opacity.ratiokR[ig] * interp2d(opacity.Ti, opacity.rhoj, @view(opacity.kRmg[:, :, ig]), Th, rho)
end

function kappaP(opacity::OpLosAlamos{T}, Th::T, rho::T, ig::Int)::T where {T<:AbstractFloat}

    if ig == -1
        return interp2d(opacity.Ti, opacity.rhoj, opacity.kP, Th, rho)
    end

    if ig < 1 || ig > ngroup(opacity)
        abort(opacity, "OpLosAlamos invalid group index", 0, "", zeros(T, 0), 0, 0, ig)
    end

    return opacity.ratiokP[ig] * interp2d(opacity.Ti, opacity.rhoj, @view(opacity.kPmg[:, :, ig]), Th, rho)
end

#===================================================================
                              Plot
===================================================================#

if Sys.islinux() && !haskey(ENV, "QT_QPA_PLATFORM")
    ENV["QT_QPA_PLATFORM"] = "xcb"
end

using Plots
using Plots.PlotMeasures

function plotop(opacity::OpLosAlamos{T}) where {T<:AbstractFloat}
    isempty(opacity.Ti) && load!(opacity)

    X = log10.(opacity.rhoj)
    Y = log10.(opacity.Ti)

    KP = [[log10(max(kappaP(opacity, opacity.Ti[iT], opacity.rhoj[irho], l), eps(T))) for iT in eachindex(opacity.Ti), irho in eachindex(opacity.rhoj)] for l in 1:ngroup(opacity)]
    KR = [[log10(max(kappaR(opacity, opacity.Ti[iT], opacity.rhoj[irho], l), eps(T))) for iT in eachindex(opacity.Ti), irho in eachindex(opacity.rhoj)] for l in 1:ngroup(opacity)]

    zlimP = (minimum(minimum.(KP)), maximum(maximum.(KP)))
    zlimR = (minimum(minimum.(KR)), maximum(maximum.(KR)))

    graph = []

    pal = :inferno

    for l in 1:ngroup(opacity)
        push!(graph, heatmap(
            X,
            Y,
            KP[l],
            xlabel="",
            ylabel=l == 1 ? "log10(T [K])" : "",
            title="kP [$(l)]",
            clims=zlimP,
            color=pal,
            colorbar=false,
            aspect_ratio=:auto,
            titlefontsize=10,
            guidefontsize=10,
            tickfontsize=10,
            bottom_margin=10mm,
            left_margin=10mm
        ))
    end

    push!(graph, heatmap(
        X,
        Y,
        KP[ngroup(opacity)],
        xlabel="",
        ylabel="",
        title="",
        clims=zlimP,
        color=pal,
        colorbar=true,
        framestyle=:none,
        aspect_ratio=:auto,
        ticks=false
    ))

    for l in 1:ngroup(opacity)
        push!(graph, heatmap(
            X,
            Y,
            KR[l],
            xlabel=l == 3 ? "log10(rho [kg/m^3])" : "",
            ylabel=l == 1 ? "log10(T [K])" : "",
            title="kR [$(l)]",
            clims=zlimR,
            color=pal,
            colorbar=false,
            aspect_ratio=:auto,
            titlefontsize=10,
            guidefontsize=10,
            tickfontsize=10,
            bottom_margin=10mm,
            left_margin=10mm
        ))
    end

    push!(graph, heatmap(
        X,
        Y,
        KR[ngroup(opacity)],
        xlabel="",
        ylabel="",
        title="",
        clims=zlimR,
        color=pal,
        colorbar=true,
        framestyle=:none,
        aspect_ratio=:auto,
        ticks=false
    ))

    fig = Plots.plot(
        graph...,
        layout=(2, ngroup(opacity) + 1),
        size=(1900, 720),
        plot_title="Los Alamos Opacity Tables",
        plot_titlefontsize=12
    )

    display(fig)
    println("Press Enter to close ...")
    readline()
end

function main()
    fop = "input/opAr4g.txt"
    ratiokR = [1.0, 1.0, 1.0, 1.0]
    ratiokP = [1.0, 1.0, 1.0, 1.0]

    opacity = OpLosAlamos{Float64}(fop, ratiokR, ratiokP)
    load!(opacity)
    plotop(opacity)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
