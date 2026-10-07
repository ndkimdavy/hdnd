include("../basic/manager.jl")
include("../basic/1d/ioparaview1d.jl")

#===================================================================
                                App
===================================================================#

const REAL = let real = lowercase(get(ENV, "REAL", "float64"))
    if real == "float64" || real == "f64"
        Float64
    elseif real == "float32" || real == "f32"
        Float32
    else
        error("Invalid REAL=$(real). Use REAL=float64 or REAL=float32")
    end
end

function readconfig(file::String, ::Type{T}) where {T<:AbstractFloat}
    config = Dict{String,Any}()
    lines = readlines(file)
    iline = 1

    while iline <= length(lines)
        line = strip(lines[iline])

        if isempty(line) || startswith(line, "#") || !occursin("=", line)
            iline += 1
            continue
        end

        key, val = strip.(split(line, "=", limit=2))

        if startswith(val, "[")
            while !endswith(val, "]")
                iline += 1

                if iline > length(lines)
                    error("Invalid array for $(key) in $(file)")
                end

                nextline = strip(lines[iline])

                if !isempty(nextline) && !startswith(nextline, "#")
                    val *= " " * nextline
                end
            end

            inner = strip(val[2:(end-1)])
            config[key] = isempty(inner) ? T[] : parse.(T, strip.(split(inner, ",")))

        elseif key in ("ndim", "nx", "ny", "nz", "nprocx", "nprocy", "nprocz", "nspectral", "nangular", "retry")
            config[key] = parse(Int, val)

        elseif key in ("mode", "srcop", "fop", "wall", "mcool", "mclosure", "fclosure", "hsolver", "rsolver", "outpath", "outfile", "ioformat")
            config[key] = val

        else
            config[key] = parse(T, val)
        end

        iline += 1
    end

    return config
end

function app()
    T = REAL

    file = "-f" in ARGS ? ARGS[findfirst(ARGS .== "-f")+1] : "input/config1d.txt"

    config = readconfig(file, T)

    ndim = Int(get(config, "ndim", 1))
    nx = Int(get(config, "nx", 128))
    ny = Int(get(config, "ny", 1))
    nz = Int(get(config, "nz", 1))

    dx = T(get(config, "dx", T(7.8125e-6)))
    dy = T(get(config, "dy", T(1.0)))
    dz = T(get(config, "dz", T(1.0)))

    nprocx = Int(get(config, "nprocx", 1))
    nprocy = Int(get(config, "nprocy", 1))
    nprocz = Int(get(config, "nprocz", 1))

    tsim = T(get(config, "tsim", T(1.0e-12)))
    courant = T(get(config, "courant", T(0.5)))

    mode = String(get(config, "mode", "meso"))
    gamma = T(get(config, "gamma", T(1.6666666666666667)))
    mu = T(get(config, "mu", T(39.948)))
    tolsim = T(get(config, "tolsim", T(1.0e-3)))

    gE = T.(get(config, "gE", T[1.47e-3, 7.94e-3, 6.11e-2, 5.0e-1, 16.0]))
    ng = length(gE) - 1

    srcop = String(get(config, "srcop", "losalamos"))
    fop = String(get(config, "fop", "input/opAr4g.txt"))
    ratiokR = T.(get(config, "ratiokR", ones(T, ng)))
    ratiokP = T.(get(config, "ratiokP", ones(T, ng)))

    nspectral = Int(get(config, "nspectral", 1024))
    nangular = Int(get(config, "nangular", 1024))

    rho0 = T(get(config, "rho0", T(1.0)))
    vx0 = T(get(config, "vx0", T(0.0)))
    vy0 = T(get(config, "vy0", T(0.0)))
    vz0 = T(get(config, "vz0", T(0.0)))
    Th0 = T(get(config, "Th0", T(1.160452e4)))

    Er0 = T.(get(config, "Er0", ones(T, ng)))
    Frx0 = T.(get(config, "Frx0", zeros(T, ng)))
    Fry0 = T.(get(config, "Fry0", zeros(T, ng)))
    Frz0 = T.(get(config, "Frz0", zeros(T, ng)))

    hrho0 = T.(get(config, "hrho0", T[0.5, 0.0]))
    hvx0 = T.(get(config, "hvx0", T[0.5, 0.0]))
    hvy0 = T.(get(config, "hvy0", T[0.5, 0.0]))
    hvz0 = T.(get(config, "hvz0", T[0.5, 0.0]))
    hTh0 = T.(get(config, "hTh0", T[0.5, 0.0]))

    hEr0 = T.(get(config, "hEr0", T[0.5, 0.0]))
    hFrx0 = T.(get(config, "hFrx0", T[0.5, 0.0]))
    hFry0 = T.(get(config, "hFry0", T[0.5, 0.0]))
    hFrz0 = T.(get(config, "hFrz0", T[0.5, 0.0]))

    wall = String(get(config, "wall", "none"))

    gmax = T(get(config, "gmax", T(1.0e3)))

    mcool = String(get(config, "mcool", "none"))
    hpcool = T.(get(config, "hpcool", zeros(T, 6 * ng)))

    mclosure = String(get(config, "mclosure", "randomsearch"))
    hpclosure = T.(get(config, "hpclosure", T[1.0e-12, 2.5e-10, -0.99, 0.99, 0.5, 1.0e-12, 10.0, 5.0, 1000.0, 42.0]))
    fclosure = String(get(config, "fclosure", "none"))

    hprmodel = T.(get(config, "hprmodel", T[1.0e-6, 1.0e-3]))

    hsolver = String(get(config, "hsolver", "hllc"))
    rsolver = String(get(config, "rsolver", "hll"))
    retry = Int(get(config, "retry", 5))

    outpath = String(get(config, "outpath", "none"))
    outfile = String(get(config, "outfile", "none"))
    ioformat = String(get(config, "ioformat", "none"))
    fps = T(get(config, "fps", T(1.0)))

    managerconfig = ManagerConfig(
        ndim,
        nx,
        ny,
        nz,
        dx,
        dy,
        dz,
        nprocx,
        nprocy,
        nprocz,
        tsim,
        courant,
        mode,
        gamma,
        mu,
        tolsim,
        gE,
        srcop,
        fop,
        ratiokR,
        ratiokP,
        nspectral,
        nangular,
        rho0,
        vx0,
        vy0,
        vz0,
        Th0,
        Er0,
        Frx0,
        Fry0,
        Frz0,
        hrho0,
        hvx0,
        hvy0,
        hvz0,
        hTh0,
        hEr0,
        hFrx0,
        hFry0,
        hFrz0,
        wall,
        gmax,
        mcool,
        hpcool,
        mclosure,
        hpclosure,
        fclosure,
        hprmodel,
        hsolver,
        rsolver,
        retry,
        outpath,
        outfile,
        ioformat,
        fps
    )

    manager = init(managerconfig)

    if lowercase(strip(outpath)) == "none" &&
       lowercase(strip(outfile)) == "none" &&
       lowercase(strip(ioformat)) == "none"

        io = nothing

    else
        io = IOParaView1D(outpath, outfile, ioformat)
    end

    run(manager, io)
end

if abspath(PROGRAM_FILE) == @__FILE__
    app()
end
