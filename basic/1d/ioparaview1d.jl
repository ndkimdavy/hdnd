include("../../interface/iio.jl")
include("../mem.jl")
include("../log.jl")

using Printf
using WriteVTK
using HDF5

#===================================================================
                           IO ParaView 1D
===================================================================#

struct IOParaView1D <: IIO
    outpath::String
    outfile::String
    ioformat::String
end

function abort(io::IOParaView1D, msg::String)

    log("ERROR: IOParaView1D")
    log("io       = $(typeof(io))")
    log("outpath  = $(io.outpath)")
    log("outfile  = $(io.outfile)")
    log("ioformat = $(io.ioformat)")

    error("$(msg)")
end

function fvtr(io::IOParaView1D, it::Int, t::T, x::Array{T,1}, nx::Int, dbmem::DBMem{T}) where {T<:AbstractFloat}

    outpath = io.outpath
    outfile = io.outfile
    data = dbmem.data

    xnode = zeros(T, nx + 1)

    for j in 2:nx
        xnode[j] = T(0.5) * (x[j-1] + x[j])
    end

    xnode[1] = x[1] - (xnode[2] - x[1])
    xnode[end] = x[end] + (x[end] - xnode[end-1])

    dy = xnode[2] - xnode[1]
    ynode = T[zero(T), dy]

    fsnap = @sprintf("%s_%08d", outfile, it)
    fvtk = joinpath(outpath, fsnap)
    fpvd = joinpath(outpath, "$(outfile).pvd")

    vtk_grid(fvtk, xnode, ynode) do vgrid
        vgrid["time", VTKFieldData()] = T[t]
        vgrid["iteration", VTKFieldData()] = Int[it]

        for (key, value) in data

            if key in (:it, :t, :dt)
                continue
            end

            vgrid[string(key), VTKCellData()] = reshape(value, nx, 1)
        end
    end

    line = "    <DataSet timestep=\"$(Float64(t))\" group=\"\" part=\"0\" file=\"$(fsnap).vtr\"/>"

    if isfile(fpvd)
        lines = readlines(fpvd)

        if !any(item -> occursin("$(fsnap).vtr", item), lines)
            iend = findfirst(item -> occursin("</Collection>", item), lines)

            if iend === nothing
                abort(io, "IOParaView1D invalid pvd file")
            end

            insert!(lines, iend, line)

            open(fpvd, "w") do out
                for item in lines
                    println(out, item)
                end
            end
        end

    else
        open(fpvd, "w") do out
            println(out, "<?xml version=\"1.0\"?>")
            println(out, "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\">")
            println(out, "  <Collection>")
            println(out, line)
            println(out, "  </Collection>")
            println(out, "</VTKFile>")
        end
    end

    return
end

function fvtkhdf(io::IOParaView1D, it::Int, t::T, x::Array{T,1}, nx::Int, dbmem::DBMem{T}) where {T<:AbstractFloat}

    outpath = io.outpath
    outfile = io.outfile
    data = dbmem.data

    xnode = zeros(T, nx + 1)

    for j in 2:nx
        xnode[j] = T(0.5) * (x[j-1] + x[j])
    end

    xnode[1] = x[1] - (xnode[2] - x[1])
    xnode[end] = x[end] + (x[end] - xnode[end-1])

    dy = xnode[2] - xnode[1]

    points = zeros(T, 3, 2 * (nx + 1))
    connectivity = zeros(Int64, 4 * nx)
    offsets = zeros(Int64, nx + 1)
    types = fill(UInt8(9), nx)

    for j in 1:(nx+1)
        points[1, j] = xnode[j]
        points[2, j] = zero(T)
        points[3, j] = zero(T)

        points[1, nx+1+j] = xnode[j]
        points[2, nx+1+j] = dy
        points[3, nx+1+j] = zero(T)
    end

    offsets[1] = 0

    for j in 1:nx
        connectivity[4*j-3] = j - 1
        connectivity[4*j-2] = j
        connectivity[4*j-1] = nx + 1 + j
        connectivity[4*j] = nx + j

        offsets[j+1] = 4 * j
    end

    fvtkhdf = joinpath(outpath, "$(outfile).vtkhdf")

    h5open(fvtkhdf, "cw") do h5

        if !haskey(h5, "VTKHDF")
            root = create_group(h5, "VTKHDF")

            attributes(root)["Version"] = Int64[1, 0]
            attributes(root)["Type"] = "UnstructuredGrid"

            HDF5.write(root, "NumberOfPoints", Int64[2*(nx+1)])
            HDF5.write(root, "NumberOfCells", Int64[nx])
            HDF5.write(root, "NumberOfConnectivityIds", Int64[4*nx])

            HDF5.write(root, "Points", points)
            HDF5.write(root, "Connectivity", connectivity)
            HDF5.write(root, "Offsets", offsets)
            HDF5.write(root, "Types", types)

            create_group(root, "CellData")

            steps = create_group(root, "Steps")

            HDF5.write(steps, "Values", T[])
            HDF5.write(steps, "PartOffsets", Int64[])
            HDF5.write(steps, "NumberOfParts", Int64[])
            HDF5.write(steps, "PointOffsets", Int64[])
            HDF5.write(steps, "CellOffsets", reshape(Int64[], 1, 0))
            HDF5.write(steps, "ConnectivityIdOffsets", reshape(Int64[], 1, 0))

            create_group(steps, "CellDataOffsets")
        end

        root = h5["VTKHDF"]
        celldata = root["CellData"]
        steps = root["Steps"]
        celldataoffsets = steps["CellDataOffsets"]

        values = read(steps["Values"])
        partoffsets = read(steps["PartOffsets"])
        numberofparts = read(steps["NumberOfParts"])
        pointoffsets = read(steps["PointOffsets"])
        celloffsets = read(steps["CellOffsets"])
        connectivityidoffsets = read(steps["ConnectivityIdOffsets"])

        nstep = length(values)

        values = vcat(values, T[t])
        partoffsets = vcat(partoffsets, Int64(0))
        numberofparts = vcat(numberofparts, Int64(1))
        pointoffsets = vcat(pointoffsets, Int64(0))
        celloffsets = hcat(celloffsets, Int64[0])
        connectivityidoffsets = hcat(connectivityidoffsets, Int64[0])

        delete_object(steps, "Values")
        delete_object(steps, "PartOffsets")
        delete_object(steps, "NumberOfParts")
        delete_object(steps, "PointOffsets")
        delete_object(steps, "CellOffsets")
        delete_object(steps, "ConnectivityIdOffsets")

        HDF5.write(steps, "Values", values)
        HDF5.write(steps, "PartOffsets", partoffsets)
        HDF5.write(steps, "NumberOfParts", numberofparts)
        HDF5.write(steps, "PointOffsets", pointoffsets)
        HDF5.write(steps, "CellOffsets", celloffsets)
        HDF5.write(steps, "ConnectivityIdOffsets", connectivityidoffsets)

        if "NSteps" in keys(attributes(steps))
            delete_attribute(steps, "NSteps")
        end

        attributes(steps)["NSteps"] = Int64(nstep + 1)

        for (key, value) in data

            if key in (:it, :t, :dt)
                continue
            end

            skey = string(key)

            if haskey(celldata, skey)
                dataall = read(celldata[skey])
                delete_object(celldata, skey)
            else
                dataall = T[]
            end

            if haskey(celldataoffsets, skey)
                offsetall = read(celldataoffsets[skey])
                delete_object(celldataoffsets, skey)
            else
                offsetall = Int64[]
            end

            offsetall = vcat(offsetall, Int64(length(dataall)))
            dataall = vcat(dataall, value)

            HDF5.write(celldata, skey, dataall)
            HDF5.write(celldataoffsets, skey, offsetall)
        end
    end

    return
end

function write(io::IOParaView1D, dbmem::DBMem{T}) where {T<:AbstractFloat}

    outpath = io.outpath
    ioformat = io.ioformat
    data = dbmem.data

    it = Int(round(data[:it][1]))
    t = data[:t][1]

    x = data[:x]
    nx = length(x)

    mkpath(outpath)

    if lowercase(ioformat) == "vtr"
        fvtr(io, it, t, x, nx, dbmem)

    elseif lowercase(ioformat) == "vtkhdf"
        fvtkhdf(io, it, t, x, nx, dbmem)

    else
        abort(io, "IOParaView1D invalid ioformat")
    end

    return nothing
end

function release(io::IOParaView1D)
    return nothing
end
