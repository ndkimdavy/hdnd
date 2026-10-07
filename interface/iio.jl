abstract type IIO end

function abort(io::IIO, msg::String, args...)
    error("io=$(typeof(io)), msg=$(msg), args=$(args)")
end

function write(io::IIO, args...)
    error("write function not implemented for io=$(typeof(io)), args=$(typeof(args))")
end

function release(io::IIO)
    error("release function not implemented for io=$(typeof(io))")
end
