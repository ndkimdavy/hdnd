#===================================================================
                                Log
===================================================================#

const LOG = get(ENV, "LOG", "all")

function log(msg)
    if lowercase(LOG) != "none"
        println(msg)
        flush(stdout)
    end
end

const CLOCK_MONOTONIC = 1

struct Clock
    sec::Int64
    nsec::Int64
end

function clock_ns()
    clock = Ref(Clock(0, 0))

    ret = ccall(
        :clock_gettime,
        Cint,
        (Cint, Ref{Clock}),
        CLOCK_MONOTONIC,
        clock
    )

    ret != 0 && error("clock_gettime failed")

    return UInt64(clock[].sec) * UInt64(1_000_000_000) + UInt64(clock[].nsec)
end

function progress(label::String, ratio::Float64, elapsed::Float64; width::Int=30, height::Int=1)

    lwlog = lowercase(LOG)
    lwlabel = lowercase(label)

    if lwlog == "none"
        return
    end

    tags = strip.(split(lwlog, ","))

    if lwlog != "all" && !any(tag -> tag == lwlabel, tags)
        return
    end

    time = (sec::Float64) -> begin
        sec = max(sec, 0.0)

        h = Int(floor(sec / 3600.0))
        sec -= 3600.0 * h

        m = Int(floor(sec / 60.0))
        sec -= 60.0 * m

        s = Int(floor(sec))
        ms = Int(floor(1000.0 * (sec - s)))

        return string(
            lpad(h, 2, '0'), ":",
            lpad(m, 2, '0'), ":",
            lpad(s, 2, '0'), ":",
            lpad(ms, 3, '0')
        )
    end

    _ratio = clamp(ratio, 0.0, 1.0)

    nfill = Int(floor(width * _ratio))
    bar = repeat("█", nfill) * repeat("░", width - nfill)

    pct = _ratio >= 1.0 ? "100" : string(round(100.0 * _ratio, digits=2))
    eta = _ratio > 0.0 ? elapsed * (1.0 - _ratio) / _ratio : 0.0

    print(
        rpad(label, 32), " [", bar, "] ",
        lpad(pct, 6), "% ",
        "T+=", time(elapsed), " ",
        "ETA=", time(eta),
        repeat("\n", height)
    )
    flush(stdout)
end
