include("log.jl")

#===================================================================
                              Reject
===================================================================#

const REJECT = get(ENV, "REJECT", "all")

mutable struct Reject
    flag::Bool
    msg::String
end

function reject!(r::Reject, msg::String)
    r.flag = true
    r.msg = msg

    lwreject = lowercase(REJECT)
    lwwords = lowercase.(split(msg))

    if lwreject == "none"
        return
    end

    tags = strip.(split(lwreject, ","))

    if lwreject != "all" && !any(tag -> tag in lwwords, tags)
        return
    end

    log("REJECT: $(msg)")
end

function reset!(r::Reject)
    r.flag = false
    r.msg = ""
end
