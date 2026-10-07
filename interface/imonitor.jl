abstract type IMonitor end

function monitor!(monitor::IMonitor, args...)
    error("monitor! function not implemented for monitor=$(typeof(monitor)), args=$(typeof(args))")
end
