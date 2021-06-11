
export filter_events_by_duration

# ============================================================================ #
function filter_events_by_duration(evt::Vector{<:Real}, dur::Real; blank::Real=NaN)
    df = diff(evt)
    evt_ts = Vector{Float64}(undef, 0)
    evt_off = Vector{Float64}(undef, 0)
    nblank = 0
    warning = false
    for k in eachindex(df)
        if dur * 0.95 <= df[k] <= dur * 1.05
            push!(evt_ts, evt[k])
        elseif blank * 0.95 <= df[k] <= blank * 1.05
            nblank += 1
            push!(evt_off, evt[k])
        elseif !isempty(evt_ts) || nblank > 0
            # we have seen correct events and now we got one of odd duration
            if k < length(df)
                @warn("Some events that *should* be valid have invalid duration")
                warning = true
            end
        end
    end
    return evt_ts, warning, evt_off
end
# ============================================================================ #
function filter_events_by_duration(evt::Vector{<:Real},
            dur::AbstractFloat, nevt::Integer, limits::Tuple{<:Real,<:Real}=(-0.05, 0.05))

    warning = false

    if length(evt) < nevt
        @warn("Number of event timestamps is < the given number of events")
        warning = true
    end

    pad = 1.0 .+ limits
    df = diff(evt)
    # bgood = (df .>= dur*pad[1]) .& (df .<= dur*pad[2])
    bgood = (dur*pad[1]) .<= df .<= (dur*pad[2])
    ngood = sum(bgood)

    evt_off = Vector{Float64}(undef, 0)

    if ngood == nevt-1
        evt_ts = evt[[bgood; true]]

    elseif ngood == nevt
        evt_ts = evt[[bgood; false]]

    else
        # proportion of events of duration dur
        pdur = ngood ./ length(df)

        if pdur <= .1
            # few / no events are of length dur: events *MUST* only mark stim
            # onset
            total_dur = nanmedian(df)
            # buse = (df .>= total_dur*pad[1]) .& (df .<= total_dur*pad[2])
            buse = (total_dur*pad[1]) .<= df .<= (total_dur*pad[2])
            nuse = sum(buse)

            if nuse == nevt
                buse = [buse; false]
            elseif nuse == nevt-1
                buse = [buse; true]
            elseif nuse < nevt-1
                buse = [buse; true]
                @warn("Failed to locate the correct number of events")
                @show(length(evt), nevt, sum(buse))
                warning = true
            end

            evt_ts = evt[buse]

        elseif pdur >= .9
            if ngood + 1 >= nevt*2
                # b/c of the number of events, how we mark the final event
                # makes no difference
                push!(bgood, true)

                # most / all events are of length dur: events mark on and offset
                # of stim and blank_dur == stim_dur
                evt_ts = evt[bgood]

                #index of first stim on event (by definition precedes stim off)
                # buse = falses(length(evt_ts))
                # buse[1:2:(nevt*2)] .= true
                # evt_ts = evt_ts[buse]
                kon = 1:2:(nevt*2)
                evt_ts = evt_ts[kon]
                evt_off = evt_ts[kon .+ 1]

            elseif ngood > nevt
                # b/c of the number of events, how we mark the final event
                # makes no difference
                push!(bgood, true)

                # events mark stim on but we have extra events of approx. the
                # correct duration so just take the first nevt...
                @warn("Correct events could not be unambiguously identified...")
                @warn("falling back to heuristic method 1. Check the results!")
                evt_ts = evt[bgood]
                evt_ts = evt_ts[1:nevt]
                warning = true
            else
                error("Failed to locate the correct number of events")
            end

        elseif 0.4 <= pdur <= 0.6
            # approx half of events are of the correct length, so just go with that...
            @warn("Correct events could not be unambiguously identified...")
            @warn("falling back to heuristic method 2. Check the results!")
            warning = true

            # event duration is our only diagnostic info, so the final event
            # cannot mark a valid event onset
            push!(bgood, false)

            @show(length(evt), length(bgood), sum(bgood), nevt)
            evt_ts = evt[bgood]
            if length(evt_ts) > nevt
                evt_ts = evt_ts[1:nevt]
            end
        else
            error("Event durations are too inconsistent to estimate blank duration")
        end
    end

    # if warning
    #     @info("Requested dur: $(dur), median(diff(evt_ts)) = $(median(diff(evt_ts)))")
    # end

    return evt_ts, warning, evt_off
end
# ============================================================================ #
nanmedian(x::Vector{<:AbstractFloat}) = median(filter(!isnan, x))
# ============================================================================= #
