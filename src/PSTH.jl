
export psth, segment, NoGrp, Trialwise, Binwise

# enumeration for grouping dimention
@enum GrpDim NoGrp Trialwise Binwise

# ============================================================================ #
function psth(ts::Vector{Float64}, evt::Vector{Float64},
    kbin::AbstractArray{Int,1}, bin_size::Float64, grp_dim::GrpDim = NoGrp)

    # default pad is 1 to correctly convert timestamps to 1-based indicies
    pad = kbin[1] < 0 ? -kbin[1] : 1

    # make sure all time stamps are sufficiently > 0, otherwise we end up
    # trying to use 0 as an index...
    if (ts[1] < bin_size) || (evt[1] < bin_size)
        ts .+= bin_size
        evt .+= bin_size
    end

    # convert timestamps to indicies based on bin size
    # pad indicies to make sure the <kbin> index vector can be subtracted from
    # event indicies without producing indicies < 1
    kspk = floor.(Int, ts ./ bin_size) .+ pad

    # if ts and evt refer to the same memory address (i.e. we're doing an
    # autocorrelation) then just have kevt and kspk refer to the same array
    if ts === evt
        kevt = kspk
    else
        kevt = floor.(Int, evt ./ bin_size) .+ pad
    end

    # maximum time that we'll be dealing with, gives us the length of the
    # "binary spike vector"
    max_time = kevt[end] + kbin[end] + 1
    ke = searchsortedlast(kspk, max_time)

    # Int8 is acceptable so long as we are assiging 1 to the <kspk> elements
    # (see warning below)
    bsv = zeros(Int, max_time)

    # WARNING: this step ensures that bsv is binary even if multiple spikes
    # would fall into the same bin (i.e. multiple instances of the same index in
    # <kspk>)
    # bsv[kspk[1:ke]] .= 1
    @inbounds @fastmath for k = 1:ke
        bsv[kspk[k]] += 1
    end

    data = zeros(length(kbin), length(kevt))
    @inbounds @fastmath for k in eachindex(evt)
        data[:,k] .= bsv[kevt[k] .+ kbin]
    end

    if grp_dim != NoGrp
        # -------------------------------------------------------------------- #
        # NOTE: this is a crappy, repetitive, inelegant solution
        # -------------------------------------------------------------------- #
        tmp = zeros(Int, max_time)
        tmp[kspk[1:ke]] = 1:ke

        tmp_data = zeros(Int, length(kbin), length(kevt))

        @inbounds @fastmath for k in eachindex(kevt)
            tmp_data[:,k] = tmp[kevt[k] .+ kbin]
        end

        # collapse across the non-requested dimention
        if grp_dim == Trialwise
            grp_length = length(evt)
            kdim = 2
        else
            grp_length = length(kbin)
            kdim = 1
        end

        grp = Vector{Vector{Int64}}(undef, size(tmp_data, kdim))
        @inbounds for k in eachindex(grp)
            grp[k] = rmzeros(vec(selectdim(tmp_data, kdim, k)))
        end
    else
        grp = Vector{Vector{Int64}}(undef, 0)
    end

    return data, grp
end
# ============================================================================ #
rmzeros(x::AbstractVector{<:Real}) = x[findall(!iszero, x)]
# ============================================================================ #
function segment(ts::Vector{Float64}, evt::Vector{Float64}, label::Vector{L},
    dur::Float64, bin_size::Float64) where L
    return segment(ts, evt, label, 0.0, dur, bin_size)
end
# ---------------------------------------------------------------------------- #
function segment(ts::Vector{Float64}, evt::Vector{Float64}, label::Vector{L},
    pre::Float64, post::Float64, bin_size::Float64) where L

    if length(label) != length(evt)
        error("Length of event label and event timestamp vectors must be equal")
    end

    kbins = round.(Int, -abs(pre)/bin_size) : round(Int, post/bin_size)
    data = psth(ts, evt, kbins, bin_size)[1]

    unique_label::Vector{L} = sort(unique(label))

    nlabel = length(unique_label)
    grouped_data = Vector{Matrix{Float64}}(undef, nlabel)

    @inbounds for k = 1:nlabel
        group = findall(isequal(unique_label[k]), label)
        grouped_data[k] = data[:, group]
    end

    return grouped_data, unique_label
end
# ============================================================================ #
