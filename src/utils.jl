# Utility functions.

function _mean!(ū::Vector{Float64}, data::DNSData{<:Any, Nz}, snap_times::Vector{Float64}) where {Nz}
    # loop over the time window of the data and compute mean
    for t in snap_times
        ū .+= dropdims(sum(data[t][1], dims=2), dims=2)./Nz
    end
    ū ./= length(snap_times)

    # add back the laminar profile
    ū .+= data.y

    return ū
end

function mean!(ū, data::DNSData; window::NTuple{2, Real}=(firstindex(data), lastindex(data)))
    # find range of snapshots inside the provided window
    snapshot_times = data.snaps
    start_ti = findfirst(x->window[1]<=x, snapshot_times)
    end_ti = findlast(x->window[2]>=x, snapshot_times)

    # overwrite the mean with zeros
    ū .= zero(Float64)

    # compute mean
    return _mean!(ū, data, snapshot_times[start_ti:end_ti])
end

mean(data::DNSData{Ny}; window::NTuple{2, Real}=(firstindex(data), lastindex(data))) where {Ny} = (ū = zeros(Float64, Ny); mean!(ū, data, window=window))


# -----------------------------------------------------------------------------
# Methods to convert simulation directory directly into fields which we know
# how to manipulate
# -----------------------------------------------------------------------------

# TODO: fix these
dns2field(loc::AbstractString; fft_flag::UInt32=ESTIMATE, times::Union{Nothing, NTuple{2, Real}}=nothing, skip_step::Int=1) = dns2field(DNSData(loc)[times, skip_step], fft_flag=fft_flag)

function dns2field(data::DNSData{Ny, Nz, Nt}; fft_flag::UInt32=ESTIMATE) where {Ny, Nz, Nt}
    # initialise grid, fields, and, and FFT plan
    grid = Grid(data.y, Nz, Nt, zeros(Ny, Ny), zeros(Ny, Ny), zeros(Ny), 2π/((data.snaps[2] - data.snaps[1])*Nt), data.β)
    u = VectorField(grid; field_type=:physical)
    U = VectorField(grid)
    FFT! = FFTPlan!(grid, flags=fft_flag)

    return dns2field!(U, u, FFT!, data)
end

dns2field!(U::VectorField{3, S},
            u::VectorField{3, P},
            FFT!::FFTPlan!{Ny, Nz, Nt},
            data::DNSData{Ny, Nz, Nt}) where {Ny, Nz, Nt, S<:SpectralField{Ny, Nz, Nt}, P<:PhysicalField{Ny, Nz, Nt}} = FFT!(U, dns2field!(u, data))

function dns2field!(U::VectorField{3, P}, data::DNSData{Ny, Nz, Nt}) where {Ny, Nz, Nt, P<:PhysicalField{Ny, Nz, Nt}}
    # loop over snaps and assign each velocity component
    for (i, snaps) in enumerate(data)
        U[1][:, :, i] .= snaps[1]
        U[2][:, :, i] .= snaps[2]
        U[3][:, :, i] .= snaps[3]
    end

    return correct_mean!(data, U)
end

correct_mean!(data::DNSData{Ny, Nz, Nt}, u::VectorField{3, S}) where {Ny, Nz, Nt, S<:SpectralField{Ny, Nz, Nt}} = (u[1][:, 1, 1] .+= data.y; return u)
function correct_mean!(data::DNSData{Ny, Nz, Nt}, u::VectorField{3, P}) where {Ny, Nz, Nt, P<:PhysicalField{Ny, Nz, Nt}}
    for nt in 1:Nt, nz in 1:Nz
        u[1][:, nz, nt] .+= data.y
    end
    
    return u
end
