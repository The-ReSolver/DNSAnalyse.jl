# Functions to convert a DNSData object into an array of spectral field

function data2PhysicalField(data::DNSData{Ny, Nz, Nt}, grid::Grid{S}) where {Ny, Nz, Nt, S}
    S[1] == Ny || throw(ArgumentError("Wall normal discretisation is not compatible!"))
    u = VectorField(grid, fieldType=PhysicalField)
    for n in 1:3, (i, snap) in enumerate(data)
        u[n][:, :, i] .= snap[n]
    end
    return u
end

function data2SpectralField(data::DNSData{Ny, Nz, Nt}, grid::Grid{S}) where {Ny, Nz, Nt, S}
    S[1] == Ny || throw(ArgumentError("Wall normal discretisation is not compatible!"))
    u = VectorField(grid)
    A = [Array{Float64, 3}(undef, data.Ny, data.Nz, data.Nt) for _ in 1:3]
    for n in 1:3, (i, snap) in enumerate(data)
        A[n][:, :, i] .= snap[n]
    end
    B = [rfft(A[i], [2, 3])./(data.Nz*data.Nt) for i in 1:3]
    for n in 1:3, nz in 1:((minimum([S[2], Nz]) >> 1) + 1), nt in 1:((minimum([S[3], Nt]) >> 1) + 1)
        u[n][:, nz, nt] .= B[n][:, nz, nt]
        u[n][:, nz, end-nt+1] .= B[n][:, nz, end-nt+1]
    end
    return u
end

data2Coefficients(data::DNSData{Ny, Nz, Nt}, grid::Grid{S}, modes) where {Ny, Nz, Nt, S} = return expand!(SpectralField(grid, modes), data2SpectralField(data, grid), grid.ws, modes)


correct_mean!(u::VectorField{3, S}, y::Vector{Float64}) where {S<:SpectralField} = (u[1][:, 1, 1] .+= y; return u)
function correct_mean!(A::Array{Float64, 3}, y::Vector{Float64})
    for nt in axes(A, 2), nz in axes(A, 3)
        A[1][:, nz, nt] .+= y
    end
    return A
end
