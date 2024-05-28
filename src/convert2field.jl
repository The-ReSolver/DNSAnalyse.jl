# Functions to convert a DNSData object into an array of spectral field

function dsn2array!(A::Vector{Array{Float64, 3}}, data::DNSData)
    for (i, snap) in enumerate(data)
        A[1][:, :, i] = snap[1]
        A[2][:, :, i] = snap[2]
        A[3][:, :, i] = snap[3]
    end
    return correct_mean!(A, data)
end
dns2array(data::DNSData{Ny, Nz, Nt}) where {Ny, Nz, Nt} = dns2array!([Array{Float64, 3}(undef, Ny, Nz, Nt)], data)

function dns2spectralfield(data::DNSData{Ny, Nz, Nt}, size::NTuple{2, Int}=((Nz >> 1) + 1, Nt); diffmatrix_order::Int=5, quadweights_order=3) where {Ny, Nz, Nt}
    grid = Grid(data.y, size..., DiffMatrix(data.y, diffmatrix_order, 1), DiffMatrix(data.y, diffmatrix_order, 2), quadweights(data.y, quadweights_order), 2π/data.L, 2π/(data.snaps[end] - data.snaps[1]))
    u = VectorField(grid)
    A = dns2array(data)
    if size[1] <= Nz && size[2] <= Nt
        for i in 1:3
            u[i] .= rfft(A[i], [2, 3])[:, 1:((size[1] >> 1) + 1), 1:size[2]]
        end
    elseif size[1] <= Nz && size[2] > Nt
        for i in 1:3
            u[i][:, :, 1:Nt] .= rfft(A[i], [2, 3])[:, 1:((size[1] >> 1) + 1), :]
        end
    elseif size[1] > Nz && size[2] <= Nt
        for i in 1:3
            u[i][:, 1:((Nz >> 1) + 1), :] .= rfft(A[i], [2, 3])[:, :, 1:size[2]]
        end
    else
        for i in 1:3
            u[i][:, 1:((Nz >> 1) + 1), 1:Nt] .= rfft(A[i], [2, 3])[:, :, :]
        end
    end
    return correct_mean!(u, data)
end

correct_mean!(u::VectorField{3, S}, data::DNSData{Ny, Nz, Nt}) where {Ny, Nz, Nt, S<:SpectralField{Ny, Nz, Nt}} = (u[1][:, 1, 1] .+= data.y; return u)
function correct_mean!(A::Array{Float64, 3}, data::DNSData{Ny, Nz, Nt}) where {Ny, Nz, Nt}
    for nt in 1:Nt, nz in 1:Nz
        A[1][:, nz, nt] .+= data.y
    end
    return A
end
