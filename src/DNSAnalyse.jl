module DNSAnalyse

using IniFile, FFTW

using Fields, FDGrids

export DNSData, loadDNS
export dns2field!, dns2field, correct_mean!
export mean!, mean

include("snapshoterror.jl")
include("snapshot.jl")
include("dnsdata.jl")
include("convert2field.jl")
include("mean.jl")

end
