module DNSAnalyse

# This purpose of this package is to provide a set of tools to easily analyse
# a set of DNS data.

# The main tool provided will be a generic interface to be able to generate an
# array based on the DNS data, relying on a user-provided method to unpack and
# interpret a single snapshot, as well as a list of the snapshots provided.

# In addition, there will be a large quantity of extra functionality to allow
# spectral analysis of the data.

# One of the targets of this package to remove the DNS data unpacking code from
# the Fields.jl package since it is becoming monolithic.

end
