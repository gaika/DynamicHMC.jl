# include("setup-and-utilities.jl")
# # include("test-utilities.jl")
# # include("test-basics.jl")
# # include("test-Hamiltonian-leapfrog.jl")
# # include("test-buildingblocks.jl")
# # include("test-stepsize.jl")
# # include("test-sample-dummy.jl")
# # include("test-sample-normal.jl")
# include("test-z.jl")

# println("FINAL RNG $(@__FILE__) @ $(@__LINE__)")
# print_rng()

RNG = srand(UInt32[0x23ef614d, 0x8332e05c, 0x3c574111, 0x121aa2f4])
peekrand(rng) = println("rng" * hex(rand(copy(rng), UInt64)))

peekrand(RNG)
for i in 1:1000
    rand(RNG, 100000);
end
peekrand(RNG)
