using Distributions
function serhash(x)
    b = IOBuffer()
    serialize(b, x)
    hash(take!(b))
end
ℓ = MvNormal([-1.42646, 0.94423, 0.852379, -1.12906, 0.0868619, 0.948781, -0.875067, 1.07243],
             [14.8357 2.42526 -2.97011 2.08363 -1.67358 4.02846 5.57947 7.28634;
              2.42526 10.8874 -1.08992 1.99358 1.85011 -2.29754 -0.0540131 1.79718;
              -2.97011 -1.08992 3.05794 0.0321187 1.8052 -1.5309 1.78163 -0.0821483;
              2.08363 1.99358 0.0321187 2.38112 -0.252784 0.666474 1.73862 2.55874;
              -1.67358 1.85011 1.8052 -0.252784 12.3109 -2.3913 -2.99741 -1.95031;
              4.02846 -2.29754 -1.5309 0.666474 -2.3913 4.89957 3.6118 5.22626;
              5.57947 -0.0540131 1.78163 1.73862 -2.99741 3.6118 10.215 9.60671;
              7.28634 1.79718 -0.0821483 2.55874 -1.95031 5.22626 9.60671 11.5554])
println(hex(serhash(ℓ)))
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
