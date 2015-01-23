set_cppo(mode = "fast")
# Example 2pl model

# Simulate some data
Y <- sim_irt(100, 10, seed = 234234)

out <- par2pl(Y, iter = 1000, chains = 1)
