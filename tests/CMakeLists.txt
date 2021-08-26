### Pre-process: .fypp -> .f90 via Fypp

# Create a list of the files to be preprocessed
set(fyppFiles
    test_uniform.fypp
    test_normal.fypp
    test_exponential.fypp
    test_gamma.fypp
    test_beta.fypp
    test_binomial.fypp
    test_geometric.fypp
    test_hypergeometric.fypp
    test_logarithmic.fypp
    test_negbinomial.fypp
    test_poisson.fypp
)



fypp_f90("${fyppFiles}" outFiles)

addtest(uniform)
addtest(normal)
addtest(exponential)
addtest(gamma)
addtest(beta)
addtest(binomial)
addtest(geometric)
addtest(hypergeometric)
addtest(logarithmic)
addtest(negbinomial)
addtest(poisson)
