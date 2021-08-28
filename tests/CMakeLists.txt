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



ADDTEST(uniform)
ADDTEST(normal)
ADDTEST(exponential)
ADDTEST(gamma)
ADDTEST(beta)
ADDTEST(binomial)
ADDTEST(geometric)
ADDTEST(hypergeometric)
ADDTEST(logarithmic)
ADDTEST(negbinomial)
ADDTEST(poisson)
