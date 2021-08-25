---
title: statslib
---

# Statistical Distribution -- Binomial Distribution

[TOC]


## `binomial_distribution_rvs` - binomial distribution random variates

### Status

Release

### Description

A binomial discrete random variate distribution B(n, p). It is used to characterize the number of successes in a sequence of n Bernoulli trials, each with the same probablity of p.

For a single trial, binomial distribution is Bernoulli distribution.

### Syntax

`result = [[stdlib_stats_distribution_binomial(module):binomial_distribution_rvs(interface)]](n, p, [array_size])`

### Class

Elemental function

### Arguments

`n`: has `intent(in)` and is a scalar of type `integer`.

`p`: has `intent(in)` and is a scalar of type `real` with single precision.

`array_size`: optional argument has `intent(in)` and is a scalar of type `integer`.

### Return value

The result is a scalar or rank one array with a size of `array_size` of type `integer`.

### Example

```fortran
program demo_binomial_rvs
    use stdlib_stats_distribution_PRNG, only : random_seed
    use statislib_binomial, only: rbinom => binomial_distribution_rvs

    implicit none
    integer :: n(2,3,4)
    real :: p(2,3,4), rv(10)
	integer :: seed_put, seed_get

    seed_put = 1234567
    call random_seed(seed_put, seed_get)

    n(:,:,:) = 20
    p(:,:,:) = 0.3
    print *, rbinom(n, p)        ! a rank 3 array of 24 binomial random variate
	
! 4  6  11  9  4  8  8  6  5  7  5  5  6  6  4  7  4  3  7  6  7  4  2  8

    print *, rbinom(20, 0.3, 10) ! an array of 10 binomial random variates 
	
! 7  4  5  8  5  4  3  8  6  4

end program demo_binomial_rvs
```

## `binomial_distribution_pmf` - Binomial probability mass function

### Status

Release

### Description

The probability mass function of the discrete binomial distribution.

p(k) = C(n, k) p<sup>k</sup>q<sup>n-k</sup>, k = 0, 1, 2, &hellip;, n

where C(n, k) is the binomial coefficent

### Syntax

`result = [[stdlib_stats_distribution_binomial(module):binomial_distribution_pmf(interface)]](k, n, p)`

### Class

Elemental function

### Arguments

`k`: has `intent(in)` and is a scalar of type `integer`.

`n`: has `intent(in)` and is a scalar of type `integer`.

`p`: has `intent(in)` and is a scalar of type `real`, single precision.

Arguments `k` and `n` must have the same type.

### Return value

The result is a scalar of type `real` with a shape conformable to auguments.

### Example

```fortran
program demo_binom_pmf
    use stdlib_stats_distribution_PRNG, only : random_seed
    use statislib_binomial, only:                               &
	                                     binom_pmf => binomial_distribution_pmf,&
                                         rbinom => binomial_distribution_rvs
    implicit none
    real :: p(2,3,4)
    integer :: n(2,3,4), m(2,3,4)
	integer :: seed_put, seed_get

    seed_put = 1234567
    call random_seed(seed_put, seed_get)

    n(:,:,:) = 20
    p(:,:,:) = 0.73
    m = rbinom(n,p)
    print *, binom_pmf(5, 20, 0.4)    ! a probability density for 5 in binomial
	
! 7.46470168E-02

    print *, binom_pmf(m, n, p)   ! a rank 3 array of binomial probability density

! 0.106533177  0.183267891  1.63479988E-02  4.01819535E-02  0.106533177 
! 0.135568008  8.14800784E-02  0.198200837  0.167461380  0.135568008 
! 0.167461380  0.167461380  0.183267891  0.183267891  0.167461380 
! 0.183267891  0.167461380  0.106533177  0.183267891  0.183267891 
! 0.183267891  0.167461380  4.80056927E-02  8.14800784E-02

end program demo_binom_pmf
```

## `binomial_distribution_cdf` - Binomial cumulative distribution function

### Status

Release

### Description

The cumuative distribution function of the discrete binomial distribution.

F(k) = &sum;<sup>k</sup> C(n, i) p<sup>i</sup>q<sup>n-i</sup>, k = 0, 1, 2, &hellip;, n

where C(n, i) is the binomial coefficient

### Syntax

`result = [[stdlib_stats_distribution_binomial(module):binomial_distribution_cdf(interface)]](k, n, p)`

### Class

Elemental function

### Arguments

`k`: has `intent(in)` and is a scalar or an array of type `integer`.

`n`: has `intent(in)` and is a scalar of type `integer`.

`p`: has `intent(in)` and is a scalar of type `real`, single precision.

Arguments `k` and `n` must have the same type.

### Return value

The result is a scalar of type `real` with a shape conformable to auguments.

### Example

```fortran
program demo_binom_cdf
    use stdlib_stats_distribution_PRNG, only : random_seed
    use statislib_binomial, only:                               &
	                                     binom_cdf => binomial_distribution_cdf,&
	                                     rbinom => binomial_distribution_rvs
    implicit none
    real :: p(2,3,4,5)
    integer :: n(2,3,4), m(2,3,4)
	integer :: seed_put, seed_get

    seed_put = 1234567
    call random_seed(seed_put, seed_get)

    n(:,:,:) = 30
    p(:,:,:) = 0.43
    m = rbinom(n,p)
    print *, binom_cdf(5, 20, 0.4) ! total probability for k not greater than 5

! 0.125598967

    print *, binom_cdf(m, n, p) 
! a rank 3 array of binomial cumulative probability

! 0.590547740  0.444910079  0.188751101  0.590547740  0.907233179 
! 0.907233179  0.188751101  2.16887915E-03  0.305481762  0.444910079 
! 0.831309795  0.907233179  0.831309795  0.188751101  0.305481762 
! 0.444910079  0.992344320  0.992344320  0.444910079  0.188751101 
! 0.103646301  0.723957717  0.590547740  0.907233179

end program demo_binom_cdf
```
