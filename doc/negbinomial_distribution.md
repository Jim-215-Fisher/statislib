---
title: Statislib
---

# Statistical Distribution - Negative Binomial Distribution

[TOC]

## `negbinomial_distribution_rvs` - Negative binomial distribution random variates

### Status

Release

### Description

A negative binomial discrete random variate distribution NB(n, p), also known as Pascal distribution, is used to characterize the number of failure in a sequence of Bernoullie trials before a specified number (n) of success occurs.

### Syntax

`result = [[statislib(module):negbinomial_distribution_rvs(interface)]](n, p, [array_size])`

### Class

Elemental function

### Arguments

`n`: has `intent in` and is a scalar of type `integer`.

`p`: has `intent in` and is a scalar of type `real` with single precision.

`array_size`: optional argument has `intent in` and is a scalar of type `integer`.

### Return value

The result is a scalar or rank one array with a size of `array_size` of type `integer`.

### Example

```fortran
program demo_nb_rvs
    use stdlib_stats_distribution_PRNG, only : rand_seed
    use statislib, only: nb_rvs => negbinomial_distribution_rvs

    implicit none
    integer :: n(2,3,4)
    real :: p(2,3,4)
    integer :: put, get

    put = 1234567
    call random_seed(put, get)

    n(:,:,:) = 20
    p(:,:,:) = 0.3
    print *, nb_rvs(n, p)       ! a rank 3 array of 60 nb random variate

!53          49          39          38          36          37          !43          55          34          37          62          52          !63          51          37          31          37          55          !38          33          60          33          53          39

    print *, nb_rvs(20, 0.3, 10)    ! an array of 10 nb random variates

!54          45          38          49          50          65          !30          50          38          32

end program demo_nb_rvs
```

## `negbinomial_distribution_pmf` - Negative binomial probability mass function

### Status

Release

### Description

The probability mass function of the discrete negative binomial distribution.

f(k; n, p) = C(k+n-1, n-1) p<sup>n</sup> (1-p)<sup>k</sup>

where C is the binomial coefficient

### Syntax

`result = [[statislib(module):negbinomial_distribution_pmf(interface)]](k, n, p)`

### Class

Elemental function

### Arguments

`k`: has `intent in` and is a scalar of type `integer`.

`n`: has `intent in` and is a scalar of type `integer`.

`p`: has `intent in` and is a scalar of type `real`.

Arguments `k` and `n` must have the same type.

### Return value

The result is a scalar of type `real` with a shape conformable to auguments.

### Example

```fortran
program demo_nb_pmf
    use stdlib_stats_distribution_PRNG, only : rand_seed
    use statislib, only: nb_pmf => negbinomial_distribution_pmf,         &
                         nb_rvs => negbinomial_distribution_rvs

    implicit none
    real :: p(2,3,4)
    integer :: n(2,3,4)
    integer :: put, get

    put = 1234567
    call random_seed(put, get)

    n(:,:,:) = 20
    p(:,:,:) = 0.7
    print *, nb_pmf(5, 20, 0.7)    ! a probability density for 5 in nb

!8.24131072E-02

    print *, nb_pmf(nb_rvs(n, p), n, p)
    ! a rank 3 array of nb probability density

!7.72150382E-02   9.43747014E-02  0.103016488      0.103016488      !0.103016488      0.114789672       7.72150382E-02   5.98418787E-02   !8.24131072E-02  0.103016488       2.12494712E-02   9.43747014E-02   !3.12492102E-02  0.116224766      0.103016488      0.103016488      !0.103016488       5.98418787E-02  0.114789672       8.24131072E-02   !5.98418787E-02   8.24131072E-02   4.41910028E-02   8.24131072E-02

end program demo_nb_pmf
```

## `negbinomial_distribution_cdf` - Negative binomial cumulative distribution function

### Status

Release

### Description

The cumuative distribution function of the discrete negative binomial distribution.

F(k; n, p) = 1 - B(p, K+1, n) / B(k+1, n)

where B is the regularized incomplet beta function$$

### Syntax

`result = [[statislib(module):negbinomial_distribution_cdf(interface)]](k, n, p)`

### Class

Elemental function

### Arguments

`k`: has `intent in` and is a scalar of type `integer`.

`n`: has `intent in` and is a scalar of type `integer`.

`p`: has `intent in` and is a scalar of type `real`.

Arguments `k` and `n` must have the same type.

### Return value

The result is a scalar of type `real` with a shape conformable to auguments.

### Example

```fortran
program demo_nb_cdf
    use stdlib_stats_distribution_PRNG, only : rand_seed
    use statislib, only: nb_cdf => negbinomial_distribution_cdf,         &
                         nb_rvs => negbinomial_distribution_rvs

    implicit none
    real :: p(2,3,4)
    integer :: n(2,3,4)
    integer :: put, get

    put = 1234567
    call random_seed(put, get)

    n(:,:,:) = 10
    p(:,:,:) = 0.83
    print *, nb_cdf(5, 10, 0.83) ! total probability for k not greater than 5

!4.41051535E-02

    print *, nb_cdf(nb_rvs(n,p), n, p)
    ! a rank 3 array of nb cumulative probability

!1.87446903E-02  0.167706698      0.263772994      0.263772994      !0.155160442      0.155160442      0.246627286       9.26578715E-02  !0.263772994      0.263772994      0.246627286       4.41051535E-02   !9.26578715E-02  0.246627286      0.263772994      0.155160442      !0.263772994      0.167706698      0.155160442      0.155160442       !9.26578715E-02  0.263772994      0.246627286      0.263772994

end program demo_nb_cdf
```
