---
title: Statislib
---

# Statistical Distribution -- Poisson Distribution

[TOC]


## `poisson_distribution_rvs` - poisson distribution random variates

### Status

Release

### Description

A poisson discrete random variate distribution Pois(&lambda;) is probability distribution of number of events with fixed mean rate &lambda; in a specified interval, such as time, distance, area or volume. The distribution is supported on the set {0,1,2,3...}.

### Syntax

`result = [[statislib(module):poisson_distribution_rvs(interface)]](lamdba, [array_size])`

### Class

Elemental function

### Arguments

`lambda`: has `intent in` and is a scalar of type `real`.

`array_size`: optional argument has `intent in` and is a scalar of type `integer`.

### Return value

The result is a scalar or rank one array with a size of `array_size` of type `integer`.

### Example

```fortran
program demo_poisson_rvs
    use stdlib_stats_distribution_PRNG, only : random_seed
    use statislib, only: rpois => poisson_distribution_rvs

    implicit none
    real :: p(2,3,4)
    integer :: put, get

    put = 1234567
    call random_seed(put, get)

    p(:,:,:) = 3.4
    print *, pois(p)       ! a rank 3 array of 60 poisson random variate

! [3, 0, 4, 4, 3, 5, 4, 1, 2, 5, 2, 2, 8, 7, 3, 7, 3, 3, 6, 7, 5, 2, 3, 4]

    print *, pois(3.4, 10)     ! an array of 10 poisson random variates

! [5, 3, 2, 5, 2, 3, 3, 5, 7, 3]

end program demo_poisson_rvs
```

## `poisson_distribution_pmf` - poisson probability mass function

### Status

Release

### Description

The probability mass function of the discrete poisson distribution.

f(x) = &lambda;<sup>k</sup> e<sup> - &lambda;</sup> / k!

### Syntax

`result = [[statislib(module):poisson_distribution_pmf(interface)]](k, lamda)`

### Class

Elemental function

### Arguments

`k`: has `intent in` and is a scalar of type `integer`.

`lambda`: has `intent in` and is a scalar of type `real`.

### Return value

The result is a scalar of type `real` with a shape conformable to auguments.

### Example

```fortran
program demo_pois_pmf
    use stdlib_stats_distribution_PRNG, only : rand_seed
    use statislib, only:         pois_pmf => poisson_distribution_pmf,  &
                                 pois_rvs => poisson_distribution_rvs

    implicit none
    real :: p(2,3,4)
    integer :: put, get

    put = 1234567
    call random_seed(put, get)

    p(:,:,:) = 5.3
    print *, pois_pmf(5, 5.3)    ! a probability density for 5 in poisson

! 0.173955172

    print *, pois_pmf(pois_rvs(p), p)
    ! a rank 3 array of poisson probability density

! [0.173955172, 2.64554471E-02, 0.153660342, 0.153660342, 0.173955172,
!  0.116342872, 0.153660342, 7.01069236E-02, 0.164108634, 0.116342872,
!  0.164108634, 0.164108634, 5.11932047E-03, 2.40566339E-02, 0.164108634,
!  4.53898683E-02, 0.164108634, 0.173955172, 4.53898683E-02, 2.40566339E-02,
!  7.70771652E-02, 0.164108634, 0.173955172, 0.153660342]

end program demo_pois_pmf
```

## `poisson_distribution_cdf` - poisson cumulative distribution function

### Status

Release

### Description

The cumuative distribution function of the discrete poisson distribution.

F(x) = e<sup>&lambda;</sup> &sum;<sup>x</sup> &lambda;<sup>i</sup> / i!

### Syntax

`result = [[statislib(module):poisson_distribution_cdf(interface)]](k, lamda)`

### Class

Elemental function

### Arguments

`k`: has `intent in` and is a scalar or an array of type `integer`.

`lambda`: has `intent in` and is a scalar of type `real`.

### Return value

The result is a scalar of type `real` with a shape conformable to auguments.

### Example

```fortran
program demo_pois_cdf
    use stdlib_stats_distribution_PRNG, only : rand_seed
    use statislib, only:         pois_cdf => poisson_distribution_cdf, &
                                 pois_rvs => poisson_distribution_rvs
    implicit none
    real :: p(2,3,4)
    integer :: put, get

    put = 1234567
    call random_seed(put, get)

    p(:,:,:) = 6.7
    print *, pois_cdf(5, 6.7)  ! total probability for k not greater than 5

! 0.340649486

    print *, pois_cdf(pois_rvs(p), p)
    ! a rank 3 array of poisson cumulative probability

! [0.495297134, 9.88079906E-02, 0.767283678, 0.767283678, 0.643316984,
!  0.921401739, 0.859569907, 0.202159047, 0.495297134, 0.921401739,
!  0.495297134, 0.340649486, 9.88079906E-02, 9.88079906E-02, 0.495297134,
!  3.71058583E-02, 0.495297134, 0.643316984, 9.47803259E-03, 3.71058583E-02,
!  0.959062934, 0.495297134, 0.643316984, 0.859569907]

end program demo_pois_cdf
```
