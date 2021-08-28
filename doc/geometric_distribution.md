---
title: Statislib
---

# Statistical Distribution -- Geometric Distribution

[TOC]


## `geometric_distribution_rvs` - geometric distribution random variates

### Status

Release

### Description

A geometric discrete random variate distribution `G(p)` is to model the number of failure before the first success. It is probability distribution supported on the set {0,1,2,3...}.

### Syntax

`result = [[statislib(module):geometric_distribution_rvs(interface)]](p, [array_size])`

### Class

Elemental function

### Arguments

`p`: has `intent in` and is a scalar of type `real`.

`array_size`: optional argument has `intent in` and is a scalar of type `integer`.

### Return value

The result is a scalar or rank one array with a size of `array_size` of type `integer`.

### Example

```fortran
program demo_geometric_rvs
    use stdlib_stats_distribution_PRNG, only : rand_seed
    use statislib, only: geo => geometric_distribution_rvs

    implicit none
    real :: p(3,4,5), rv(10)
    integer :: put, get

    put = 1234567
    call random_seed(put, get)

    p(:,:,:) = 0.13
    print *, geo(p)       ! a rank 3 array of 60 geometric random variate

! [13, 4, 0, 0, 15, 1, 1, 5, 9, 2, 8, 7, 4, 3, 12, 3, 11, 18, 3, 3, 2,
!  10, 28, 1, 2, 14, 7, 1, 7, 11, 16, 1, 3, 14, 0, 0, 9, 18, 1, 8, 4, 6,
!   0, 18, 8, 2, 13, 16, 1, 9, 8, 15, 4, 2, 1, 8, 28, 3, 2, 11]

    print *, geo(0.13, 10)     ! an array of 10 geometric random variates

! [6, 15, 1, 1, 1, 7, 0, 6, 3, 0]

end program demo_geometric_rvs
```

## `geometric_distribution_pmf` - geometric probability mass function

### Status

Release

### Description

The probability mass function of the discrete geometric distribution.

f(x) = (1 - p)<sup>k</sup>

### Syntax

`result = [[statislib(module):geometric_distribution_pmf(interface)]](k, p)`

### Class

Elemental function

### Arguments

`k`: has `intent in` and is a scalar of type `integer`.

`p`: has `intent in` and is a scalar of type `real`.


### Return value

The result is a scalar of type `real` with a shape conformable to auguments.

### Example

```fortran
program demo_geo_pmf
    use stdlib_stats_distribution_PRNG, only : rand_seed
    use statislib, only: geo_pmf => geometric_distribution_ pmf,&
                         geo => geometric_distribution_rvs

    implicit none
    real :: p(2,3,4)
    integer :: m(2,3,4)
    integer :: put, get

    put = 1234567
    call random_seed(put, get)

    p(:,:,:) = 0.1
    m = geo(p)
    print *, geo_pmf(5, 0.1)    ! a probability density for 5 in geometric

! 5.90489879E-02

    print *, geo_pmf(m, p)      ! a rank 3 array of geometric probability density

! [1.66771747E-02, 5.90489879E-02, 9.99999940E-02, 9.99999940E-02,
!  1.21576609E-02, 8.10000002E-02, 8.99999961E-02, 4.78296801E-02,
!  3.13810445E-02, 8.10000002E-02, 3.13810445E-02, 3.48678306E-02,
!  5.90489879E-02, 6.56099841E-02, 1.85301956E-02, 6.56099841E-02,
!  2.28767842E-02, 8.86293128E-03, 6.56099841E-02, 6.56099841E-02,
!  7.28999972E-02, 2.54186485E-02, 2.02755351E-03, 8.99999961E-02]

end program demo_geo_pmf
```

## `geometric_distribution_cdf` - geometric cumulative distribution function

### Status

Release

### Description

The cumuative distribution function of the discrete geometric distribution.

F(x) = 1 - (1 - p)<sup>k + 1</sup>

### Syntax

`result = [[statislib(module):geometric_distribution_cdf(interface)]](k, p)`

### Class

Elemental function

### Arguments

`k`: has `intent in` and is a scalar or an array of type `integer`.

`p`: has `intent in` and is a scalar of type `real`.

### Return value

The result is a scalar of type `real` with a shape conformable to auguments.

### Example

```fortran
program demo_geo_cdf
    use stdlib_stats_distribution_PRNG, only : rand_seed
    use statislib, only: geo_cdf => geometric_distribution_cdf, &
                         geo => geometric_distribution_rvs

    implicit none
    real :: p(2,3,4)
    integer :: m(2,3,4)
    integer :: put, get

    put = 1234567
    call random_seed(put, get)

    p(:,:,:) = 0.12
    m = geo(p)
    print *, geo_cdf(5, 0.12)  ! total probability for k not greater than 5

! 0.535595953

    print *, geo_cdf(m, p)  ! a rank 3 array of geometric cumulative probability

! [2.00418849E-02, 7.19634444E-02, 0.119999997, 0.119999997, 1.55204320E-02,
!  0.105599999, 0.105599999, 5.57284914E-02, 3.79774049E-02, 9.29279998E-02,
!  3.79774049E-02, 4.31561433E-02, 7.19634444E-02, 7.19634444E-02,
!  2.27748696E-02, 8.17766413E-02, 2.58805323E-02, 1.05767399E-02,
!  8.17766413E-02, 8.17766413E-02, 9.29279998E-02, 2.94097029E-02,
!  2.59215664E-03, 0.105599999]

end program demo_geo_cdf
```
