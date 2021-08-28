---
title: Statislib
---

# Statistical Distribution -- Hypergeometric Distribution

[TOC]

## `hypergeometric_distribution_rvs` - Hypergeometric distribution random variates

### Status

Release

### Description

An hypergeometric distribution describes the probability of appearance of certain attributes (like success, color, shape, etc.) in `k` random draws, without replacement, from a finit populatin containing objects (size `a`) with that attributes and objects (size `b`) without such attributes. In contrast, binomial distribution describes `k` draws from n population with replacement.

### Syntax

`result = [[statislib(module):hypergeometric_distribution_rvs(interface)]](k, a, b, [array_size])`

### Class

Elemental function

### Arguments

`k`: has `intent in` and is a scalar of type `integer`.

`a`: has `intent in` and is a scalar of type `integer`.

`b`: has `intent in` and is a scalar of type `integer`.

`array_size`: optional argument has `intent in` and is a scalar of type `integer`.

### Return value

The result is a scalar or rank one array with a size of `array_size` of type `integer`.

### Example

```fortran
program demo_hypgeo_rvs
    use stdlib_stats_distribution_PRNG, only : rand_seed
    use statislib, only: hpg_rvs => hypergeometric_distribution_rvs

    implicit none
    integer :: k(2,3,4), a(2,3,4), b(2,3,4)
    integer :: put, get

    put = 12345678
    call random_seed(put, get)

    k(:,:,:) = 25; a(:,:,:) = 30; b(:,:,:) = 40
    ! 25 draws from a sample size of (30+40) in which 30 samples with desired
    ! attributes

    print *, hpg_rvs(k, a, b)       ! a rank 3 array of 60 hypgeo random variate

!11          14          12          10          10          13
!12          13          11          10          11           8
!15          13          11          11          12           9
!13          10          12          14          10           8

    k(1,1,1) = 40; a(1,1,1) = 35; b(1,1,1) = 30
    ! 40 draws from a sample size of (35 + 30) in which 35 sampes with the
    ! desired attributes

    print *, hpg_rvs(k(1,1,1), a(1,1,1), b(1,1,1), 10)
    ! an array of 10 hypergeometric random variates

!21          25          23          20          22          21
!22          21          19          24

end program demo_hypgeo_rvs
```

## `hypergeometric_distribution_pmf` - Hypergeometric probability mass function

### Status

Release

### Description

The probability mass function of the hypergeometric distribution.

f(x; k, a, b) = C(a, x) C(b, k - x) / C(a + b, k)

where C is the binomial coefficient

### Syntax

`result = [[statislib(module):hypergeometric_distribution_pmf(interface)]](x, k, a, b)`

### Class

Elemental function

### Arguments

`x`: has `intent in` and is a scalar of type `integer`.

`k`: has `intent in` and is a scalar of type `integer`.

`a`: has `intent in` and is a scalar of type `integer`.

`b`: has `intent in` and is a scalar of type `integer`.

All arguments must have the same type and kind.

### Return value

The result is a scalar of type `real` with a shape conformable to auguments.

### Example

```fortran
program demo_hypgeo_pmf
    use stdlib_stats_distribution_PRNG, only : rand_seed
    use statislib, only: hpg_pmf => hypergeometric_distribution_pmf, &
                         hpg_rvs => hypergeometric_distribution_rvs

    implicit none
    integer :: k(2,3,4), a(2,3,4), b(2,3,4)
    integer :: put, get

    put = 12345678
    call random_seed(put, get)

    k(:,:,:) = 25; a(:,:,:) = 35; b(:,:,:) = 45

    print *, hpg_pmf(10, k(1,1,1), a(1,1,1), b(1,1,1))
    ! a probability density for random variate value of 10 in hypgeo with 25
    ! draws from a population with 35 objects with the desired attributes and 45
    ! objects without.

!0.174210683

    print *, hpg_pmf(hpg_rvs(k, a, b), k, a, b)
    ! a rank 3 array of hypgeo probability density

!0.191580668       6.47993386E-02  0.167633072      0.174210683
!0.174210683      0.116835177      0.167633072      0.116835177
!0.191580668      0.174210683      0.191580668      0.125632703
!2.85117105E-02  0.116835177      0.167633072      0.191580668
!0.167633072      0.125632703      0.116835177      0.174210683
!0.116835177       2.85117105E-02  0.191580668       7.14382008E-02

end program demo_hypgeo_pmf
```

## `hypergeometric_distribution_cdf` - Hypergeometric cumulative distribution function

### Status

Release

### Description

The cumuative distribution function of the hypergeometric distribution.

F(x; k, a, b) = &sum;<sup>x</sup>C(a, i) C(b, k - i) / C(a + b, k);  i = 0, 1, &hellip;, x

where C is the binomial coefficient

### Syntax

`result = [[statislib(module):hypergeometric_distribution_cdf(interface)]](x, k, a, b)`

### Class

Elemental function

### Arguments

`x`: has `intent in` and is a scalar of type `integer`.

`k`: has `intent in` and is a scalar of type `integer`.

`a`: has `intent in` and is a scalar of type `integer`.

`b`: has `intent in` and is a scalar of type `integer`.

All arguments must have the same type and kind.

### Return value

The result is a scalar of type `real` with a shape conformable to auguments.

### Example

```fortran
program demo_hypgeo_cdf
    use stdlib_stats_distribution_PRNG, only : rand_seed
    use statislib, only: hpg_cdf => hypergeometric_distribution_cdf,     &
                         hpg_rvs => hypergeometric_distribution_rvs

    implicit none
    integer :: k(2,3,4), a(2,3,4), b(2,3,4)
    integer :: put, get

    put = 1234567
    call random_seed(put, get)

    k(:,:,:) = 1000; a(:,:,:) = 800; b(:,:,:) = 800

    print *, hpg_cdf(500, k(1,1,1),a(1,1,1),b(1,1,1))
    ! total probability for k not greater than 700

!0.520590842

    print *, hpg_cdf(hpg_rvs(k, a, b), k, a, b)
    ! a rank 3 array of hypgeo cumulative probability

!0.479409307      0.163258627      0.714986622      0.358878732
!0.358878732       5.46886921E-02  0.251017034      0.163258627
!0.561553657      0.479409307      0.117466621      0.321060985
!0.285013378       2.19860170E-02  0.479409307      0.641121387
!0.321060985      0.955840051      0.945311368      0.978014112
!0.398132980       6.71087205E-02  0.601867199      0.748982966

end program demo_hypgeo_cdf
```
