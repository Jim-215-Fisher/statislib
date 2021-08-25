---
title: statslib
---

# Statistical Distribution -- Logarithmic Distribution

[TOC]


## `logarithmic_distribution_rvs` - logarithmic distribution random variates

### Status

Release

### Description

A logarithmic distribution, also known as logarithmic series distribution or simply log-series distribution, is a discrete probability distribution derived from the MacLaurin expansion. The distribution is only defined at integer value and is supported on the set {1,2,3...}.

### Syntax

`result = [[stdlib_stats_distribution_logarithmic(module):logarithmic_distribution_rvs(interface)]](p, [array_size])`

### Class

Elemental function

### Arguments

`p`: has `intent in` and is a scalar of type `real`.

`array_size`: optional argument has `intent in` and is a scalar of type `integer`.

### Return value

The result is a scalar or rank one array with a size of `array_size` of type `integer`.

### Example

```fortran
program demo_logarithmic_rvs
    use stdlib_stats_distribution_PRNG, only : rand_seed
    use statislib_logarithmic, only: log_rvs => logarithmic_distribution_rvs

    implicit none
    real :: p(3,4,5)
    integer :: put, get

    put = 12345678
    call random_seed(put, get)

    p(:,:,:) = 0.43
    print *, log_rvs(p)       ! a rank 3 array of 60 logarithmic random variate

!1           3           1           1           1           2
!1           2           1           1           1           1
!4           2           1           1           1           1
!2           1           2           3           1           1
!1           3           2           1           1           1
!1           1           1           2           4           1
!1           1           1           1           1           1
!2           2           2           4           1           1
!1           4           2           1           1           2
!3           3           1           1           1           1

    print *, log_rvs(0.96, 10)     ! an array of 10 logarithmic random variates

!3           8           1           3           1           3
!3          30           3           2

end program demo_logarithmic_rvs
```

## `logarithmic_distribution_pmf` - Logarithmic probability mass function

### Status

Release

### Description

The probability mass function of the discrete logarithmic distribution.

f(k)=- p<sup>k</sup> / [k&middot;ln(1-p)]

### Syntax

`result = [[stdlib_stats_distribution_logarithmic(module):logarithmic_distribution_pmf(interface)]](k, p)`

### Class

Elemental function

### Arguments

`k`: has `intent in` and is a scalar of type `integer`.

`p`: has `intent in` and is a scalar of type `real`.


### Return value

The result is a scalar or array of type `real` with array shape conformable to auguments.

### Example

```fortran
program demo_log_pmf
    use stdlib_stats_distribution_PRNG, only : rand_seed
    use statislib_logarithmic, only:   log_pmf => logarithmic_distribution_ pmf,&
                                       log_rvs => logarithmic_distribution_rvs

    implicit none
    real :: p(2,3,4)
    integer :: m(2,3,4)
    integer :: put, get

    put = 12345678
    call random_seed(put, get)

    p(:,:,:) = 0.98
    m = log_rvs(p)
    print *, log_pmf(5, 0.67)    ! a probability density for 5 in logarithmic

!2.43559238E-02

    print *, log_pmf(m, p)   ! a rank 3 array of logarithmic probability density

!4.40803869E-03  0.122749761       3.27161234E-03   2.08861958E-02
!8.01965147E-02  0.250509709      0.250509709       8.01965147E-02
!1.37605444E-02   1.93856540E-03   7.44988676E-03  0.250509709
!1.15636736E-02  0.250509709       1.06657883E-02   1.67159177E-02
!0.250509709       5.89444377E-02   4.62124422E-02   5.89444377E-02
!5.89444377E-02   1.27288382E-04   4.62124422E-02  0.250509709

end program demo_log_pmf
```

## `logarithmic_distribution_cdf` - Logarithmic cumulative distribution function

### Status

Release

### Description

The cumuative distribution function of the discrete logarithmic distribution.

F(k) = 1 + B(p; k + 1, 0) / ln(1 - p)

where B is the incomplete beta fucntion

### Syntax

`result = [[stdlib_stats_distribution_logarithmic(module):logarithmic_distribution_cdf(interface)]](k, p)`

### Class

Elemental function

### Arguments

`k`: has `intent in` and is a scalar or an array of type `integer`.

`p`: has `intent in` and is a scalar of type `real`.

### Return value

The result is a scalar or array of type `real` with array shape conformable to auguments.

### Example

```fortran
program demo_log_cdf
    use stdlib_stats_distribution_PRNG, only : rand_seed
    use statislib_logarithmic, only:   log_cdf => logarithmic_distribution_cdf, &
                                       log_rvs => logarithmic_distribution_rvs

    implicit none
    real :: p(2,3,4)
    integer :: m(2,3,4)
    integer :: put, get

    put = 12345678
    call random_seed(put, get)

    p(:,:,:) = 0.99
    m = log_rvs(p)
    print *, log_cdf(5, 0.76)  ! total probability for k not greater than 5

!0.931417704

    print *, log_cdf(m, p)  ! a rank 3 array of logarithmic cumulative probability

!0.900935471      0.321388721      0.917033851      0.721470177
!0.443768948      0.214975730      0.214975730      0.443768948
!0.214975730      0.214975730      0.739792585      0.214975730
!0.910335481      0.214975730      0.391621292      0.443768948
!0.980208516      0.214975730      0.391621292      0.783866465
!0.214975730      0.519143164      0.632457078      0.632457078

end program demo_log_cdf
```
