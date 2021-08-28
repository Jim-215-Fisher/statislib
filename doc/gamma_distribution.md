---

title: Statislib
---

# Statistical Distribution -- Gamma Distribution

[TOC]

## `gamma_distribution_rvs` - gamma distribution random variates

### Status

Release

### Description

Gamma distribution is the maximum entropy probability distribution, and is frequently used to model waiting time.

With one augument for shape parameter &alpha;, the function returns a standard gamma distributed random variate &Gamma;(&alpha;) with rate parameter &beta; = 1.0.

With two auguments, the function returns a scalar gamma distributed random variate &Gamma;(&alpha;, &beta;), or an array of  standard gamma distributed random variates &Gamma;(&alpha;) with rate parameter &beta; = 1.0. 

With three auguments, the function returns an array of gamma distributed random variates &Gamma;(&alpha;, &beta;).

For `complex` auguments, the real and imaginary parts are independent of each other.

### Syntax

`result = [[statislib(module):gamma_distribution_rvs(interface)]](alpha [[, beta]] [, array_size])`

### Class

Elemental function

### Arguments

`alpha` : has `intent(in)` and is a scalar of type `real` or `complex`.

`beta`: optional argument has `intent(in)` and is a scalar of type `real` or `complex`. The default value of `beta` is 1.0 for the `real` and (1.0, 1.0i) for the `complex`.

`array_size`: optional argument has `intent(in)` and is a scalar of type `integer`.

### Return value

The result is a scalar or rank one array, with a size of `array_size`, of type `real` or `complex`.

### Example

```fortran
program demo_gamma_rvs
    use stdlib_stats_distribution_PRNG, only : random_seed
    use statislib, only: rgamma => gamma_distribution_rvs

    implicit none
    real ::  g(2,3,4)
    complex :: shape, scale
    integer :: put, get

    put = 1234567
    call random_seed(put, get)

    print *, rgamma(2.0)
    !single standard gamma random variate with alpha of 2.0, beta=1.0

! 2.50538206

    print *, rgamma(3.0,2.0)  !gamma random variate with alpha=3.0, beta=2.0

! 1.30591583

    g(:,:,:) = 0.5
    print *, rgamma(g)
    !a rank 3 array of 60 standard gamma random variates with beta=0.5

!  1.03841162  1.33044529  0.912742674  0.131288037  0.638593793 
!  1.03565669E-02  0.624804378  1.12179172  4.91380468E-02  6.69969944E-03 
!  6.67014271E-02  0.132111162  0.101102419  0.648416579  1.14922595 
!  2.29003578E-02  1.85964716E-04  1.21213868E-02  1.69112933 
!  7.30440915E-02  0.395139128  0.182758048  0.427981257  0.985665262

    print *, rgamma(0.5,10)
    ! an array of 10 random variates with alpha=0.5, beta=1.0

!  1.39297554E-04  0.296419382  0.352113068  2.80515051  3.65264394E-04 
!  0.197743446  5.54569438E-02  9.30598825E-02  1.02596343  1.85311246

    shape = (3.0, 4.0)
    scale = (2.0, 0.7)
    print *, rgamma(shape,scale)
    !single complex gamma random variate with real part of alpha = 3.0,
    !beta=2.0; imagainary part of alpha=4.0, beta=0.7

! (0.826188326,3.54749799)

end program demo_gamma_rvs
```

## `gamma_distribution_pdf` - gamma probability density function

### Status

Release

### Description

The probability density function of the continuous gamma distribution.

 f(x)= &beta;<sup>&alpha;</sup> / &Gamma;(&alpha;)  x<sup>&alpha;-1</sup>&middot; e<sup>-&beta;&middot;x</sup>

### Syntax

`result = [[statislib(module):gamma_distribution_pdf(interface)]](x, alpha, beta)`

### Class

Elemental function

### Arguments

`x`: has `intent(in)` and is a scalar of type `real` or `complex`.

`alpha` has `intent(in)` and is a scalar of type real` or `complex`.

`beta`: has `intent(in)` and is a scalar of type `real` or `complex`.

All arguments must have the same type.

### Return value

The result is a scalar of type `real`.

### Example

```fortran
program demo_gamma_pdf
    use stdlib_stats_distribution_PRNG, onyl : random_seed
    use statislib, only: rgamma => gamma_distribution_rvs,&
                         gamma_pdf => gamma_distribution_pdf

    implicit none
    real :: x(2,3,4),g(2,3,4),s(2,3,4)
    complex :: alpha, scale
    integer :: put, get

    put = 1234567
    call random_seed(put, get)

    print *, gamma_pdf(1.0, 1.0, 1.0)
    !a probability density at 1.0 with alpha=1.0, beta=1.0

! 0.367879450

    g(:,:,:) = 2.0
    s(:,:,:) = 1.0
    x = realpha(rgamma(2.0, 1.0, 24),[2,3,4]) ! gamma random variates array
    print *, gamma_pdf(x,g,s)     ! a rank 3 gamma probability density array

!  0.204550430  0.320178866  0.274986655  0.348611295  0.101865448 
!  0.102199331  0.358981341  0.223676488  0.254329354  0.356714427 
!  0.267390072  0.305148095  0.367848188  7.26194456E-02  1.49471285E-02 
!  0.246272027  0.360770017  0.339665830  0.101558588  0.358678699 
!  0.224196941  0.359253854  7.56355673E-02  0.251869917

    alpha = (1.0, 1.5)
    scale = (1.0, 2.)
    print *, gamma_pdf((1.5,1.0), alpha, scale)
    ! a complex expon probability density function at (1.5,1.0) with real part
    !of alpha=1.0, beta=1.0 and imaginary part of alpha=1.5, beta=2.0

! 9.63761061E-02

end program demo_gamma_pdf
```

## `gamma_distribution_cdf` - gamma cumulative distribution function

### Status

Release

### Description

Cumulative distribution function of the gamma continuous distribution

F(x)= &gamma;(&alpha;, &beta;&middot;x) / &Gamma;(&alpha;)

### Syntax

`result = [[statislib(module):gamma_distribution_cdf(interface)]](x, alpha, beta)`

### Class

Elemental function

### Arguments

`x`: has `intent(in)` and is a scalar of type `real` or `complex`.

`alpha`: has `intent(in)` and is a scalar of type `real` or `complex`.

`beta`: has `intent(in)` and is a scalar of type `real` or `complex`.

All arguments must have the same type.

### Return value

The result is a scalar of type `real`.

### Example

```fortran
program demo_gamma_cdf
    use stdlib_stats_distribution_PRNG, onyl : random_seed
    use statislib, only: rgamma => gamma_distribution_rvs,               &
                         gamma_cdf => gamma_distribution_cdf

    implicit none
    real :: x(2,3,4),g(2,3,4),s(2,3,4)
    complex :: alpha, scale
    integer :: seed_put, seed_get

    seed_put = 1234567
    call random_seed(seed_put, seed_get)

    print *, gamma_cdf(1.0, 0.5,1.0)
    ! a standard gamma cumulative at 1.0 with a alpha=0.5, beta=1.0

! 0.842700839

    print *, gamma_cdf(2.0, 1.5,2.0)
    ! a cumulative at 2.0 with a alpha=1.5, beta=2.0

! 0.953988254

    g(:,:,:) = 1.0
    s(:,:,:) = 1.0
    x = reshape(rgamma(1.0, 1.0, 24),[2,3,4])
    !gamma random variates array with a alpha=1.0, beta=1.0
    print *, gamma_cdf(x,g,s)        ! a rank 3 standard gamma cumulative array

!  0.710880339  0.472411335  0.578345954  0.383050948  0.870905757 
!  0.870430350  0.170215249  0.677347481  0.620089889  0.161825046 
!  4.17549349E-02  0.510665894  0.252201647  0.911497891  0.984424412 
!  0.635621786  0.177783430  0.414842933  0.871342421  0.338317066 
!  2.06879266E-02  0.335232288  0.907408893  0.624871135

    alpha = (.7, 2.1)
    scale = (0.5,1.0)
    print *, gamma_cdf((0.5,0.5),alpha,scale)
    !complex gamma cumulative distribution at (0.5,0.5) with real part of
    !alpha=0.7,beta=0.5 and imaginary part of alpha=2.1,beta=1.0

! 2.87349485E-02

end program demo_gamma_cdf

```
