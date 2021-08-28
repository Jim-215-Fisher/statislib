---

title: Statislib
---

# Statistical Distribution -- Beta Distribution

[TOC]

## `beta_distribution_rvs` - beta distribution random variates

### Status

Release

### Description

Beta distribution is used to model random variables in a finit range, also known as the beta distribution of the first kind.

With two auguments for shape parameters &alpha;>0, &beta;>0, the function returns a beta distributed random variate Beta(&alpha;,&beta;).

With three auguments, the function returns a rank one array of beta distributed random variate Beta(&alpha;, &beta;).

For `complex` auguments, the real and imaginary parts are independent of each other.

### Syntax

`result = [[statislib(module):beta_distribution_rvs(interface)]](alpha, beta [, array_size])`

### Class

Elemental function

### Arguments

`alpha` : has `intent(in)` ans is a scalar of type `real` or `complex`.

`beta`:  has `intent(in)` and is a scalar of type `real` or `complex`.

`array_size`: optional argument has `intent(in)` and is a scalar of type `integer`.

### Return value

The result is a scalar or rank one array, with a size of `array_size`, of type `real` or `complex`.

### Example

```fortran
program demo_beta_rvs
    use stdlib_stats_distribution_PRNG, only : random_seed
    use statislib, only: rbeta => beta_distribution_rvs

    implicit none
    real ::  aa(2,3,4), bb(2,3,4)
    complex :: a, b
    integer :: put, get

    put = 1234567
    call random_seed(put, get)

    print *, rbeta(0.5, 2.0)
    !single standard beta random variate with shape alpha=0.5, beta=2.0

! 3.38860527E-02

    print *, rbeta(3.0,2.0)  !beta random variate with alpha=3.0, beta=2.0

! 0.570277154

    aa(:,:,:) = 0.8; bb(:,:,:)=0.6
    print *, rbeta(aa, bb)
    !a rank 3 array of 24 beta random variates with alpha=0.8, beta=0.6

!  0.251766384  0.578202426  0.539138556  0.210130826  0.908130825 
!  0.880996943  9.49194748E-03  0.945992589  0.290732056  0.803920329 
!  7.64303207E-02  0.943150401  0.927998245  0.831781328  0.671169102 
!  0.983966410  0.289062619  0.801237404  0.891931713  0.897902310 
!  0.845606744  1.50359496E-02  0.913162351  0.915781260

    print *, rbeta(1.2,0.7,10)
    ! an array of 10 random variates with alpha=1.2, beta=0.7

!  3.69704030E-02  0.748639643  0.896924615  0.350249082  0.813054740 
!  0.448072791  9.39368978E-02  0.475665420  0.567661405  0.893835664

    a = (3.0, 4.0)
    b = (2.0, 0.7)
    print *, rbeta(a, b)
    !single complex beta random variate with real part of alpha = 3.0, beta=2.0;
    !imagainary part of alpha=4.0, beta=0.7

! (0.730923295,0.878909111)

end program demo_beta_rvs
```

## `beta_distribution_pdf` - beta probability density function

### Status

Release

### Description

The probability density function of the continuous beta distribution.

f(x) = x<sup>&alpha;-1</sup>(1-x)<sup>&beta;-1</sup> / B(&alpha;,&beta;)

where B(&alpha;,&beta;) = &Gamma;(&alpha;)&Gamma;(&beta;) / &Gamma;(&alpha; + &beta;)

x is supported in [0, 1]

### Syntax

`result = [[statislib(module):beta_distribution_pdf(interface)]](x, alpha, beta)`

### Class

Elemental function

### Arguments

`x`: has `intent(in)` and is a scalar of type `real` or `complex`.

`alpha` has `intent(in)` and is a scalar of type real` or `complex`.

`beta`: has `intent(in)` and is a scalar of type `real` or `complex`.

All arguments must have the same type.

### Return value

The result is a scalar or an array, with a shape conformable to auguments, of type `real`.

### Example

```fortran
program demo_beta_pdf
    use stdlib_stats_distribution_PRNG, onyl : random_seed
    use statislib, only: rbeta => beta_distribution_rvs,&
                         beta_pdf => beta_distribution_pdf

    implicit none
    real :: x(2,3,4),aa(2,3,4),bb(2,3,4)
    complex :: a, b
    integer :: put, get

    put = 1234567
    call random_seed(put, get)

    print *, beta_pdf(0.65, 1.0, 1.0)
    !a probability density at 0.65 with alpha=1.0, beta=1.0

! 1.00000000

    aa(:,:,:) = 2.0
    bb(:,:,:) = 1.0
    x = reshape(rbeta(2.0, 1.0, 24),[2,3,4]) ! beta random variates array
    print *, beta_pdf(x,aa,bb)     ! a rank 3 beta probability density array

!  1.59333873  1.60591197  1.26951313  0.825298846  1.84633350  0.715566635 
!  0.589380264  1.71169996  1.20676148  1.79188335  0.853198409  1.60287392 
!  0.818829894  1.05774701  0.608810604  1.40428901  1.45229220  1.92566359 
!  1.81786251  1.69419682  1.60652530  1.87064970  1.78721440  0.851465702

    a = (1.0, 1.5)
    b = (1.0, 2.)
    print *, beta_pdf((0.5,0.3), a, b)
    ! a complex expon probability density function at (1.5,1.0) with real part of
    !alpha=1.0, beta=1.0 and imaginary part of alpha=1.5, beta=2.0

! 1.43777180

end program demo_beta_pdf
```

## `beta_distribution_cdf` - beta cumulative distribution function

### Status

Release

### Description

Cumulative distribution function of the beta continuous distribution

F(x)=B(x; &alpha;, &beta;) / B(&alpha;, &beta;) 

where B(x; &alpha;, &beta;) = &int;<sup>x</sup> t<sup>&alpha;-1</sup>(1-t)<sup>&beta;-1</sup> dt

x is supported in [0, 1]

### Syntax

`result = [[statislib(module):beta_distribution_cdf(interface)]](x, alpha, beta)`

### Class

Elemental function

### Arguments

`x`: has `intent(in)` and is a scalar of type `real` or `complex`.

`alpha`: has `intent(in)` and is a scalar of type `real` or `complex`.

`beta`: has `intent(in)` and is a scalar of type `real` or `complex`.

All arguments must have the same type.

### Return value

The result is a scalar of type `real` with a shape conformable to auguments.

### Example

```fortran
program demo_beta_cdf
    use stdlib_stats_distribution_PRNG, onyl : random_seed
    use statislib, only: rbeta => beta_distribution_rvs,&
                         beta_cdf => beta_distribution_cdf

    implicit none
    real :: x(2,3,4),aa(2,3,4),bb(2,3,4)
    complex :: a, b
    integer :: seed_put, seed_get

    seed_put = 1234567
    call random_seed(seed_put, seed_get)

    print *, beta_cdf(0.4, 0.5,1.0)
    ! a standard beta cumulative at 0.4 with alpha=0.5, beta=1.0

! 0.632455528

    print *, beta_cdf(0.8, 1.5,2.0)
    ! a cumulative at 0.8 with alpha=1.5, beta=2.0

! 0.930204272

    aa(:,:,:) = 2.0
    bb(:,:,:) = 3.0
    x = reshape(rbeta(2.0, 3.0, 24),[2,3,4])
    !beta random variates array with alpha=2.0, beta=3.0
    print *, beta_cdf(x,aa,bb)        ! a rank 3 standard beta cumulative array

!  0.671738267  0.630592465  0.557911158  0.157377750  0.808665335 
!  8.00214931E-02  0.118469201  0.868274391  0.268321663  0.850062668 
!  7.99631923E-02  0.756867588  0.201488778  0.228500694  0.100621507 
!  0.412083119  0.480990469  0.831927001  0.791745305  0.654913783 
!  0.528246164  0.275093734  0.757254481  0.212538764

    a = (.7, 2.1)
    b = (0.5,1.0)
    print *, beta_cdf((0.5,0.5),a,b)
    !complex beta cumulative distribution at (0.5,0.5) with real part of
    !alpha=0.7,beta=0.5 and imaginary part of alpha=2.1,beta=1.0

! 9.32183489E-02

end program demo_beta_cdf

```
