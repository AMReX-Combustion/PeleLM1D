
! Element section

Elements
     H O N AR HE C
End

! Species section

Species
     H H2 O OH H2O O2 HO2 H2O2 N2 AR HE
End

! Thermo section

Thermo All
            300            1000            5000 
!
H                       H   1               G   300.000  5000.000   300.000    1
 2.50000000e+00 0.00000000e+00 0.00000000e+00 0.00000000e+00 0.00000000e+00    2
 2.54716300e+04-4.60117600e-01 2.50000000e+00 0.00000000e+00 0.00000000e+00    3
 0.00000000e+00 0.00000000e+00 2.54716300e+04-4.60117600e-01                   4
!
H2                      H   2               G   300.000  5000.000   300.000    1
 2.99142300e+00 7.00064400e-04-5.63382900e-08-9.23157800e-12 1.58275200e-15    2
-8.35034000e+02-1.35511000e+00 3.29812400e+00 8.24944200e-04-8.14301500e-07    3
-9.47543400e-11 4.13487200e-13-1.01252100e+03-3.29409400e+00                   4
!
O                       O   1               G   300.000  5000.000   300.000    1
 2.54206000e+00-2.75506200e-05-3.10280300e-09 4.55106700e-12-4.36805200e-16    2
 2.92308000e+04 4.92030800e+00 2.94642900e+00-1.63816600e-03 2.42103200e-06    3
-1.60284300e-09 3.89069600e-13 2.91476400e+04 2.96399500e+00                   4
!
OH                      O   1H   1          G   200.000  6000.000   200.000    1
 2.86472886e+00 1.05650448e-03-2.59082758e-07 3.05218674e-11-1.33195876e-15    2
 3.68362875e+03 5.70164073e+00 4.12530561e+00-3.22544939e-03 6.52764691e-06    3
-5.79853643e-09 2.06237379e-12 3.34630913e+03-6.90432960e-01                   4
!
H2O                     H   2O   1          G   300.000  5000.000   300.000    1
 2.67214600e+00 3.05629300e-03-8.73026000e-07 1.20099600e-10-6.39161800e-15    2
-2.98992100e+04 6.86281700e+00 3.38684200e+00 3.47498200e-03-6.35469600e-06    3
 6.96858100e-09-2.50658800e-12-3.02081100e+04 2.59023300e+00                   4
!
O2                      O   2               G   300.000  5000.000   300.000    1
 3.69757800e+00 6.13519700e-04-1.25884200e-07 1.77528100e-11-1.13643500e-15    2
-1.23393000e+03 3.18916600e+00 3.21293600e+00 1.12748600e-03-5.75615000e-07    3
 1.31387700e-09-8.76855400e-13-1.00524900e+03 6.03473800e+00                   4
!
HO2                     H   1O   2          G   200.000  3500.000   200.000    1
 4.01721090e+00 2.23982013e-03-6.33658150e-07 1.14246370e-10-1.07908535e-14    2
 1.11856713e+02 3.78510215e+00 4.30179801e+00-4.74912051e-03 2.11582891e-05    3
-2.42763894e-08 9.29225124e-12 2.94808040e+02 3.71666245e+00                   4
!
H2O2                    H   2O   2          G   300.000  5000.000   300.000    1
 4.57316700e+00 4.33613600e-03-1.47468900e-06 2.34890400e-10-1.43165400e-14    2
-1.80069600e+04 5.01137000e-01 3.38875400e+00 6.56922600e-03-1.48501300e-07    3
-4.62580600e-09 2.47151500e-12-1.76631500e+04 6.78536300e+00                   4
!
N2                      N   2               G   300.000  5000.000   300.000    1
 2.92664000e+00 1.48797700e-03-5.68476100e-07 1.00970400e-10-6.75335100e-15    2
-9.22797700e+02 5.98052800e+00 3.29867700e+00 1.40824000e-03-3.96322200e-06    3
 5.64151500e-09-2.44485500e-12-1.02090000e+03 3.95037200e+00                   4
!
AR                      AR  1               G   300.000  5000.000   300.000    1
 2.50000000e+00 0.00000000e+00 0.00000000e+00 0.00000000e+00 0.00000000e+00    2
-7.45375000e+02 4.36600100e+00 2.50000000e+00 0.00000000e+00 0.00000000e+00    3
 0.00000000e+00 0.00000000e+00-7.45375000e+02 4.36600100e+00                   4
!
HE                      HE  1               G   300.000  5000.000   300.000    1
 2.50000000e+00 0.00000000e+00 0.00000000e+00 0.00000000e+00 0.00000000e+00    2
-7.45375000e+02 9.15348900e-01 2.50000000e+00 0.00000000e+00 0.00000000e+00    3
 0.00000000e+00 0.00000000e+00-7.45375000e+02 9.15348800e-01                   4

End

! Reaction section

Reactions
H + O2 <=> O + OH                              1.04e+14              0          15286     !     1
O + H2 <=> H + OH                             3.818e+12              0           7948     !     2
    DUPLICATE
O + H2 <=> H + OH                             8.792e+14              0          19170     !     3
    DUPLICATE
H2 + OH <=> H2O + H                            2.16e+08           1.51           3430     !     4
OH + OH <=> O + H2O                               33400           2.42          -1930     !     5
H2 + M <=> H + H + M                          4.577e+19           -1.4         104380     !     6
    H2 / 3.50 / H2O / 13.00 / AR / 1.00 / HE / 1.00 / 
H2 + AR <=> H + H + AR                         5.84e+18           -1.1         104380     !     7
H2 + HE <=> H + H + HE                         5.84e+18           -1.1         104380     !     8
O + O + M <=> O2 + M                          6.165e+15           -0.5              0     !     9
    H2 / 3.50 / H2O / 13.00 / AR / 1.00 / HE / 1.00 / 
O + O + AR <=> O2 + AR                        1.886e+13              0          -1788     !    10
O + O + HE <=> O2 + HE                        1.886e+13              0          -1788     !    11
O + H + M <=> OH + M                          4.714e+18             -1              0     !    12
    H2 / 3.50 / H2O / 13.00 / AR / 1.75 / HE / 1.75 / 
H2O + M <=> H + OH + M                        6.064e+27         -3.322         120790     !    13
    H2 / 4.00 / H2O / 1.00 / HE / 2.10 / N2 / 3.00 / O2 / 2.50 / 
H2O + H2O <=> H + OH + H2O                    1.006e+26          -2.44         120180     !    14
H + O2 (+M) <=> HO2 (+M)                    4.65084e+12           0.44              0     !    15
    H2 / 3.00 / H2O / 15.00 / O2 / 1.78 / AR / 1.67 / HE / 1.80 / 
    LOW /      6.366e+20           -1.72           524.8 /
    TROE /            0.5           1e-30           1e+30 /
HO2 + H <=> H2 + O2                             2750000           2.09          -1451     !    16
HO2 + H <=> OH + OH                           7.079e+13              0            295     !    17
HO2 + O <=> O2 + OH                            2.85e+10              1        -723.93     !    18
HO2 + OH <=> H2O + O2                          2.89e+13              0           -497     !    19
HO2 + HO2 <=> H2O2 + O2                         4.2e+14              0          11982     !    20
    DUPLICATE
HO2 + HO2 <=> H2O2 + O2                         1.3e+11              0        -1629.3     !    21
    DUPLICATE
H2O2 (+M) <=> OH + OH (+M)                        2e+12            0.9          48749     !    22
    H2O / 8.50 / N2 / 2.50 / O2 / 2.20 / HE / 1.65 / H2O2 / 8.70 / H2 / 4.70 / 
    LOW /       2.49e+24            -2.3           48749 /
    TROE /           0.43           1e-30           1e+30 /
H2O2 + H <=> H2O + OH                          2.41e+13              0           3970     !    23
H2O2 + H <=> HO2 + H2                          4.82e+13              0           7950     !    24
H2O2 + O <=> OH + HO2                           9550000              2           3970     !    25
H2O2 + OH <=> HO2 + H2O                        1.74e+12              0            318     !    26
    DUPLICATE
H2O2 + OH <=> HO2 + H2O                        7.59e+13              0           7270     !    27
    DUPLICATE
HO2 + H <=> O + H2O                            3.97e+12              0            671     !    28
O + OH + M <=> HO2 + M                            8e+15              0              0     !    29
    H2 / 3.00 / H2O / 13.00 / AR / 1.70 / HE / 1.70 / 

End

! version
! $Id$

! Generated automatically by ChemkinPickler on Fri Aug 19 14:40:05 2016

! End of file 