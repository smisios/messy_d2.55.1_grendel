SUBROUTINE fubrad_ini_param_dyn(status)
  !
  ! Purpose: Initialize the indices of the Chappuis fluxes in input data array.
  ! --------
  ! Author: Markus Kunze, FU-Berlin, January, 2012.
  ! -------
  INTEGER, INTENT(out) :: status
  !
  ! LOCAL
  INTEGER,               PARAMETER :: nherz_def   = 15
  INTEGER,               PARAMETER :: nhart_def   = 10
  INTEGER,               PARAMETER :: nhug_def    = 18
  !
  ! ntr_km  -  maximum truncation for Koppers and Murtagh (1996)
  ! ntr_k   -  size of parameter arrays for Kockarts (1994)
  INTEGER,               PARAMETER :: ntr_km = 20
  INTEGER,               PARAMETER :: ntr_k  = 12
  CHARACTER(LEN=*),      PARAMETER :: substr = 'fubrad_ini_param'
  status = 1 ! ERROR
  nherz = nherz_def
  nhart = nhart_def
  nhug  = nhug_def
  !
  SELECT CASE(nbands)
  CASE( 10)
    nherz =   1
    nhart =   1
    nhug  =   1
    nchap =   1
    IDX_SRC  = (/  3,  5 /)
    IDX_SRB  = (/  6,  6 /)
    IDX_HRZ  = (/  7,  7 /)
    IDX_HAR  = (/  8,  8 /)
    IDX_HUG  = (/  9,  9 /)
    IDX_CHA  = (/ 11, 11 /)
    IDX_ADD  =  10
  CASE( 14)
    nherz =   2
    nhart =   2
    nhug  =   2
    nchap =   2
    IDX_SRC  = (/  3,  5 /)
    IDX_SRB  = (/  6,  6 /)
    IDX_HRZ  = (/  7,  8 /)
    IDX_HAR  = (/  9, 10 /)
    IDX_HUG  = (/ 11, 12 /)
    IDX_CHA  = (/ 14, 15 /)
    IDX_ADD  =  13
  CASE( 18)
    nherz =   3
    nhart =   3
    nhug  =   3
    nchap =   3
    IDX_SRC  = (/  3,  5 /)
    IDX_SRB  = (/  6,  6 /)
    IDX_HRZ  = (/  7,  9 /)
    IDX_HAR  = (/ 10, 12 /)
    IDX_HUG  = (/ 13, 15 /)
    IDX_CHA  = (/ 17, 19 /)
    IDX_ADD  =  16
  CASE( 22)
    nherz =   4
    nhart =   4
    nhug  =   4
    nchap =   4
    IDX_SRC  = (/  3,  5 /)
    IDX_SRB  = (/  6,  6 /)
    IDX_HRZ  = (/  7, 10 /)
    IDX_HAR  = (/ 11, 14 /)
    IDX_HUG  = (/ 15, 18 /)
    IDX_CHA  = (/ 20, 23 /)
    IDX_ADD  =  19
  CASE( 26)
    nherz =   5
    nhart =   5
    nhug  =   5
    nchap =   5
    IDX_SRC  = (/  3,  5 /)
    IDX_SRB  = (/  6,  6 /)
    IDX_HRZ  = (/  7, 11 /)
    IDX_HAR  = (/ 12, 16 /)
    IDX_HUG  = (/ 17, 21 /)
    IDX_CHA  = (/ 23, 27 /)
    IDX_ADD  =  22
  CASE( 30)
    nherz =   6
    nhart =   6
    nhug  =   6
    nchap =   6
    IDX_SRC  = (/  3,  5 /)
    IDX_SRB  = (/  6,  6 /)
    IDX_HRZ  = (/  7, 12 /)
    IDX_HAR  = (/ 13, 18 /)
    IDX_HUG  = (/ 19, 24 /)
    IDX_CHA  = (/ 26, 31 /)
    IDX_ADD  =  25
  CASE( 34)
    nherz =   7
    nhart =   7
    nhug  =   7
    nchap =   7
    IDX_SRC  = (/  3,  5 /)
    IDX_SRB  = (/  6,  6 /)
    IDX_HRZ  = (/  7, 13 /)
    IDX_HAR  = (/ 14, 20 /)
    IDX_HUG  = (/ 21, 27 /)
    IDX_CHA  = (/ 29, 35 /)
    IDX_ADD  =  28
  CASE( 38)
    nherz =   8
    nhart =   8
    nhug  =   8
    nchap =   8
    IDX_SRC  = (/  3,  5 /)
    IDX_SRB  = (/  6,  6 /)
    IDX_HRZ  = (/  7, 14 /)
    IDX_HAR  = (/ 15, 22 /)
    IDX_HUG  = (/ 23, 30 /)
    IDX_CHA  = (/ 32, 39 /)
    IDX_ADD  =  31
  CASE( 42)
    nherz =   9
    nhart =   9
    nhug  =   9
    nchap =   9
    IDX_SRC  = (/  3,  5 /)
    IDX_SRB  = (/  6,  6 /)
    IDX_HRZ  = (/  7, 15 /)
    IDX_HAR  = (/ 16, 24 /)
    IDX_HUG  = (/ 25, 33 /)
    IDX_CHA  = (/ 35, 43 /)
    IDX_ADD  =  34
  CASE( 46)
    nherz =  10
    nhart =  10
    nhug  =  10
    nchap =  10
    IDX_SRC  = (/  3,  5 /)
    IDX_SRB  = (/  6,  6 /)
    IDX_HRZ  = (/  7, 16 /)
    IDX_HAR  = (/ 17, 26 /)
    IDX_HUG  = (/ 27, 36 /)
    IDX_CHA  = (/ 38, 47 /)
    IDX_ADD  =  37
  CASE( 50)
    nherz =  11
    nhart =  11
    nhug  =  11
    nchap =  11
    IDX_SRC  = (/  3,  5 /)
    IDX_SRB  = (/  6,  6 /)
    IDX_HRZ  = (/  7, 17 /)
    IDX_HAR  = (/ 18, 28 /)
    IDX_HUG  = (/ 29, 39 /)
    IDX_CHA  = (/ 41, 51 /)
    IDX_ADD  =  40
  CASE( 54)
    nherz =  12
    nhart =  12
    nhug  =  12
    nchap =  12
    IDX_SRC  = (/  3,  5 /)
    IDX_SRB  = (/  6,  6 /)
    IDX_HRZ  = (/  7, 18 /)
    IDX_HAR  = (/ 19, 30 /)
    IDX_HUG  = (/ 31, 42 /)
    IDX_CHA  = (/ 44, 55 /)
    IDX_ADD  =  43
  CASE( 58)
    nherz =  13
    nhart =  13
    nhug  =  13
    nchap =  13
    IDX_SRC  = (/  3,  5 /)
    IDX_SRB  = (/  6,  6 /)
    IDX_HRZ  = (/  7, 19 /)
    IDX_HAR  = (/ 20, 32 /)
    IDX_HUG  = (/ 33, 45 /)
    IDX_CHA  = (/ 47, 59 /)
    IDX_ADD  =  46
  CASE( 62)
    nherz =  14
    nhart =  14
    nhug  =  14
    nchap =  14
    IDX_SRC  = (/  3,  5 /)
    IDX_SRB  = (/  6,  6 /)
    IDX_HRZ  = (/  7, 20 /)
    IDX_HAR  = (/ 21, 34 /)
    IDX_HUG  = (/ 35, 48 /)
    IDX_CHA  = (/ 50, 63 /)
    IDX_ADD  =  49
  CASE( 66)
    nherz =  15
    nhart =  15
    nhug  =  15
    nchap =  15
    IDX_SRC  = (/  3,  5 /)
    IDX_SRB  = (/  6,  6 /)
    IDX_HRZ  = (/  7, 21 /)
    IDX_HAR  = (/ 22, 36 /)
    IDX_HUG  = (/ 37, 51 /)
    IDX_CHA  = (/ 53, 67 /)
    IDX_ADD  =  52
  CASE( 70)
    nherz =  16
    nhart =  16
    nhug  =  16
    nchap =  16
    IDX_SRC  = (/  3,  5 /)
    IDX_SRB  = (/  6,  6 /)
    IDX_HRZ  = (/  7, 22 /)
    IDX_HAR  = (/ 23, 38 /)
    IDX_HUG  = (/ 39, 54 /)
    IDX_CHA  = (/ 56, 71 /)
    IDX_ADD  =  55
  CASE( 74)
    nherz =  17
    nhart =  17
    nhug  =  17
    nchap =  17
    IDX_SRC  = (/  3,  5 /)
    IDX_SRB  = (/  6,  6 /)
    IDX_HRZ  = (/  7, 23 /)
    IDX_HAR  = (/ 24, 40 /)
    IDX_HUG  = (/ 41, 57 /)
    IDX_CHA  = (/ 59, 75 /)
    IDX_ADD  =  58
  CASE( 78)
    nherz =  18
    nhart =  18
    nhug  =  18
    nchap =  18
    IDX_SRC  = (/  3,  5 /)
    IDX_SRB  = (/  6,  6 /)
    IDX_HRZ  = (/  7, 24 /)
    IDX_HAR  = (/ 25, 42 /)
    IDX_HUG  = (/ 43, 60 /)
    IDX_CHA  = (/ 62, 79 /)
    IDX_ADD  =  61
  CASE( 82)
    nherz =  19
    nhart =  19
    nhug  =  19
    nchap =  19
    IDX_SRC  = (/  3,  5 /)
    IDX_SRB  = (/  6,  6 /)
    IDX_HRZ  = (/  7, 25 /)
    IDX_HAR  = (/ 26, 44 /)
    IDX_HUG  = (/ 45, 63 /)
    IDX_CHA  = (/ 65, 83 /)
    IDX_ADD  =  64
  CASE( 86)
    nherz =  20
    nhart =  20
    nhug  =  20
    nchap =  20
    IDX_SRC  = (/  3,  5 /)
    IDX_SRB  = (/  6,  6 /)
    IDX_HRZ  = (/  7, 26 /)
    IDX_HAR  = (/ 27, 46 /)
    IDX_HUG  = (/ 47, 66 /)
    IDX_CHA  = (/ 68, 87 /)
    IDX_ADD  =  67
  CASE( 90)
    nherz =  21
    nhart =  21
    nhug  =  21
    nchap =  21
    IDX_SRC  = (/  3,  5 /)
    IDX_SRB  = (/  6,  6 /)
    IDX_HRZ  = (/  7, 27 /)
    IDX_HAR  = (/ 28, 48 /)
    IDX_HUG  = (/ 49, 69 /)
    IDX_CHA  = (/ 71, 91 /)
    IDX_ADD  =  70
  CASE( 94)
    nherz =  22
    nhart =  22
    nhug  =  22
    nchap =  22
    IDX_SRC  = (/  3,  5 /)
    IDX_SRB  = (/  6,  6 /)
    IDX_HRZ  = (/  7, 28 /)
    IDX_HAR  = (/ 29, 50 /)
    IDX_HUG  = (/ 51, 72 /)
    IDX_CHA  = (/ 74, 95 /)
    IDX_ADD  =  73
  CASE( 98)
    nherz =  23
    nhart =  23
    nhug  =  23
    nchap =  23
    IDX_SRC  = (/  3,  5 /)
    IDX_SRB  = (/  6,  6 /)
    IDX_HRZ  = (/  7, 29 /)
    IDX_HAR  = (/ 30, 52 /)
    IDX_HUG  = (/ 53, 75 /)
    IDX_CHA  = (/ 77, 99 /)
    IDX_ADD  =  76
  CASE(102)
    nherz =  24
    nhart =  24
    nhug  =  24
    nchap =  24
    IDX_SRC  = (/  3,  5 /)
    IDX_SRB  = (/  6,  6 /)
    IDX_HRZ  = (/  7, 30 /)
    IDX_HAR  = (/ 31, 54 /)
    IDX_HUG  = (/ 55, 78 /)
    IDX_CHA  = (/ 80,103 /)
    IDX_ADD  =  79
  CASE(106)
    nherz =  25
    nhart =  25
    nhug  =  25
    nchap =  25
    IDX_SRC  = (/  3,  5 /)
    IDX_SRB  = (/  6,  6 /)
    IDX_HRZ  = (/  7, 31 /)
    IDX_HAR  = (/ 32, 56 /)
    IDX_HUG  = (/ 57, 81 /)
    IDX_CHA  = (/ 83,107 /)
    IDX_ADD  =  82
  CASE(110)
    nherz =  26
    nhart =  26
    nhug  =  26
    nchap =  26
    IDX_SRC  = (/  3,  5 /)
    IDX_SRB  = (/  6,  6 /)
    IDX_HRZ  = (/  7, 32 /)
    IDX_HAR  = (/ 33, 58 /)
    IDX_HUG  = (/ 59, 84 /)
    IDX_CHA  = (/ 86,111 /)
    IDX_ADD  =  85
  CASE(114)
    nherz =  27
    nhart =  27
    nhug  =  27
    nchap =  27
    IDX_SRC  = (/  3,  5 /)
    IDX_SRB  = (/  6,  6 /)
    IDX_HRZ  = (/  7, 33 /)
    IDX_HAR  = (/ 34, 60 /)
    IDX_HUG  = (/ 61, 87 /)
    IDX_CHA  = (/ 89,115 /)
    IDX_ADD  =  88
  CASE(118)
    nherz =  28
    nhart =  28
    nhug  =  28
    nchap =  28
    IDX_SRC  = (/  3,  5 /)
    IDX_SRB  = (/  6,  6 /)
    IDX_HRZ  = (/  7, 34 /)
    IDX_HAR  = (/ 35, 62 /)
    IDX_HUG  = (/ 63, 90 /)
    IDX_CHA  = (/ 92,119 /)
    IDX_ADD  =  91
  CASE(122)
    nherz =  29
    nhart =  29
    nhug  =  29
    nchap =  29
    IDX_SRC  = (/  3,  5 /)
    IDX_SRB  = (/  6,  6 /)
    IDX_HRZ  = (/  7, 35 /)
    IDX_HAR  = (/ 36, 64 /)
    IDX_HUG  = (/ 65, 93 /)
    IDX_CHA  = (/ 95,123 /)
    IDX_ADD  =  94
  CASE(126)
    nherz =  30
    nhart =  30
    nhug  =  30
    nchap =  30
    IDX_SRC  = (/  3,  5 /)
    IDX_SRB  = (/  6,  6 /)
    IDX_HRZ  = (/  7, 36 /)
    IDX_HAR  = (/ 37, 66 /)
    IDX_HUG  = (/ 67, 96 /)
    IDX_CHA  = (/ 98,127 /)
    IDX_ADD  =  97
  CASE(130)
    nherz =  31
    nhart =  31
    nhug  =  31
    nchap =  31
    IDX_SRC  = (/  3,  5 /)
    IDX_SRB  = (/  6,  6 /)
    IDX_HRZ  = (/  7, 37 /)
    IDX_HAR  = (/ 38, 68 /)
    IDX_HUG  = (/ 69, 99 /)
    IDX_CHA  = (/101,131 /)
    IDX_ADD  = 100
  CASE(134)
    nherz =  32
    nhart =  32
    nhug  =  32
    nchap =  32
    IDX_SRC  = (/  3,  5 /)
    IDX_SRB  = (/  6,  6 /)
    IDX_HRZ  = (/  7, 38 /)
    IDX_HAR  = (/ 39, 70 /)
    IDX_HUG  = (/ 71,102 /)
    IDX_CHA  = (/104,135 /)
    IDX_ADD  = 103
  CASE(138)
    nherz =  33
    nhart =  33
    nhug  =  33
    nchap =  33
    IDX_SRC  = (/  3,  5 /)
    IDX_SRB  = (/  6,  6 /)
    IDX_HRZ  = (/  7, 39 /)
    IDX_HAR  = (/ 40, 72 /)
    IDX_HUG  = (/ 73,105 /)
    IDX_CHA  = (/107,139 /)
    IDX_ADD  = 106
  CASE(142)
    nherz =  34
    nhart =  34
    nhug  =  34
    nchap =  34
    IDX_SRC  = (/  3,  5 /)
    IDX_SRB  = (/  6,  6 /)
    IDX_HRZ  = (/  7, 40 /)
    IDX_HAR  = (/ 41, 74 /)
    IDX_HUG  = (/ 75,108 /)
    IDX_CHA  = (/110,143 /)
    IDX_ADD  = 109
  CASE(145)
    nherz =  35
    nhart =  34
    nhug  =  35
    nchap =  35
    IDX_SRC  = (/  3,  5 /)
    IDX_SRB  = (/  6,  6 /)
    IDX_HRZ  = (/  7, 41 /)
    IDX_HAR  = (/ 42, 75 /)
    IDX_HUG  = (/ 76,110 /)
    IDX_CHA  = (/112,146 /)
    IDX_ADD  = 111
  CASE(148)
    nherz =  36
    nhart =  34
    nhug  =  36
    nchap =  36
    IDX_SRC  = (/  3,  5 /)
    IDX_SRB  = (/  6,  6 /)
    IDX_HRZ  = (/  7, 42 /)
    IDX_HAR  = (/ 43, 76 /)
    IDX_HUG  = (/ 77,112 /)
    IDX_CHA  = (/114,149 /)
    IDX_ADD  = 113
  CASE(151)
    nherz =  37
    nhart =  34
    nhug  =  37
    nchap =  37
    IDX_SRC  = (/  3,  5 /)
    IDX_SRB  = (/  6,  6 /)
    IDX_HRZ  = (/  7, 43 /)
    IDX_HAR  = (/ 44, 77 /)
    IDX_HUG  = (/ 78,114 /)
    IDX_CHA  = (/116,152 /)
    IDX_ADD  = 115
  CASE(153)
    nherz =  37
    nhart =  34
    nhug  =  38
    nchap =  38
    IDX_SRC  = (/  3,  5 /)
    IDX_SRB  = (/  6,  6 /)
    IDX_HRZ  = (/  7, 43 /)
    IDX_HAR  = (/ 44, 77 /)
    IDX_HUG  = (/ 78,115 /)
    IDX_CHA  = (/117,154 /)
    IDX_ADD  = 116
  CASE(155)
    nherz =  37
    nhart =  34
    nhug  =  39
    nchap =  39
    IDX_SRC  = (/  3,  5 /)
    IDX_SRB  = (/  6,  6 /)
    IDX_HRZ  = (/  7, 43 /)
    IDX_HAR  = (/ 44, 77 /)
    IDX_HUG  = (/ 78,116 /)
    IDX_CHA  = (/118,156 /)
    IDX_ADD  = 117
  CASE(157)
    nherz =  37
    nhart =  34
    nhug  =  40
    nchap =  40
    IDX_SRC  = (/  3,  5 /)
    IDX_SRB  = (/  6,  6 /)
    IDX_HRZ  = (/  7, 43 /)
    IDX_HAR  = (/ 44, 77 /)
    IDX_HUG  = (/ 78,117 /)
    IDX_CHA  = (/119,158 /)
    IDX_ADD  = 118
  CASE(159)
    nherz =  37
    nhart =  34
    nhug  =  41
    nchap =  41
    IDX_SRC  = (/  3,  5 /)
    IDX_SRB  = (/  6,  6 /)
    IDX_HRZ  = (/  7, 43 /)
    IDX_HAR  = (/ 44, 77 /)
    IDX_HUG  = (/ 78,118 /)
    IDX_CHA  = (/120,160 /)
    IDX_ADD  = 119
  CASE(161)
    nherz =  37
    nhart =  34
    nhug  =  42
    nchap =  42
    IDX_SRC  = (/  3,  5 /)
    IDX_SRB  = (/  6,  6 /)
    IDX_HRZ  = (/  7, 43 /)
    IDX_HAR  = (/ 44, 77 /)
    IDX_HUG  = (/ 78,119 /)
    IDX_CHA  = (/121,162 /)
    IDX_ADD  = 120
  CASE(163)
    nherz =  37
    nhart =  34
    nhug  =  43
    nchap =  43
    IDX_SRC  = (/  3,  5 /)
    IDX_SRB  = (/  6,  6 /)
    IDX_HRZ  = (/  7, 43 /)
    IDX_HAR  = (/ 44, 77 /)
    IDX_HUG  = (/ 78,120 /)
    IDX_CHA  = (/122,164 /)
    IDX_ADD  = 121
  CASE(165)
    nherz =  37
    nhart =  34
    nhug  =  44
    nchap =  44
    IDX_SRC  = (/  3,  5 /)
    IDX_SRB  = (/  6,  6 /)
    IDX_HRZ  = (/  7, 43 /)
    IDX_HAR  = (/ 44, 77 /)
    IDX_HUG  = (/ 78,121 /)
    IDX_CHA  = (/123,166 /)
    IDX_ADD  = 122
  CASE(167)
    nherz =  37
    nhart =  34
    nhug  =  45
    nchap =  45
    IDX_SRC  = (/  3,  5 /)
    IDX_SRB  = (/  6,  6 /)
    IDX_HRZ  = (/  7, 43 /)
    IDX_HAR  = (/ 44, 77 /)
    IDX_HUG  = (/ 78,122 /)
    IDX_CHA  = (/124,168 /)
    IDX_ADD  = 123
  CASE(169)
    nherz =  37
    nhart =  34
    nhug  =  46
    nchap =  46
    IDX_SRC  = (/  3,  5 /)
    IDX_SRB  = (/  6,  6 /)
    IDX_HRZ  = (/  7, 43 /)
    IDX_HAR  = (/ 44, 77 /)
    IDX_HUG  = (/ 78,123 /)
    IDX_CHA  = (/125,170 /)
    IDX_ADD  = 124
  CASE(171)
    nherz =  37
    nhart =  34
    nhug  =  47
    nchap =  47
    IDX_SRC  = (/  3,  5 /)
    IDX_SRB  = (/  6,  6 /)
    IDX_HRZ  = (/  7, 43 /)
    IDX_HAR  = (/ 44, 77 /)
    IDX_HUG  = (/ 78,124 /)
    IDX_CHA  = (/126,172 /)
    IDX_ADD  = 125
  CASE(173)
    nherz =  37
    nhart =  34
    nhug  =  48
    nchap =  48
    IDX_SRC  = (/  3,  5 /)
    IDX_SRB  = (/  6,  6 /)
    IDX_HRZ  = (/  7, 43 /)
    IDX_HAR  = (/ 44, 77 /)
    IDX_HUG  = (/ 78,125 /)
    IDX_CHA  = (/127,174 /)
    IDX_ADD  = 126
  CASE(175)
    nherz =  37
    nhart =  34
    nhug  =  49
    nchap =  49
    IDX_SRC  = (/  3,  5 /)
    IDX_SRB  = (/  6,  6 /)
    IDX_HRZ  = (/  7, 43 /)
    IDX_HAR  = (/ 44, 77 /)
    IDX_HUG  = (/ 78,126 /)
    IDX_CHA  = (/128,176 /)
    IDX_ADD  = 127
  CASE(177)
    nherz =  37
    nhart =  34
    nhug  =  50
    nchap =  50
    IDX_SRC  = (/  3,  5 /)
    IDX_SRB  = (/  6,  6 /)
    IDX_HRZ  = (/  7, 43 /)
    IDX_HAR  = (/ 44, 77 /)
    IDX_HUG  = (/ 78,127 /)
    IDX_CHA  = (/129,178 /)
    IDX_ADD  = 128
  CASE(179)
    nherz =  37
    nhart =  34
    nhug  =  51
    nchap =  51
    IDX_SRC  = (/  3,  5 /)
    IDX_SRB  = (/  6,  6 /)
    IDX_HRZ  = (/  7, 43 /)
    IDX_HAR  = (/ 44, 77 /)
    IDX_HUG  = (/ 78,128 /)
    IDX_CHA  = (/130,180 /)
    IDX_ADD  = 129
  CASE(181)
    nherz =  37
    nhart =  34
    nhug  =  52
    nchap =  52
    IDX_SRC  = (/  3,  5 /)
    IDX_SRB  = (/  6,  6 /)
    IDX_HRZ  = (/  7, 43 /)
    IDX_HAR  = (/ 44, 77 /)
    IDX_HUG  = (/ 78,129 /)
    IDX_CHA  = (/131,182 /)
    IDX_ADD  = 130
  CASE(183)
    nherz =  37
    nhart =  34
    nhug  =  53
    nchap =  53
    IDX_SRC  = (/  3,  5 /)
    IDX_SRB  = (/  6,  6 /)
    IDX_HRZ  = (/  7, 43 /)
    IDX_HAR  = (/ 44, 77 /)
    IDX_HUG  = (/ 78,130 /)
    IDX_CHA  = (/132,184 /)
    IDX_ADD  = 131
  CASE(185)
    nherz =  37
    nhart =  34
    nhug  =  54
    nchap =  54
    IDX_SRC  = (/  3,  5 /)
    IDX_SRB  = (/  6,  6 /)
    IDX_HRZ  = (/  7, 43 /)
    IDX_HAR  = (/ 44, 77 /)
    IDX_HUG  = (/ 78,131 /)
    IDX_CHA  = (/133,186 /)
    IDX_ADD  = 132
  CASE(187)
    nherz =  37
    nhart =  34
    nhug  =  55
    nchap =  55
    IDX_SRC  = (/  3,  5 /)
    IDX_SRB  = (/  6,  6 /)
    IDX_HRZ  = (/  7, 43 /)
    IDX_HAR  = (/ 44, 77 /)
    IDX_HUG  = (/ 78,132 /)
    IDX_CHA  = (/134,188 /)
    IDX_ADD  = 133
  CASE(189)
    nherz =  37
    nhart =  34
    nhug  =  56
    nchap =  56
    IDX_SRC  = (/  3,  5 /)
    IDX_SRB  = (/  6,  6 /)
    IDX_HRZ  = (/  7, 43 /)
    IDX_HAR  = (/ 44, 77 /)
    IDX_HUG  = (/ 78,133 /)
    IDX_CHA  = (/135,190 /)
    IDX_ADD  = 134
  CASE(191)
    nherz =  37
    nhart =  34
    nhug  =  57
    nchap =  57
    IDX_SRC  = (/  3,  5 /)
    IDX_SRB  = (/  6,  6 /)
    IDX_HRZ  = (/  7, 43 /)
    IDX_HAR  = (/ 44, 77 /)
    IDX_HUG  = (/ 78,134 /)
    IDX_CHA  = (/136,192 /)
    IDX_ADD  = 135
  CASE(193)
    nherz =  37
    nhart =  34
    nhug  =  58
    nchap =  58
    IDX_SRC  = (/  3,  5 /)
    IDX_SRB  = (/  6,  6 /)
    IDX_HRZ  = (/  7, 43 /)
    IDX_HAR  = (/ 44, 77 /)
    IDX_HUG  = (/ 78,135 /)
    IDX_CHA  = (/137,194 /)
    IDX_ADD  = 136
  CASE(195)
    nherz =  37
    nhart =  34
    nhug  =  59
    nchap =  59
    IDX_SRC  = (/  3,  5 /)
    IDX_SRB  = (/  6,  6 /)
    IDX_HRZ  = (/  7, 43 /)
    IDX_HAR  = (/ 44, 77 /)
    IDX_HUG  = (/ 78,136 /)
    IDX_CHA  = (/138,196 /)
    IDX_ADD  = 137
  CASE(197)
    nherz =  37
    nhart =  34
    nhug  =  60
    nchap =  60
    IDX_SRC  = (/  3,  5 /)
    IDX_SRB  = (/  6,  6 /)
    IDX_HRZ  = (/  7, 43 /)
    IDX_HAR  = (/ 44, 77 /)
    IDX_HUG  = (/ 78,137 /)
    IDX_CHA  = (/139,198 /)
    IDX_ADD  = 138
  CASE(199)
    nherz =  37
    nhart =  34
    nhug  =  61
    nchap =  61
    IDX_SRC  = (/  3,  5 /)
    IDX_SRB  = (/  6,  6 /)
    IDX_HRZ  = (/  7, 43 /)
    IDX_HAR  = (/ 44, 77 /)
    IDX_HUG  = (/ 78,138 /)
    IDX_CHA  = (/140,200 /)
    IDX_ADD  = 139
  CASE(201)
    nherz =  37
    nhart =  34
    nhug  =  62
    nchap =  62
    IDX_SRC  = (/  3,  5 /)
    IDX_SRB  = (/  6,  6 /)
    IDX_HRZ  = (/  7, 43 /)
    IDX_HAR  = (/ 44, 77 /)
    IDX_HUG  = (/ 78,139 /)
    IDX_CHA  = (/141,202 /)
    IDX_ADD  = 140
  CASE(203)
    nherz =  37
    nhart =  34
    nhug  =  63
    nchap =  63
    IDX_SRC  = (/  3,  5 /)
    IDX_SRB  = (/  6,  6 /)
    IDX_HRZ  = (/  7, 43 /)
    IDX_HAR  = (/ 44, 77 /)
    IDX_HUG  = (/ 78,140 /)
    IDX_CHA  = (/142,204 /)
    IDX_ADD  = 141
  CASE(205)
    nherz =  37
    nhart =  34
    nhug  =  64
    nchap =  64
    IDX_SRC  = (/  3,  5 /)
    IDX_SRB  = (/  6,  6 /)
    IDX_HRZ  = (/  7, 43 /)
    IDX_HAR  = (/ 44, 77 /)
    IDX_HUG  = (/ 78,141 /)
    IDX_CHA  = (/143,206 /)
    IDX_ADD  = 142
  CASE(207)
    nherz =  37
    nhart =  34
    nhug  =  65
    nchap =  65
    IDX_SRC  = (/  3,  5 /)
    IDX_SRB  = (/  6,  6 /)
    IDX_HRZ  = (/  7, 43 /)
    IDX_HAR  = (/ 44, 77 /)
    IDX_HUG  = (/ 78,142 /)
    IDX_CHA  = (/144,208 /)
    IDX_ADD  = 143
  CASE(209)
    nherz =  37
    nhart =  34
    nhug  =  66
    nchap =  66
    IDX_SRC  = (/  3,  5 /)
    IDX_SRB  = (/  6,  6 /)
    IDX_HRZ  = (/  7, 43 /)
    IDX_HAR  = (/ 44, 77 /)
    IDX_HUG  = (/ 78,143 /)
    IDX_CHA  = (/145,210 /)
    IDX_ADD  = 144
  CASE(211)
    nherz =  37
    nhart =  34
    nhug  =  67
    nchap =  67
    IDX_SRC  = (/  3,  5 /)
    IDX_SRB  = (/  6,  6 /)
    IDX_HRZ  = (/  7, 43 /)
    IDX_HAR  = (/ 44, 77 /)
    IDX_HUG  = (/ 78,144 /)
    IDX_CHA  = (/146,212 /)
    IDX_ADD  = 145
  CASE(213)
    nherz =  37
    nhart =  34
    nhug  =  68
    nchap =  68
    IDX_SRC  = (/  3,  5 /)
    IDX_SRB  = (/  6,  6 /)
    IDX_HRZ  = (/  7, 43 /)
    IDX_HAR  = (/ 44, 77 /)
    IDX_HUG  = (/ 78,145 /)
    IDX_CHA  = (/147,214 /)
    IDX_ADD  = 146
  CASE(215)
    nherz =  37
    nhart =  34
    nhug  =  69
    nchap =  69
    IDX_SRC  = (/  3,  5 /)
    IDX_SRB  = (/  6,  6 /)
    IDX_HRZ  = (/  7, 43 /)
    IDX_HAR  = (/ 44, 77 /)
    IDX_HUG  = (/ 78,146 /)
    IDX_CHA  = (/148,216 /)
    IDX_ADD  = 147
  CASE(217)
    nherz =  37
    nhart =  34
    nhug  =  70
    nchap =  70
    IDX_SRC  = (/  3,  5 /)
    IDX_SRB  = (/  6,  6 /)
    IDX_HRZ  = (/  7, 43 /)
    IDX_HAR  = (/ 44, 77 /)
    IDX_HUG  = (/ 78,147 /)
    IDX_CHA  = (/149,218 /)
    IDX_ADD  = 148
  CASE(219)
    nherz =  37
    nhart =  34
    nhug  =  71
    nchap =  71
    IDX_SRC  = (/  3,  5 /)
    IDX_SRB  = (/  6,  6 /)
    IDX_HRZ  = (/  7, 43 /)
    IDX_HAR  = (/ 44, 77 /)
    IDX_HUG  = (/ 78,148 /)
    IDX_CHA  = (/150,220 /)
    IDX_ADD  = 149
  CASE(221)
    nherz =  37
    nhart =  34
    nhug  =  72
    nchap =  72
    IDX_SRC  = (/  3,  5 /)
    IDX_SRB  = (/  6,  6 /)
    IDX_HRZ  = (/  7, 43 /)
    IDX_HAR  = (/ 44, 77 /)
    IDX_HUG  = (/ 78,149 /)
    IDX_CHA  = (/151,222 /)
    IDX_ADD  = 150
  CASE(223)
    nherz =  37
    nhart =  34
    nhug  =  73
    nchap =  73
    IDX_SRC  = (/  3,  5 /)
    IDX_SRB  = (/  6,  6 /)
    IDX_HRZ  = (/  7, 43 /)
    IDX_HAR  = (/ 44, 77 /)
    IDX_HUG  = (/ 78,150 /)
    IDX_CHA  = (/152,224 /)
    IDX_ADD  = 151
  CASE(225)
    nherz =  37
    nhart =  34
    nhug  =  74
    nchap =  74
    IDX_SRC  = (/  3,  5 /)
    IDX_SRB  = (/  6,  6 /)
    IDX_HRZ  = (/  7, 43 /)
    IDX_HAR  = (/ 44, 77 /)
    IDX_HUG  = (/ 78,151 /)
    IDX_CHA  = (/153,226 /)
    IDX_ADD  = 152
  CASE(227)
    nherz =  37
    nhart =  34
    nhug  =  75
    nchap =  75
    IDX_SRC  = (/  3,  5 /)
    IDX_SRB  = (/  6,  6 /)
    IDX_HRZ  = (/  7, 43 /)
    IDX_HAR  = (/ 44, 77 /)
    IDX_HUG  = (/ 78,152 /)
    IDX_CHA  = (/154,228 /)
    IDX_ADD  = 153
  CASE(229)
    nherz =  37
    nhart =  34
    nhug  =  76
    nchap =  76
    IDX_SRC  = (/  3,  5 /)
    IDX_SRB  = (/  6,  6 /)
    IDX_HRZ  = (/  7, 43 /)
    IDX_HAR  = (/ 44, 77 /)
    IDX_HUG  = (/ 78,153 /)
    IDX_CHA  = (/155,230 /)
    IDX_ADD  = 154
  CASE(231)
    nherz =  37
    nhart =  34
    nhug  =  77
    nchap =  77
    IDX_SRC  = (/  3,  5 /)
    IDX_SRB  = (/  6,  6 /)
    IDX_HRZ  = (/  7, 43 /)
    IDX_HAR  = (/ 44, 77 /)
    IDX_HUG  = (/ 78,154 /)
    IDX_CHA  = (/156,232 /)
    IDX_ADD  = 155
  CASE(233)
    nherz =  37
    nhart =  34
    nhug  =  78
    nchap =  78
    IDX_SRC  = (/  3,  5 /)
    IDX_SRB  = (/  6,  6 /)
    IDX_HRZ  = (/  7, 43 /)
    IDX_HAR  = (/ 44, 77 /)
    IDX_HUG  = (/ 78,155 /)
    IDX_CHA  = (/157,234 /)
    IDX_ADD  = 156
  CASE(235)
    nherz =  37
    nhart =  34
    nhug  =  79
    nchap =  79
    IDX_SRC  = (/  3,  5 /)
    IDX_SRB  = (/  6,  6 /)
    IDX_HRZ  = (/  7, 43 /)
    IDX_HAR  = (/ 44, 77 /)
    IDX_HUG  = (/ 78,156 /)
    IDX_CHA  = (/158,236 /)
    IDX_ADD  = 157
  CASE(237)
    nherz =  37
    nhart =  34
    nhug  =  80
    nchap =  80
    IDX_SRC  = (/  3,  5 /)
    IDX_SRB  = (/  6,  6 /)
    IDX_HRZ  = (/  7, 43 /)
    IDX_HAR  = (/ 44, 77 /)
    IDX_HUG  = (/ 78,157 /)
    IDX_CHA  = (/159,238 /)
    IDX_ADD  = 158
  CASE(239)
    nherz =  37
    nhart =  34
    nhug  =  81
    nchap =  81
    IDX_SRC  = (/  3,  5 /)
    IDX_SRB  = (/  6,  6 /)
    IDX_HRZ  = (/  7, 43 /)
    IDX_HAR  = (/ 44, 77 /)
    IDX_HUG  = (/ 78,158 /)
    IDX_CHA  = (/160,240 /)
    IDX_ADD  = 159
  CASE(241)
    nherz =  37
    nhart =  34
    nhug  =  82
    nchap =  82
    IDX_SRC  = (/  3,  5 /)
    IDX_SRB  = (/  6,  6 /)
    IDX_HRZ  = (/  7, 43 /)
    IDX_HAR  = (/ 44, 77 /)
    IDX_HUG  = (/ 78,159 /)
    IDX_CHA  = (/161,242 /)
    IDX_ADD  = 160
  CASE(243)
    nherz =  37
    nhart =  34
    nhug  =  83
    nchap =  83
    IDX_SRC  = (/  3,  5 /)
    IDX_SRB  = (/  6,  6 /)
    IDX_HRZ  = (/  7, 43 /)
    IDX_HAR  = (/ 44, 77 /)
    IDX_HUG  = (/ 78,160 /)
    IDX_CHA  = (/162,244 /)
    IDX_ADD  = 161
  CASE(245)
    nherz =  37
    nhart =  34
    nhug  =  84
    nchap =  84
    IDX_SRC  = (/  3,  5 /)
    IDX_SRB  = (/  6,  6 /)
    IDX_HRZ  = (/  7, 43 /)
    IDX_HAR  = (/ 44, 77 /)
    IDX_HUG  = (/ 78,161 /)
    IDX_CHA  = (/163,246 /)
    IDX_ADD  = 162
  CASE(247)
    nherz =  37
    nhart =  34
    nhug  =  85
    nchap =  85
    IDX_SRC  = (/  3,  5 /)
    IDX_SRB  = (/  6,  6 /)
    IDX_HRZ  = (/  7, 43 /)
    IDX_HAR  = (/ 44, 77 /)
    IDX_HUG  = (/ 78,162 /)
    IDX_CHA  = (/164,248 /)
    IDX_ADD  = 163
  CASE DEFAULT
    WRITE(*,*) substr, ': ERROR IN CHOICE OF NBANDS '// &
      ' ! NBANDS = ',nbands,' NOT VALID; ONLY 49, 55, OR 106.'
  END SELECT
  !
  status = 0
  RETURN
END SUBROUTINE fubrad_ini_param_dyn
