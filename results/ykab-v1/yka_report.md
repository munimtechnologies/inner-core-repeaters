# inner-core-repeaters 0.2.0: YKA run report

Processed 406 of 7365 catalogue events in 26.8 min.

Coverage by year (tried/ok): 1998: 2/2, 2000: 7/5, 2001: 10/10, 2002: 12/12, 2003: 15/15, 2004: 4/4, 2005: 15/15, 2006: 9/9, 2007: 17/17, 2008: 5/5, 2009: 9/9, 2010: 12/12, 2011: 10/10, 2012: 8/8, 2013: 8/8, 2014: 8/8, 2015: 9/9, 2016: 15/15, 2017: 11/11, 2018: 8/8, 2019: 4/4, 2020: 15/15, 2021: 61/61, 2022: 49/49, 2023: 25/24, 2024: 20/20, 2025: 19/19, 2026: 23/22

Element-years with |median static| >= 0.1 s (SNR >= 4):

* 2014 R4 0.116 n=7
* 2015 B6 0.214 n=8
* 2015 B7 0.101 n=8
* 2015 R3 0.126 n=8
* 2015 R8 -0.156 n=8
* 2017 B7 0.164 n=10
* 2017 B9 0.115 n=10
* 2017 R3 0.119 n=10
* 2018 B1 0.161 n=7
* 2018 B7 0.181 n=7
* 2018 R3 0.210 n=7
* 2019 B7 0.194 n=4
* 2019 R3 0.207 n=4
* 2020 B4 0.277 n=13
* 2020 B7 0.306 n=13
* 2020 R8 0.151 n=13

## Predictions on YKA

```
{
 "n_rows": 2972,
 "n_ok": 479,
 "n_in_sample": 406,
 "n_out_of_sample": 73,
 "P1_in_sample": {
  "n_pairs": 406,
  "n_cross_turn": 131,
  "n_multiplets": 30,
  "rho_B": -0.014998877067842586,
  "p_B": 0.7631825760272349,
  "rho_M": 0.019628495490564456,
  "p_M": 0.6933434653439753,
  "bootstrap": {
   "n": 1000,
   "frac_B_better": 0.124,
   "median": -0.039231564851247736,
   "p05": -0.12003434637413117,
   "p95": 0.05211856103570714
  },
  "testable": true,
  "holds": false
 },
 "P2_out_of_sample": {
  "n_pairs": 73,
  "n_cross_turn": 13,
  "n_multiplets": 11,
  "rho_B": 0.1452239911144021,
  "p_B": 0.2202402651784992,
  "rho_M": 0.08154387263976307,
  "p_M": 0.4928206271241679,
  "bootstrap": {
   "n": 1000,
   "frac_B_better": 0.85,
   "median": 0.04842604246719792,
   "p05": -0.033080228315967894,
   "p95": 0.26225908225606
  },
  "testable": true,
  "holds": false
 },
 "P3_new_event_best_match": {
  "cases": [
   {
    "ev2": "us6000rk5n",
    "t2": 2025.826,
    "n_partners": 5,
    "best_ev1": "usp000abph",
    "best_t1": 2001.21,
    "B_pick_t1": 2001.21,
    "M_pick_t1": 2007.224,
    "discriminating": true,
    "matches_B": true,
    "matches_M": false
   },
   {
    "ev2": "us6000tmzw",
    "t2": 2026.643,
    "n_partners": 5,
    "best_ev1": "usp000dsnt",
    "best_t1": 2005.44,
    "B_pick_t1": 2001.21,
    "M_pick_t1": 2007.224,
    "discriminating": true,
    "matches_B": false,
    "matches_M": false
   }
  ],
  "n_discriminating": 2,
  "n_match_B": 1,
  "n_match_M": 0,
  "testable": true,
  "holds": true
 },
 "P4_calibration": {
  "n_near": 101,
  "n_far": 20,
  "median_D_near": 0.11658280708757829,
  "median_D_far": 0.15433740941628948,
  "testable": true,
  "holds": true
 },
 "P5_fit_in_sample": {
  "n_pairs": 406,
  "best_rho": 0.08199197854370269,
  "best_T": 2009.75,
  "best_ratio": 5.0,
  "rho_M": 0.019628495490564456
 },
 "P5_fit_all": {
  "n_pairs": 479,
  "best_rho": 0.08348022816013138,
  "best_T": 2009.75,
  "best_ratio": 5.0,
  "rho_M": 0.03580397620565858
 },
 "robust_excess_D_in": {
  "n_pairs": 406,
  "n_cross_turn": 131,
  "n_multiplets": 30,
  "rho_B": -0.025722184342874,
  "p_B": 0.6053131030789903,
  "rho_M": 0.010940721285548873,
  "p_M": 0.8260453653506185,
  "bootstrap": {
   "n": 1000,
   "frac_B_better": 0.13,
   "median": -0.04077916466685373,
   "p05": -0.10663686215245434,
   "p95": 0.05717087097503383
  },
  "testable": true,
  "holds": false
 },
 "robust_excess_D_out": {
  "n_pairs": 73,
  "n_cross_turn": 13,
  "n_multiplets": 11,
  "rho_B": 0.14096630877452795,
  "p_B": 0.23421165163917165,
  "rho_M": 0.0884857460199926,
  "p_M": 0.4566128855809481,
  "bootstrap": {
   "n": 1000,
   "frac_B_better": 0.79,
   "median": 0.038589186422117096,
   "p05": -0.0365251059429394,
   "p95": 0.24567986025132524
  },
  "testable": true,
  "holds": false
 },
 "robust_station_D_in": {
  "n_pairs": 378,
  "n_cross_turn": 122,
  "n_multiplets": 28,
  "rho_B": -0.007308343856428653,
  "p_B": 0.8873787912924765,
  "rho_M": 0.02177461435035814,
  "p_M": 0.6730287876407172,
  "bootstrap": {
   "n": 1000,
   "frac_B_better": 0.131,
   "median": -0.03048747772597367,
   "p05": -0.23802217026498107,
   "p95": 0.05928838827593771
  },
  "testable": true,
  "holds": false
 },
 "robust_station_D_out": {
  "n_pairs": 73,
  "n_cross_turn": 13,
  "n_multiplets": 11,
  "rho_B": 0.1468591879550784,
  "p_B": 0.21503321709459325,
  "rho_M": 0.12051092188078492,
  "p_M": 0.3098287724626936,
  "bootstrap": {
   "n": 1000,
   "frac_B_better": 0.677,
   "median": 0.028571575490065455,
   "p05": -0.041400732669352173,
   "p95": 0.29540243806938293
  },
  "testable": true,
  "holds": false
 },
 "sensitivity_windows": {
  "D_c13": {
   "P1": {
    "n_pairs": 406,
    "n_cross_turn": 131,
    "n_multiplets": 30,
    "rho_B": 0.020293916845640986,
    "p_B": 0.6834995665031922,
    "rho_M": 0.045351666041321216,
    "p_M": 0.36205181733330644,
    "bootstrap": {
     "n": 300,
     "frac_B_better": 0.25333333333333335,
     "median": -0.024678969359938102,
     "p05": -0.05345356848482217,
     "p95": 0.1004664777315107
    },
    "testable": true,
    "holds": false
   },
   "P2": {
    "n_pairs": 73,
    "n_cross_turn": 13,
    "n_multiplets": 11,
    "rho_B": 0.15161051462421327,
    "p_B": 0.20039771416853241,
    "rho_M": 0.07204121930149329,
    "p_M": 0.5447242840067283,
    "bootstrap": {
     "n": 300,
     "frac_B_better": 0.8966666666666666,
     "median": 0.06275420760294412,
     "p05": -0.06020776198786664,
     "p95": 0.26553943940299146
    },
    "testable": true,
    "holds": false
   },
   "P4": {
    "n_near": 101,
    "n_far": 20,
    "median_D_near": 0.057601711448791004,
    "median_D_far": 0.10802591626597724,
    "testable": true,
    "holds": true
   }
  },
  "D_c310": {
   "P1": {
    "n_pairs": 406,
    "n_cross_turn": 131,
    "n_multiplets": 30,
    "rho_B": -0.0396420603317155,
    "p_B": 0.4256741895499337,
    "rho_M": -0.024014430910982636,
    "p_M": 0.6294827065335148,
    "bootstrap": {
     "n": 300,
     "frac_B_better": 0.18666666666666668,
     "median": -0.017631009721071008,
     "p05": -0.12213974992997313,
     "p95": 0.045955111369134105
    },
    "testable": true,
    "holds": false
   },
   "P2": {
    "n_pairs": 73,
    "n_cross_turn": 13,
    "n_multiplets": 11,
    "rho_B": 0.15389361964704434,
    "p_B": 0.19362408020520802,
    "rho_M": 0.12143650499814886,
    "p_M": 0.3061000329690909,
    "bootstrap": {
     "n": 300,
     "frac_B_better": 0.68,
     "median": 0.02792488851717057,
     "p05": -0.050068405446212536,
     "p95": 0.2688967631777041
    },
    "testable": true,
    "holds": false
   },
   "P4": {
    "n_near": 101,
    "n_far": 20,
    "median_D_near": 0.6210782372840865,
    "median_D_far": 0.5864268663969314,
    "testable": true,
    "holds": false
   }
  },
  "D_full": {
   "P1": {
    "n_pairs": 406,
    "n_cross_turn": 131,
    "n_multiplets": 30,
    "rho_B": -0.006301778715571819,
    "p_B": 0.8992672196686563,
    "rho_M": 0.030443251132906305,
    "p_M": 0.5407594037105667,
    "bootstrap": {
     "n": 300,
     "frac_B_better": 0.12333333333333334,
     "median": -0.03987521491104207,
     "p05": -0.10833196013637603,
     "p95": 0.06155733988639404
    },
    "testable": true,
    "holds": false
   },
   "P2": {
    "n_pairs": 73,
    "n_cross_turn": 13,
    "n_multiplets": 11,
    "rho_B": 0.15201160064173766,
    "p_B": 0.1991956558352233,
    "rho_M": 0.08009379242255955,
    "p_M": 0.5005696038973105,
    "bootstrap": {
     "n": 300,
     "frac_B_better": 0.89,
     "median": 0.054614820528278346,
     "p05": -0.06933518948921857,
     "p95": 0.254057859703021
    },
    "testable": true,
    "holds": false
   },
   "P4": {
    "n_near": 101,
    "n_far": 20,
    "median_D_near": 0.10408934270365133,
    "median_D_far": 0.13269807688838536,
    "testable": true,
    "holds": true
   }
  }
 },
 "cross_array_D_spearman": {
  "n": 29,
  "rho": 0.34926108374384224,
  "p": 0.06330177715735706
 }
}
```
