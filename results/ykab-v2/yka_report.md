# inner-core-repeaters 0.2.0: YKA run report

Processed 406 of 7366 catalogue events in 0.6 min.

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
 "n_rows": 2971,
 "n_ok": 1650,
 "n_in_sample": 1463,
 "n_out_of_sample": 187,
 "P1_in_sample": {
  "n_pairs": 1463,
  "n_cross_turn": 506,
  "n_multiplets": 44,
  "rho_B": -0.016110139375218827,
  "p_B": 0.5380840013856762,
  "rho_M": 0.017121189409024223,
  "p_M": 0.5128800790871428,
  "bootstrap": {
   "n": 1000,
   "frac_B_better": 0.131,
   "median": -0.03204925447073753,
   "p05": -0.059331309300628285,
   "p95": 0.016435793260838476
  },
  "testable": true,
  "holds": false
 },
 "P2_out_of_sample": {
  "n_pairs": 187,
  "n_cross_turn": 29,
  "n_multiplets": 22,
  "rho_B": 0.06011546691428802,
  "p_B": 0.41376039489284355,
  "rho_M": 0.02933101861197465,
  "p_M": 0.6902675517345951,
  "bootstrap": {
   "n": 1000,
   "frac_B_better": 0.848,
   "median": 0.025027231825958236,
   "p05": -0.020301022058498867,
   "p95": 0.07851635977499519
  },
  "testable": true,
  "holds": false
 },
 "P3_new_event_best_match": {
  "cases": [
   {
    "ev2": "us6000rii9",
    "t2": 2025.805,
    "n_partners": 2,
    "best_ev1": "usp000f0ed",
    "best_t1": 2006.961,
    "B_pick_t1": 2002.491,
    "M_pick_t1": 2006.961,
    "discriminating": true,
    "matches_B": false,
    "matches_M": true
   },
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
    "ev2": "us6000runq",
    "t2": 2025.95,
    "n_partners": 2,
    "best_ev1": "usp000bt5x",
    "best_t1": 2003.218,
    "B_pick_t1": 2003.218,
    "M_pick_t1": 2008.131,
    "discriminating": true,
    "matches_B": true,
    "matches_M": false
   },
   {
    "ev2": "us6000t7ng",
    "t2": 2026.476,
    "n_partners": 4,
    "best_ev1": "usp0008jpz",
    "best_t1": 1998.265,
    "B_pick_t1": 2001.078,
    "M_pick_t1": 2006.255,
    "discriminating": true,
    "matches_B": false,
    "matches_M": false
   },
   {
    "ev2": "us6000thwz",
    "t2": 2026.589,
    "n_partners": 4,
    "best_ev1": "usp000by6r",
    "best_t1": 2003.399,
    "B_pick_t1": 2001.078,
    "M_pick_t1": 2006.255,
    "discriminating": true,
    "matches_B": false,
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
  "n_discriminating": 6,
  "n_match_B": 2,
  "n_match_M": 1,
  "testable": true,
  "holds": true
 },
 "P4_calibration": {
  "n_near": 321,
  "n_far": 92,
  "median_D_near": 0.09827435544036078,
  "median_D_far": 0.08383037697471835,
  "testable": true,
  "holds": false
 },
 "P5_fit_in_sample": {
  "n_pairs": 1463,
  "best_rho": 0.08088764712918278,
  "best_T": 2010.25,
  "best_ratio": 5.0,
  "rho_M": 0.017121189409024223
 },
 "P5_fit_all": {
  "n_pairs": 1650,
  "best_rho": 0.07879111200274588,
  "best_T": 2010.25,
  "best_ratio": 5.0,
  "rho_M": 0.01276285567714876
 },
 "robust_excess_D_in": {
  "n_pairs": 1463,
  "n_cross_turn": 506,
  "n_multiplets": 44,
  "rho_B": -0.022865270345880668,
  "p_B": 0.38214841262958066,
  "rho_M": 0.015429555507562753,
  "p_M": 0.555392766677755,
  "bootstrap": {
   "n": 1000,
   "frac_B_better": 0.148,
   "median": -0.0379040355374477,
   "p05": -0.07543273154264714,
   "p95": 0.02365235244500904
  },
  "testable": true,
  "holds": false
 },
 "robust_excess_D_out": {
  "n_pairs": 187,
  "n_cross_turn": 29,
  "n_multiplets": 22,
  "rho_B": 0.06111011198015115,
  "p_B": 0.40606300301797016,
  "rho_M": 0.036350423366451466,
  "p_M": 0.6213663324564905,
  "bootstrap": {
   "n": 1000,
   "frac_B_better": 0.822,
   "median": 0.02087465265620423,
   "p05": -0.019779696096966677,
   "p95": 0.0744877963322671
  },
  "testable": true,
  "holds": false
 },
 "robust_station_D_in": {
  "n_pairs": 593,
  "n_cross_turn": 201,
  "n_multiplets": 42,
  "rho_B": -0.025702883628684642,
  "p_B": 0.5321754037772946,
  "rho_M": 0.012134988737012348,
  "p_M": 0.7680751348278118,
  "bootstrap": {
   "n": 1000,
   "frac_B_better": 0.025,
   "median": -0.042315043141301294,
   "p05": -0.16452097027852816,
   "p95": -0.007010177500533931
  },
  "testable": true,
  "holds": false
 },
 "robust_station_D_out": {
  "n_pairs": 118,
  "n_cross_turn": 20,
  "n_multiplets": 21,
  "rho_B": 0.039686800404646866,
  "p_B": 0.669606643269999,
  "rho_M": 0.04923325262308313,
  "p_M": 0.5965015795476748,
  "bootstrap": {
   "n": 1000,
   "frac_B_better": 0.379,
   "median": -0.011575332683372592,
   "p05": -0.08099248434321257,
   "p95": 0.0691616591086879
  },
  "testable": true,
  "holds": false
 },
 "sensitivity_windows": {
  "D_c13": {
   "P1": {
    "n_pairs": 1463,
    "n_cross_turn": 506,
    "n_multiplets": 44,
    "rho_B": -0.0012196291614252994,
    "p_B": 0.9628240916321398,
    "rho_M": 0.05142148430621796,
    "p_M": 0.049246697589438754,
    "bootstrap": {
     "n": 300,
     "frac_B_better": 0.013333333333333334,
     "median": -0.051772880317850434,
     "p05": -0.07327086947390116,
     "p95": -0.010666791547293495
    },
    "testable": true,
    "holds": false
   },
   "P2": {
    "n_pairs": 187,
    "n_cross_turn": 29,
    "n_multiplets": 22,
    "rho_B": 0.03218649411471084,
    "p_B": 0.6618895912427597,
    "rho_M": 0.00856275623121277,
    "p_M": 0.9074061236514479,
    "bootstrap": {
     "n": 300,
     "frac_B_better": 0.7433333333333333,
     "median": 0.02263090455158155,
     "p05": -0.05123920824008559,
     "p95": 0.06353513905413669
    },
    "testable": true,
    "holds": false
   },
   "P4": {
    "n_near": 321,
    "n_far": 92,
    "median_D_near": 0.04377619335177785,
    "median_D_far": 0.046575415574863166,
    "testable": true,
    "holds": true
   }
  },
  "D_c310": {
   "P1": {
    "n_pairs": 1463,
    "n_cross_turn": 506,
    "n_multiplets": 44,
    "rho_B": -0.005992848374010965,
    "p_B": 0.8188473986745369,
    "rho_M": -0.009571999217065663,
    "p_M": 0.7145019399592102,
    "bootstrap": {
     "n": 300,
     "frac_B_better": 0.69,
     "median": 0.008481265909077223,
     "p05": -0.022254910431985592,
     "p95": 0.054634956454639776
    },
    "testable": true,
    "holds": false
   },
   "P2": {
    "n_pairs": 187,
    "n_cross_turn": 29,
    "n_multiplets": 22,
    "rho_B": 0.08907028213419269,
    "p_B": 0.22541216823862967,
    "rho_M": 0.08331345266627271,
    "p_M": 0.2569509456471421,
    "bootstrap": {
     "n": 300,
     "frac_B_better": 0.5033333333333333,
     "median": 0.00018125385174661696,
     "p05": -0.04199781656895157,
     "p95": 0.09533211668356387
    },
    "testable": true,
    "holds": false
   },
   "P4": {
    "n_near": 321,
    "n_far": 92,
    "median_D_near": 0.70221471739656,
    "median_D_far": 0.5864268663969314,
    "testable": true,
    "holds": false
   }
  },
  "D_full": {
   "P1": {
    "n_pairs": 1463,
    "n_cross_turn": 506,
    "n_multiplets": 44,
    "rho_B": -0.01186507816858053,
    "p_B": 0.6502192689421328,
    "rho_M": 0.021050768362002988,
    "p_M": 0.4210646139144129,
    "bootstrap": {
     "n": 300,
     "frac_B_better": 0.12666666666666668,
     "median": -0.03235046408071887,
     "p05": -0.05933236026428273,
     "p95": 0.01455439301536699
    },
    "testable": true,
    "holds": false
   },
   "P2": {
    "n_pairs": 187,
    "n_cross_turn": 29,
    "n_multiplets": 22,
    "rho_B": 0.06095045493083364,
    "p_B": 0.4072928158902004,
    "rho_M": 0.027064622567065138,
    "p_M": 0.7131052005022936,
    "bootstrap": {
     "n": 300,
     "frac_B_better": 0.89,
     "median": 0.028418006671578493,
     "p05": -0.016467537568087214,
     "p95": 0.07470602904496777
    },
    "testable": true,
    "holds": false
   },
   "P4": {
    "n_near": 321,
    "n_far": 92,
    "median_D_near": 0.08679385486371394,
    "median_D_far": 0.07136985536027901,
    "testable": true,
    "holds": false
   }
  }
 },
 "cross_array_D_spearman": {
  "n": 132,
  "rho": 0.4456282969742189,
  "p": 8.580415865354038e-08
 }
}
```
