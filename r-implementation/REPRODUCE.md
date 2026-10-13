# Reproducing the tables and figures

Scripts read `MIO_DIR` (the working directory holding `scs.R`, `src_scs.cpp`, `sim/`, `data/`); set it
before running. Heavy simulations were run as SLURM arrays (`run_*.sh`); each script documents its
environment variables and seeds.

| Paper object | Script |
|---|---|
| Fig. 2 (phase map), Fig. 3 (selection), Fig. 4 (scaling), Supp. S3 figures and Table S2 | `sim/run_sims.R`, `sim/figures.R`, `sim/tables.R` |
| Table 1 (identical-conditions benchmark) | `sim/bench_table.R`, `sim/bench_table_aggregate.R` |
| Supp. S2 brute force and L0Learn head-to-head | `sim/bench_l0learn.R` |
| Supp. S2 Table S1 (greedy / marginal / SCS vs brute force) | `sim/greedy_vs_brute.R` |
| Supp. S3 harder regimes (Table S3) | `sim/run_hard.R`, `sim/dgp_hard.R` |
| Supp. S3.8 ridge and restart sensitivity | `sim/restart_robustcore_sim.R`, `scs_eps_sensitivity()` in `scs.R` |
| Supp. S6 consistency check, subsampling scheme | `sim/consistency_sim.R`, `sim/subsample_scheme_sim.R`, `sim/aggregate_strengthen.R` |
| Supp. S4 null calibration of the CPSS reference (Table S6) | `sim/null_pfer_calibration.R` |
| Supp. S4 planted-signal calibration | `sim/planted_pfer_calibration.R`, `sim/aggregate_planted.R` |
| Table 2, Supp. S3.9 (per-cluster z-test comparators) | `sim/bh_comparator_sim.R`, `sim/aggregate_bh.R`, `sim/table_bh.R` |
| Sec. 6 NY cohort, Figs. 5-6 | `pull_sparcs_hf.R`, `figures_realdata.R` |
| Sec. 6 TX cohort | `tx_hf_analysis.R` |
| Supp. S4 risk-adjustment sensitivity (Tables S4-S5) | `realdata_sensitivity.R`, `tx_sensitivity.R` |
| Supp. S4 sample-split audit check | `audit_curve_splitval.R` |
| Flagged-hospital table | `flagged_table.R` |
