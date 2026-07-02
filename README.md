# RSVinea — ODE Model of RSV Transmission and Disease

> **Work in progress — adult vaccination branch**
> This branch is under active development and extends the model to incorporate adult vaccination. The code, parameters, and results here are **not yet validated** and should not be used for policy or clinical decisions. Please refer to the `main` branch for the current stable release.

RSVinea is an age-structured, compartmental ordinary differential equation (ODE) model for simulating Respiratory Syncytial Virus (RSV) transmission and disease burden. It was developed as part of the [RespiCompass](https://github.com/european-modelling-hubs/RespiCompass) European respiratory virus modelling hub. The model is country-agnostic — each analysis is driven by a per-country configuration file — and ships with a worked example (`IE`).

---

## Model overview

The model tracks RSV transmission through **31 fine-grained age groups** (0–1 month to 65+ years) and **four infection-history tiers** (naive, one, two, three or more prior infections). Each age group × tier combination has the following compartments:

| Compartment | Description |
|-------------|-------------|
| `V0` | Vaccinated (waning maternal immunity / mAbs), naive tier only |
| `S0`–`S3` | Susceptible, by number of prior infections |
| `E0`–`E3` | Exposed (latent, pre-infectious) |
| `I0`–`I3` | Infectious |
| `H0`–`H3` | Hospitalised |
| `R0`–`R3` | Recovered (temporary immunity, wanes back to susceptible) |
| `D0`–`D3` | Deceased (cumulative absorbing state) |

### Key features

- **Seasonality** — cosine wave on the transmission rate with configurable amplitude, peak day, and shape exponent
- **Age-structured mixing** — Prem et al. 2021 synthetic contact matrices rebanded to the fine model age groups
- **Vaccination** — age-dependent waning protection; configurable season windows and a one-shot catch-up campaign
- **Year-to-year season scalars** — `season2_effect` and `season3_effect` allow between-season transmission differences
- **Adaptive calibration** — likelihood-based iterative parameter fitting against RespiCompass hospital admission data using Latin Hypercube Sampling
- **Scenario analysis** — parallel simulation across calibrated parameter sets, uncertainty samples, and intervention scenarios
- **Demographic ageing** — monthly cohort shifting and birth replenishment implemented as discrete ODE events

---

## Repository structure

```
src/
├── launch.R               # Main entry point — calibrate, simulate, plot
├── results_evaluation.R   # Post-hoc evaluation and RespiCompass submission formatting
│
├── R/                     # Internal model scripts (sourced by dependencies.R)
│   ├── model.R            # ODE system, ageing events, initial conditions, output formatting
│   ├── calibration.R      # Adaptive parameter fitting pipeline
│   ├── scenarios.R        # Parallel scenario simulation
│   ├── uncertainty.R      # Latin Hypercube Sampling for parameter uncertainty
│   ├── results.R          # Post-processing and summary table generation
│   ├── plotting.R         # All visualisation functions
│   ├── load_data.R        # Data loading from RespiCompass
│   ├── parse_input.R      # YAML parsing, scenario application, parameter overrides
│   ├── options.R          # Global options, data URLs, parallelisation settings
│   ├── directories.R      # Output directory management
│   ├── auxiliary.R        # General utility functions
│   ├── unit_tests.R       # Single-simulation regression testing utilities
│   └── dependencies.R     # Package management and source loading
│
├── config/
│   └── default.yaml       # Default model parameters and calibration settings
├── input/
│   └── IE.yaml            # Example country parameter overrides
└── data/
    ├── contact_matrices/  # Prem et al. 2021 contact matrices (pre-computed .rdata)
    ├── epidemiological/   # RespiCompass hospital admissions and age-burden CSVs
    └── population/        # Monthly birth data (Eurostat)
```

---

## Quick start

### Requirements

- **R** ≥ 4.3.0
- All required packages are installed automatically via [pacman](https://github.com/trinker/pacman) on first run.

### Running the example analysis

1. Open `src/RSVinea.Rproj` in RStudio — this sets the working directory correctly.
2. Open and **source** `src/launch.R`.

The script runs three sequential steps controlled by `do_step`:

| Step | Function | Description |
|------|----------|-------------|
| 1 | `run_calibration(o)` | Adaptive fitting to RespiCompass hospital admissions |
| 2 | `run_scenarios(o)` | Parallel simulation across all defined intervention scenarios |
| 3 | `run_results(o)` | Scenario comparison plots and parameter summary tables |

Outputs are written to `src/output/`:

```
src/output/
├── 0_testing/          # Unit test outputs
├── 1_calibration/IE/   # Per-round parameter samples, likelihoods, and fit figures
├── 2_scenarios/IE/     # Raw and summarised scenario simulation results
└── 3_results/IE/       # Figures (.png) and CSV summary tables
```

---

## Configuration

### Parameter files

All model parameters are initialised in `src/config/default.yaml`. Country-specific overrides are placed in `src/input/<analysis_name>.yaml` (e.g. `IE.yaml`). The launch script selects the input file via `analysis_name = "IE"` in `set_options()`.

Parameters that should vary stochastically across simulations are specified with an `uncertainty:` block in the YAML (e.g. using a normal or uniform distribution). Calibration prior bounds are set under `calibration_parameters:`.

### Key options (`options.R`)

| Option | Default | Description |
|--------|---------|-------------|
| `n_best_samples` | 10 | Number of best-fitting calibration samples to carry forward |
| `n_parameter_sets` | 10 | Number of uncertainty parameter sets per scenario |
| `quantiles` | `[0.05, 0.95]` | Credible interval bounds for summary plots |
| `overwrite_samples` | `TRUE` | Re-run calibration even if samples already exist |

### Adapting for a new country

1. Copy `src/input/IE.yaml` and rename it `<ISO2>.yaml`.
2. Update `contact_matrix_countries` (ISO-3 code), `calibration_options.country` (ISO-2 code), `data_end`, and `data_days`.
3. In the launch script, change `analysis_name = "IE"` to `analysis_name = "<ISO2>"`.
4. Provide country-specific birth data or update `o$births_url` in `options.R`.

---

## Data sources

| Data | Source | Usage |
|------|--------|-------|
| Weekly RSV hospital admissions | [RespiCompass](https://github.com/european-modelling-hubs/RespiCompass) | Primary calibration target |
| 4-weekly age-stratified hospital burden | RespiCompass | Age-resolved calibration target |
| Population estimates by age group | RespiCompass population data (GitHub URL) | Age-group sizes |
| Monthly birth counts | Eurostat (`data/population/monthly_births.csv`) | Demographic ageing events |
| Social contact matrices | [Prem et al. 2021](https://doi.org/10.1371/journal.pcbi.1009098) | Age-structured mixing |

> **Note:** The RSV hospital admission and burden files in `src/data/epidemiological/` are **mock placeholder data** (small samples in the expected format), not real surveillance data, and are not sufficient to calibrate the model. Replace them with your own data before running a real analysis — see [`src/data/README.md`](src/data/README.md).

---

## Key parameters

| Parameter | Description | Calibrated? |
|-----------|-------------|-------------|
| `beta_A` | Baseline transmission rate (strain A) | Yes |
| `amplitude` | Seasonal forcing amplitude | Yes |
| `peak_day` | Days after simulation start to the first seasonal peak (not a calendar day-of-year) | Yes |
| `seasonality_exponent` | Shape of the cosine seasonality curve | Yes |
| `season2_effect`, `season3_effect` | Year-to-year scalars on beta | Yes |
| `omega_1`–`omega_4` | Waning immunity durations (days) | Yes |
| `p_hosp_A` | Baseline hospitalisation probability (oldest age group) | Yes |
| `init_inf` | Global seeding intensity scalar | Yes |
| `vacc_coverage` | Routine vaccination coverage | Scenario-dependent |
| `vacc_IE` | Vaccine immunisation effectiveness | Uncertain (normal distribution) |
| `vaccination_catch_up_coverage` | Catch-up campaign coverage | Yes |

---

## Scenarios

Three scenarios are defined in `default.yaml` (overridable in the country YAML):

| Scenario ID | Description |
|-------------|-------------|
| `baseline` | Current vaccination programme as configured |
| `no_vacc` | No vaccination: `vacc_coverage = 0`, `vacc_IE = 0`, no catch-up |
| `high_vacc` | Increased coverage: 95% routine, 80% catch-up |

To add a new scenario, append a block under `scenarios:` in the country YAML:

```yaml
scenarios:
- id: my_scenario
  name: "My scenario description"
  vacc_coverage: 0.70
```

---

## Licence

GNU General Public License v2 (GPL-2). See [LICENSE](LICENSE).
