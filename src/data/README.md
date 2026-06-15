# Data

Input data used by the model. This is a mix of **public reference data** and **mock
placeholder data** included so the example (`IE`) analysis is runnable out of the box.

> [!IMPORTANT]
> The files in `epidemiological/` are **MOCK DATA** — tiny, illustrative samples in the
> expected format, **not** real surveillance data. They are **not** sufficient to
> calibrate the model. **Replace them with your own data in the same format** before
> running a real analysis.

## `epidemiological/` — MOCK DATA (replace with your own)

- `RSV_weekly_counts.csv` — weekly RSV hospital admission counts.
  Columns: `date_wk_floor` (week start, `YYYY-MM-DD`), `season_name`, `case_counts`.
- `RSV_monthly_prop_age.csv` — 4-weekly age-group proportions of RSV hospital burden.
  Columns: `date_28days_floor` (**`YYYY-MM-DD`**), `age_gp_modelling`, `proportion`, `season_name`.

> Dates in both epidemiological files must be ISO `YYYY-MM-DD`; the loader parses them
> with `as.Date()` and other formats (e.g. `DD/MM/YYYY`) will fail.

## `population/`

- `monthly_births.csv` — monthly live births (Eurostat). Used for demographic ageing.
  Public reference data; replace with your country's births when adapting the model.

## `contact_matrices/`

- `contact_all.rdata` — synthetic social contact matrices (Prem et al. 2021), all countries.
  Public reference data.
