# This script identifies calprotectin remission periods (FCP < 150 sustained for >= 1 year) per patient
import pandas as pd

REMISSION_CUTOFF = 150
MIN_REMISSION_DAYS = 365
DATE = "calprotectin_observation_start_date"
CALPROTECTIN_NUMERIC_COL = "calprotectin_value" #standardized value

fcp = pd.read_csv("data/combined_fcp_data_standardised.csv", parse_dates=[DATE])

fcp = fcp.dropna(subset=["patient_id", DATE, CALPROTECTIN_NUMERIC_COL])

## One measurement per patient per day; if a day has several, keep the highest (conservative)
daily = (
    fcp.groupby(["patient_id", DATE])[CALPROTECTIN_NUMERIC_COL].max()
    .reset_index()
    .sort_values(["patient_id", DATE])
)
daily["in_remission"] = daily[CALPROTECTIN_NUMERIC_COL] < REMISSION_CUTOFF

## Consecutive measurements with the same status form one run; a new run starts on every status change or new patient
daily["run_id"] = (daily["in_remission"] != daily.groupby("patient_id")["in_remission"].shift()).cumsum()

## Remission periods: first to last consecutive measurement < cutoff
periods = (
    daily[daily["in_remission"]]
    .groupby(["patient_id", "run_id"])
    .agg(start_date=(DATE, "min"), end_date=(DATE, "max"), n_measurements=(CALPROTECTIN_NUMERIC_COL, "size"))
    .reset_index()
    .drop(columns="run_id")
)
periods["duration_days"] = (periods["end_date"] - periods["start_date"]).dt.days

## All remission periods of at least 1 year, numbered per patient
long_periods = periods[periods["duration_days"] >= MIN_REMISSION_DAYS].copy()
long_periods["period_number"] = long_periods.groupby("patient_id").cumcount() + 1

## Per-patient summary, including patients who never reached remission
summary = pd.DataFrame({"patient_id": daily["patient_id"].unique()})
summary["longest_remission_days"] = (
    summary["patient_id"].map(periods.groupby("patient_id")["duration_days"].max()).fillna(0).astype(int)
)
summary["n_remission_periods_1yr"] = (
    summary["patient_id"].map(long_periods.groupby("patient_id").size()).fillna(0).astype(int)
)
summary["remission_1yr"] = summary["n_remission_periods_1yr"] > 0

long_periods.to_csv("data/fcp_remission_periods_1yr.csv", index=False)
summary.to_csv("data/fcp_remission_summary.csv", index=False)