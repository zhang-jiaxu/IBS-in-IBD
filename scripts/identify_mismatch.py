# This script standardises the calprotectin string and numeric columns and flags rows where they disagree
import numpy as np
import pandas as pd

STRING = "calprotectin_observation_test_value_string"
NUMERIC = "calprotectin_observation_test_value_numeric"

fcp = pd.read_csv("data/combined_fcp_data.csv")

def filter_strings(df):
    # filter out rows where the string column is not a valid number or comparator
    string_raw = df[STRING].astype("string").str.strip()
    parsed = string_raw.str.extract(r"^(?P<comparator>[<>]=?)?\s*(?P<number>\d+(?:[.,]\d+)?)$")
    valid = parsed["number"].notna()
    return df[valid]

def standardize_decimals(df):
    # replace commas with dots in the string column for decimal standardization
    df[STRING] = df[STRING].str.replace(",", ".", regex=False)
    return df

## Parse the string column into a comparator ("<", ">", "<=", ">=" or none) and a number
## Accepts e.g. "100", "100.0", "<15", "< 15", "1,5"; anything else counts as unparseable
string_raw = fcp[STRING].astype("string").str.strip()
parsed = string_raw.str.extract(r"^(?P<comparator>[<>]=?)?\s*(?P<number>\d+(?:[.,]\d+)?)$")
fcp["string_comparator"] = parsed["comparator"].fillna("")
fcp["string_value"] = pd.to_numeric(parsed["number"].str.replace(",", ".", regex=False), errors="coerce")

## Numeric column as float, so "100" and "100.0" compare equal
numeric_raw = fcp[NUMERIC].astype("string").str.strip()
fcp["numeric_value"] = pd.to_numeric(numeric_raw, errors="coerce")

## Compare the two columns
has_string = string_raw.notna() & (string_raw != "")
has_numeric = numeric_raw.notna() & (numeric_raw != "")
string_ok = fcp["string_value"].notna()
numeric_ok = fcp["numeric_value"].notna()
equal = string_ok & numeric_ok & np.isclose(fcp["string_value"], fcp["numeric_value"])
has_comparator = fcp["string_comparator"] != ""

fcp["value_check"] = np.select(
    [
        ~has_string & ~has_numeric,
        has_string & ~string_ok,
        has_numeric & ~numeric_ok,
        string_ok & ~numeric_ok,
        numeric_ok & ~string_ok,
        equal & ~has_comparator,
        equal & has_comparator,
    ],
    [
        "both_empty",
        "string_unparseable",
        "numeric_unparseable",
        "string_only",            # e.g. "<15" with empty numeric
        "numeric_only",
        "match",                  # e.g. "100" vs 100.0
        "match_with_comparator",  # e.g. "<15" vs 15.0
    ],
    default="mismatch",
)

## Standardised value; left empty for rows that need review
fcp["calprotectin_value"] = np.select(
    [fcp["value_check"].isin(["match", "match_with_comparator", "numeric_only"]),
     fcp["value_check"] == "string_only"],
    [fcp["numeric_value"], fcp["string_value"]],
    default=np.nan,
)
fcp["calprotectin_comparator"] = fcp["string_comparator"]

print(fcp["value_check"].value_counts(dropna=False))

## Rows to review by hand
review = fcp[fcp["value_check"].isin(["mismatch", "string_unparseable", "numeric_unparseable"])]
review_cols = ["pseudo_id", "calprotectin_observation_start_date", STRING, NUMERIC,
               "calprotectin_observation_test_value_unit", "value_check"]
review[review_cols].to_csv("data/fcp_value_mismatches.csv", index=False)

# drop duplicates based on "patient_id", "calprotectin_observation_start_date" and "calprotectin_value"
fcp = fcp.drop_duplicates(subset=["patient_id", "calprotectin_observation_start_date", "calprotectin_value"])
# print shape
print(f"Shape after dropping duplicates: {fcp.shape}")

fcp.drop(columns=["string_comparator", "string_value", "numeric_value"]).to_csv(
    "data/combined_fcp_data_standardised.csv", index=False
)