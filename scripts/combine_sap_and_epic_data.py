# This script just combines SAP and EPIC data into a single dataset
import pandas as pd
# rename columns for concat
def harmonise_columns(df):
    df = df.rename(columns=lambda c: c.replace("calprotectine_", "calprotectin_"))
    for col in [c for c in df.columns if "_measurement_" in c]:
        target = col.replace("_measurement_", "_observation_")
        if target in df.columns:
            df[target] = df[target].combine_first(df[col])
            df = df.drop(columns=col)
        else:
            df = df.rename(columns={col: target})
    return df

def parse_dates(col):
    s = col.astype("string").str.strip()
    is_js = s.str.contains(r"GMT[+-]\d{4}", na=False)
    js = pd.to_datetime(
        s[is_js].str.replace(r"\s*\(.*\)$", "", regex=True),
        format="%a %b %d %Y %H:%M:%S GMT%z",
        utc=True,
    ).dt.tz_convert("Europe/Amsterdam").dt.tz_localize(None)
    other = pd.to_datetime(s[~is_js], errors="coerce")
    return pd.concat([js, other]).reindex(s.index).dt.normalize()

# 1. combine crp data
sap_data_crp = pd.read_csv("data/raw/ibs-in-ibd-sap-5-crp.csv")
epic_data_crp = pd.read_csv("data/raw/ibs-in-ibd-epic-1-crp.csv")
# harmonise columns for CRP data before combining
sap_data_crp = harmonise_columns(sap_data_crp)
epic_data_crp = harmonise_columns(epic_data_crp)
## Combine and sort CRP data by pseudo_id and observation start date and within each pseudo_id by the earliest observation date
combined_crp_data = pd.concat([sap_data_crp, epic_data_crp], ignore_index=True)
combined_crp_data["crp_observation_start_date"] = pd.to_datetime(combined_crp_data["crp_observation_start_date"])
# sort by "pseudo_id" and then "crp_observation_start_date" 
combined_crp_data["_group_min"] = (
    combined_crp_data.groupby("pseudo_id")["crp_observation_start_date"].transform("min")
)

combined_crp_data = (
    combined_crp_data
    .sort_values(["_group_min", "pseudo_id", "crp_observation_start_date"])
    .drop(columns="_group_min")
    .reset_index(drop=True)
)
# save combined CRP data
combined_crp_data.to_csv("data/combined_crp_data.csv", index=False)
print(f"Combined CRP data shape: {combined_crp_data.shape}")
# 2. combine fcp data with extra QoC FCP data
sap_data_fcp = harmonise_columns(pd.read_csv("data/raw/ibs-in-ibd-sap-6-calprotectine.csv"))
epic_data_fcp = harmonise_columns(pd.read_csv("data/raw/ibs-in-ibd-epic-2-calprotectin.csv"))
combined_fcp_data = pd.concat(
    [sap_data_fcp.assign(source="sap"), epic_data_fcp.assign(source="epic")], ignore_index=True
)

# load the additional dataset of FCP measurements from QoC test results
extra_qoc_data_fcp = pd.read_excel("data/raw/QoC_test_result_mumc_2024_01_18.xlsx")
print(f"QoC rows loaded: {len(extra_qoc_data_fcp)}")
# filter "test_valid" = true and drop duplicates by "patient_code", "test_date" and "result"
extra_qoc_data_fcp = extra_qoc_data_fcp[extra_qoc_data_fcp["test_valid"] == True]
extra_qoc_data_fcp = extra_qoc_data_fcp.drop_duplicates(subset=["patient_code", "test_date", "result"])
print(f"QoC rows after validity filter and deduplication: {len(extra_qoc_data_fcp)}")

# map QoC columns to the SAP/EPIC names; "result" goes into the string column so the standardisation script parses it
extra_qoc_data_fcp = extra_qoc_data_fcp.rename(columns={
    "patient_code": "patient_id",
    "test_date": "calprotectin_observation_start_date",
    "result": "calprotectin_observation_test_value_string",
})[["patient_id", "calprotectin_observation_start_date", "calprotectin_observation_test_value_string"]]
extra_qoc_data_fcp["source"] = "qoc"

# keep only QoC patients that are in the SAP/EPIC data
extra_qoc_data_fcp = extra_qoc_data_fcp[extra_qoc_data_fcp["patient_id"].isin(combined_fcp_data["patient_id"])]
print(f"QoC rows matched to a SAP/EPIC patient_id: {len(extra_qoc_data_fcp)}")

# fill pseudo_id from the SAP/EPIC data (same patient_id -> same pseudo_id)
id_map = combined_fcp_data[["patient_id", "pseudo_id"]].dropna().drop_duplicates()
# NOTE: one patient_id might have different pseudo_ids from EPIC/SAP; we take the first one encountered
extra_qoc_data_fcp["pseudo_id"] = extra_qoc_data_fcp["patient_id"].map(
    id_map.drop_duplicates("patient_id").set_index("patient_id")["pseudo_id"]
)

# label QoC measurements
extra_qoc_data_fcp["calprotectin_observation_description"] = "Calprotectin (QoC home test)"
# fill in unit as "ug/g" in the QoC data
extra_qoc_data_fcp["calprotectin_observation_test_value_unit"] = "ug/g"

combined_fcp_data_final = pd.concat([combined_fcp_data, extra_qoc_data_fcp], ignore_index=True)
combined_fcp_data_final["calprotectin_observation_start_date"] = parse_dates(
    combined_fcp_data_final["calprotectin_observation_start_date"]
)

# sort groups by each patient's earliest date, then by date within each patient
combined_fcp_data_final["_group_min"] = (
    combined_fcp_data_final.groupby("patient_id")["calprotectin_observation_start_date"].transform("min")
)
combined_fcp_data_final = (
    combined_fcp_data_final
    .sort_values(["_group_min", "patient_id", "calprotectin_observation_start_date"])
    .drop(columns="_group_min")
    .reset_index(drop=True)
)

print(f"Final combined FCP data shape: {combined_fcp_data_final.shape}")
print(combined_fcp_data_final["source"].value_counts())
combined_fcp_data_final.to_csv("data/combined_fcp_data.csv", index=False)