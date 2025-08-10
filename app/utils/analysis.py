import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import os
import uuid

UPLOAD_FOLDER = "uploads/"
STATIC_FOLDER = "static/plots/"

os.makedirs(STATIC_FOLDER, exist_ok=True)

def save_plot():
    """Generate unique filename for plot"""
    filename = f"{uuid.uuid4().hex}.png"
    full_path = os.path.join(STATIC_FOLDER, filename)
    return full_path, f"/{full_path}"

def perform_analysis(df, analysis_type):
    """Dispatches visual analysis based on selected type"""
    if analysis_type == "growth_curve":
        return growth_curve_plot(df)
    elif analysis_type == "amr_heatmap":
        return amr_heatmap_plot(df)
    return None, None

def growth_curve_plot(df):
    try:
        if 'Time' not in df.columns or 'OD' not in df.columns:
            return "Missing columns: 'Time' and 'OD'", None

        plt.figure(figsize=(8, 5))
        sns.lineplot(data=df, x="Time", y="OD", marker='o')
        plt.title("Bacterial Growth Curve")
        plt.xlabel("Time (hours)")
        plt.ylabel("Optical Density (OD)")
        plot_path, plot_url = save_plot()
        plt.savefig(plot_path)
        plt.close()
        return "Growth curve plotted successfully.", plot_url
    except Exception as e:
        return f"Error: {str(e)}", None

def amr_heatmap_plot(df):
    try:
        if df.select_dtypes(include=['number']).shape[1] < 2:
            return "Heatmap requires at least 2 numeric columns", None

        plt.figure(figsize=(10, 6))
        sns.heatmap(df.corr(), annot=True, cmap="viridis")
        plt.title("Antimicrobial Resistance Correlation")
        plot_path, plot_url = save_plot()
        plt.savefig(plot_path)
        plt.close()
        return "AMR heatmap generated successfully.", plot_url
    except Exception as e:
        return f"Error: {str(e)}", None

def growth_stats(df):
    """Calculate growth rate and doubling time"""
    try:
        if 'Time' not in df.columns or 'OD' not in df.columns:
            return {"Error": "Missing 'Time' and 'OD' columns"}
        df = df.sort_values("Time")
        df["logOD"] = df["OD"].apply(lambda x: stats.nan if x <= 0 else np.log(x))
        slope, _, _, _, _ = stats.linregress(df["Time"], df["logOD"])
        doubling_time = np.log(2) / slope
        return {
            "Growth Rate (1/hr)": round(slope, 4),
            "Doubling Time (hr)": round(doubling_time, 2)
        }
    except Exception as e:
        return {"Error": str(e)}

def inhibition_stats(df):
    """Zone of inhibition stats + t-test"""
    try:
        antibiotics = df.columns[1:]  # assume 1st column = sample ID
        means = df[antibiotics].mean()
        stds = df[antibiotics].std()

        stats_dict = {
            f"{ab} Mean (mm)": round(means[ab], 2)
            for ab in antibiotics
        }

        # Run t-test between first two antibiotics (if available)
        if len(antibiotics) >= 2:
            t_stat, p_value = stats.ttest_ind(df[antibiotics[0]], df[antibiotics[1]], nan_policy='omit')
            ttest_dict = {
                "Antibiotic 1": antibiotics[0],
                "Antibiotic 2": antibiotics[1],
                "t-statistic": t_stat,
                "p-value": p_value
            }
        else:
            ttest_dict = {"Note": "Need at least 2 antibiotic columns"}

        return stats_dict, ttest_dict
    except Exception as e:
        return {"Error": str(e)}, {}

def cfu_stats(df):
    """Basic summary for CFU counts"""
    try:
        return df.describe()
    except Exception as e:
        return pd.DataFrame({"Error": [str(e)]})

def compute_diversity(df, method="shannon"):
    from scipy.stats import entropy
    counts = df["count"].values
    proportions = counts / counts.sum()

    if method == "shannon":
        index = entropy(proportions)
        return pd.DataFrame({"Index": ["Shannon"], "Value": [round(index, 4)]})
    elif method == "simpson":
        simpson = 1 - sum(proportions**2)
        return pd.DataFrame({"Index": ["Simpson"], "Value": [round(simpson, 4)]})
    else:
        raise ValueError("Invalid method selected.")



def perform_growth_rate(initial_pop, final_pop, time):
    """
    Calculates microbial growth rate using the formula:
    μ = (ln(Nt) - ln(N0)) / t
    """
    try:
        growth_rate = (math.log(final_pop) - math.log(initial_pop)) / time
        return f"Growth Rate (μ): {growth_rate:.4f} per hour"
    except Exception as e:
        return f"Error in growth rate calculation: {str(e)}"


def perform_cfu(colony_count, dilution_factor, volume_plated):
    """
    CFU/mL = (Number of Colonies × Dilution Factor) / Volume Plated (mL)
    """
    try:
        cfu = (colony_count * dilution_factor) / volume_plated
        return f"Colony Forming Units (CFU/mL): {cfu:.2e}"
    except Exception as e:
        return f"Error in CFU calculation: {str(e)}"


def perform_zone_of_inhibition(diameter):
    """
    Measures zone of inhibition in mm (assumes a circle).
    """
    try:
        radius = diameter / 2
        area = math.pi * (radius ** 2)
        return f"Zone of Inhibition Area: {area:.2f} mm²"
    except Exception as e:
        return f"Error in zone of inhibition calculation: {str(e)}"


def perform_diversity_index(counts, method="shannon"):
    """
    Calculate biodiversity index using either Shannon or Simpson index.
    """
    try:
        total = sum(counts)
        if method == "shannon":
            proportions = [count / total for count in counts if count > 0]
            shannon_index = -sum(p * math.log(p) for p in proportions)
            return f"Shannon Diversity Index (H'): {shannon_index:.4f}"

        elif method == "simpson":
            proportions = [count / total for count in counts]
            simpson_index = 1 - sum(p ** 2 for p in proportions)
            return f"Simpson Diversity Index (1 - D): {simpson_index:.4f}"

        else:
            return "Unknown diversity index method selected."

    except Exception as e:
        return f"Error in diversity index calculation: {str(e)}"


def calculate_reed_muench(dilutions, infected, total):
    prop_infected = [i / t if t != 0 else 0 for i, t in zip(infected, total)]
    cumulative_above = [sum(prop_infected[i:]) for i in range(len(prop_infected))]
    cumulative_below = [sum(prop_infected[:i+1]) for i in range(len(prop_infected))]

    diffs = [abs(a - b) for a, b in zip(cumulative_above, cumulative_below)]
    min_index = diffs.index(min(diffs))

    # Interpolation to estimate 50% endpoint
    if min_index + 1 < len(dilutions):
        d1, d2 = dilutions[min_index], dilutions[min_index + 1]
        r1, r2 = cumulative_above[min_index], cumulative_above[min_index + 1]
        tcid50 = d1 + ((0.5 - r1) * (d2 - d1) / (r2 - r1))
        return round(tcid50, 4)
    else:
        return "Calculation not possible with provided data."

