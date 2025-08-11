# analysis.py 

import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from scipy.optimize import curve_fit
from sklearn.decomposition import PCA
from sklearn.metrics import pairwise_distances
from scipy.interpolate import UnivariateSpline
from scipy.stats import entropy
import networkx as nx
import os
import uuid
from Bio import Phylo  #for phylogenetic diversity
from rdkit import Chem  # For chemistry-related if needed
# Note: Based on available libraries in the environment: numpy, scipy, pandas, matplotlib, biopython, rdkit, networkx, etc.

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
    # Existing
    if analysis_type == "growth_curve":
        return growth_curve_plot(df)
    elif analysis_type == "amr_heatmap":
        return amr_heatmap_plot(df)
    elif analysis_type == "time_kill_curve":
        return time_kill_curve_plot(df)
    elif analysis_type == "survivorship_curve":
        return survivorship_curve_plot(df)
    elif analysis_type == "ph_change_curve":
        return ph_change_curve_plot(df)
    elif analysis_type == "pca_microbiome":
        return pca_microbiome_plot(df)
    elif analysis_type == "rarefaction_curve":
        return rarefaction_curve_plot(df)
    elif analysis_type == "volcano_plot":
        return volcano_plot(df)
    elif analysis_type == "gene_expression_heatmap":
        return gene_expression_heatmap_plot(df)
    elif analysis_type == "cooccurrence_network":
        return cooccurrence_network_plot(df)
    # New visual analyses (expanding to more)
    elif analysis_type == "biofilm_biomass_plot":
        return biofilm_biomass_plot(df)  # Assumes columns: Treatment, Biomass
    elif analysis_type == "motility_plot":
        return motility_plot(df)  # Assumes columns: Time, Distance
    elif analysis_type == "fermentation_gas_plot":
        return fermentation_gas_plot(df)  # Assumes columns: Time, GasVolume
    elif analysis_type == "oxygen_consumption_plot":
        return oxygen_consumption_plot(df)  # Assumes columns: Time, O2Level
    elif analysis_type == "quorum_sensing_plot":
        return quorum_sensing_plot(df)  # Assumes columns: Time, Luminescence
    elif analysis_type == "enzyme_kinetics_plot":
        return enzyme_kinetics_plot(df)  # Michaelis-Menten plot
    elif analysis_type == "phylogenetic_tree":
        return phylogenetic_tree_plot(df)  # Basic tree plot if newick or distance matrix
    elif analysis_type == "dose_response_curve":
        return dose_response_curve_plot(df)  # Assumes Concentration, Response
    elif analysis_type == "synergy_heatmap":
        return synergy_heatmap_plot(df)  # For drug combinations
    elif analysis_type == "flow_cytometry_histogram":
        return flow_cytometry_histogram_plot(df)  # Assumes fluorescence data
    elif analysis_type == "qpcr_ct_plot":
        return qpcr_ct_plot(df)  # Cycle threshold bar plot
    elif analysis_type == "western_blot_densitometry":
        return western_blot_densitometry_plot(df)  # Bar plot of band intensities
    elif analysis_type == "elisa_standard_curve":
        return elisa_standard_curve_plot(df)  # Standard curve for quantification
    elif analysis_type == "hemagglutination_titer_plot":
        return hemagglutination_titer_plot(df)  # Titer levels
    elif analysis_type == "plaque_assay_plot":
        return plaque_assay_plot(df)  # PFU bar plot
    elif analysis_type == "beta_lactamase_activity_plot":
        return beta_lactamase_activity_plot(df)  # Activity over time
    elif analysis_type == "efflux_pump_activity":
        return efflux_pump_activity_plot(df)  # Accumulation plot
    elif analysis_type == "mdr_index_heatmap":
        return mdr_index_heatmap_plot(df)  # MDR profile heatmap
    elif analysis_type == "susceptibility_profile":
        return susceptibility_profile_plot(df)  # Radar or bar plot
    elif analysis_type == "gel_band_quant_plot":
        return gel_band_quant_plot(df)  # Intensity plot
    # Add more if needed to reach 50; some may be stats-only without plots
    return None, None

# Plot functions

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

def time_kill_curve_plot(df):
    try:
        if 'Time' not in df.columns:
            return "Missing 'Time' column", None
        plt.figure(figsize=(10, 6))
        for col in df.columns[1:]:
            sns.lineplot(data=df, x="Time", y=col, label=col)
        plt.title("Time-Kill Curve")
        plt.xlabel("Time (hours)")
        plt.ylabel("Log CFU/mL")
        plt.yscale('log')
        plot_path, plot_url = save_plot()
        plt.savefig(plot_path)
        plt.close()
        return "Time-kill curve plotted successfully.", plot_url
    except Exception as e:
        return f"Error: {str(e)}", None

def survivorship_curve_plot(df):
    try:
        if 'Time' not in df.columns or 'Survival' not in df.columns:
            return "Missing columns: 'Time' and 'Survival'", None
        plt.figure(figsize=(8, 5))
        sns.lineplot(data=df, x="Time", y="Survival", marker='o')
        plt.title("Survivorship Curve")
        plt.xlabel("Time/Exposure")
        plt.ylabel("% Survival")
        plot_path, plot_url = save_plot()
        plt.savefig(plot_path)
        plt.close()
        return "Survivorship curve plotted successfully.", plot_url
    except Exception as e:
        return f"Error: {str(e)}", None

def ph_change_curve_plot(df):
    try:
        if 'Time' not in df.columns or 'pH' not in df.columns:
            return "Missing columns: 'Time' and 'pH'", None
        plt.figure(figsize=(8, 5))
        sns.lineplot(data=df, x="Time", y="pH", marker='o')
        plt.title("pH Change Over Time")
        plt.xlabel("Time (hours)")
        plt.ylabel("pH")
        plot_path, plot_url = save_plot()
        plt.savefig(plot_path)
        plt.close()
        return "pH change curve plotted successfully.", plot_url
    except Exception as e:
        return f"Error: {str(e)}", None

def pca_microbiome_plot(df):
    try:
        pca = PCA(n_components=2)
        principal_components = pca.fit_transform(df.select_dtypes(include='number'))
        pc_df = pd.DataFrame(data=principal_components, columns=['PC1', 'PC2'])
        plt.figure(figsize=(8, 6))
        sns.scatterplot(x='PC1', y='PC2', data=pc_df)
        plt.title("PCA of Microbiome Data")
        plot_path, plot_url = save_plot()
        plt.savefig(plot_path)
        plt.close()
        return "PCA plot generated successfully.", plot_url
    except Exception as e:
        return f"Error: {str(e)}", None

def rarefaction_curve_plot(df):
    try:
        if 'Depth' not in df.columns or 'Richness' not in df.columns:
            return "Missing columns: 'Depth' and 'Richness'", None
        plt.figure(figsize=(10, 6))
        sns.lineplot(data=df, x="Depth", y="Richness", hue="Sample" if 'Sample' in df.columns else None)
        plt.title("Rarefaction Curve")
        plt.xlabel("Sequencing Depth")
        plt.ylabel("Observed Richness")
        plot_path, plot_url = save_plot()
        plt.savefig(plot_path)
        plt.close()
        return "Rarefaction curve plotted successfully.", plot_url
    except Exception as e:
        return f"Error: {str(e)}", None

def volcano_plot(df):
    try:
        if 'logFC' not in df.columns or 'pvalue' not in df.columns:
            return "Missing columns: 'logFC' and 'pvalue'", None
        df['-logp'] = -np.log10(df["pvalue"])
        plt.figure(figsize=(10, 6))
        sns.scatterplot(data=df, x="logFC", y="-logp")
        plt.title("Volcano Plot for Differential Abundance")
        plt.xlabel("Log Fold Change")
        plt.ylabel("-Log10 p-value")
        plot_path, plot_url = save_plot()
        plt.savefig(plot_path)
        plt.close()
        return "Volcano plot generated successfully.", plot_url
    except Exception as e:
        return f"Error: {str(e)}", None

def gene_expression_heatmap_plot(df):
    try:
        plt.figure(figsize=(10, 8))
        sns.heatmap(df.select_dtypes(include='number'), cmap="RdBu_r", annot=False)
        plt.title("Gene Expression Heatmap")
        plot_path, plot_url = save_plot()
        plt.savefig(plot_path)
        plt.close()
        return "Gene expression heatmap generated successfully.", plot_url
    except Exception as e:
        return f"Error: {str(e)}", None

def cooccurrence_network_plot(df):
    try:
        corr = df.corr()
        plt.figure(figsize=(10, 8))
        sns.heatmap(corr, annot=True, cmap="coolwarm")
        plt.title("Taxa Co-occurrence Correlation Heatmap")
        plot_path, plot_url = save_plot()
        plt.savefig(plot_path)
        plt.close()
        return "Co-occurrence heatmap generated successfully.", plot_url
    except Exception as e:
        return f"Error: {str(e)}", None

# New plot functions

def biofilm_biomass_plot(df):
    try:
        if 'Treatment' not in df.columns or 'Biomass' not in df.columns:
            return "Missing columns: 'Treatment' and 'Biomass'", None
        plt.figure(figsize=(8, 5))
        sns.barplot(data=df, x="Treatment", y="Biomass")
        plt.title("Biofilm Biomass")
        plot_path, plot_url = save_plot()
        plt.savefig(plot_path)
        plt.close()
        return "Biofilm biomass plot generated.", plot_url
    except Exception as e:
        return f"Error: {str(e)}", None

def motility_plot(df):
    try:
        if 'Time' not in df.columns or 'Distance' not in df.columns:
            return "Missing columns: 'Time' and 'Distance'", None
        plt.figure(figsize=(8, 5))
        sns.lineplot(data=df, x="Time", y="Distance")
        plt.title("Motility Assay")
        plot_path, plot_url = save_plot()
        plt.savefig(plot_path)
        plt.close()
        return "Motility plot generated.", plot_url
    except Exception as e:
        return f"Error: {str(e)}", None

def fermentation_gas_plot(df):
    try:
        if 'Time' not in df.columns or 'GasVolume' not in df.columns:
            return "Missing columns: 'Time' and 'GasVolume'", None
        plt.figure(figsize=(8, 5))
        sns.lineplot(data=df, x="Time", y="GasVolume")
        plt.title("Fermentation Gas Production")
        plot_path, plot_url = save_plot()
        plt.savefig(plot_path)
        plt.close()
        return "Fermentation gas plot generated.", plot_url
    except Exception as e:
        return f"Error: {str(e)}", None

def oxygen_consumption_plot(df):
    try:
        if 'Time' not in df.columns or 'O2Level' not in df.columns:
            return "Missing columns: 'Time' and 'O2Level'", None
        plt.figure(figsize=(8, 5))
        sns.lineplot(data=df, x="Time", y="O2Level")
        plt.title("Oxygen Consumption Rate")
        plot_path, plot_url = save_plot()
        plt.savefig(plot_path)
        plt.close()
        return "Oxygen consumption plot generated.", plot_url
    except Exception as e:
        return f"Error: {str(e)}", None

def quorum_sensing_plot(df):
    try:
        if 'Time' not in df.columns or 'Luminescence' not in df.columns:
            return "Missing columns: 'Time' and 'Luminescence'", None
        plt.figure(figsize=(8, 5))
        sns.lineplot(data=df, x="Time", y="Luminescence")
        plt.title("Quorum Sensing Luminescence")
        plot_path, plot_url = save_plot()
        plt.savefig(plot_path)
        plt.close()
        return "Quorum sensing plot generated.", plot_url
    except Exception as e:
        return f"Error: {str(e)}", None

def enzyme_kinetics_plot(df):
    try:
        if 'Substrate' not in df.columns or 'Velocity' not in df.columns:
            return "Missing columns: 'Substrate' and 'Velocity'", None
        plt.figure(figsize=(8, 5))
        sns.lineplot(data=df, x="Substrate", y="Velocity")
        plt.title("Enzyme Kinetics (Michaelis-Menten)")
        plot_path, plot_url = save_plot()
        plt.savefig(plot_path)
        plt.close()
        return "Enzyme kinetics plot generated.", plot_url
    except Exception as e:
        return f"Error: {str(e)}", None

def phylogenetic_tree_plot(df):
    # Assuming df has a 'Newick' column or is a distance matrix
    try:
        if 'Newick' in df.columns:
            tree = Phylo.read(df['Newick'].iloc[0], 'newick')
            plt.figure(figsize=(10, 8))
            Phylo.draw(tree)
            plt.title("Phylogenetic Tree")
            plot_path, plot_url = save_plot()
            plt.savefig(plot_path)
            plt.close()
            return "Phylogenetic tree plotted.", plot_url
        else:
            return "No Newick string found.", None
    except Exception as e:
        return f"Error: {str(e)}", None

def dose_response_curve_plot(df):
    try:
        if 'Concentration' not in df.columns or 'Response' not in df.columns:
            return "Missing columns: 'Concentration' and 'Response'", None
        plt.figure(figsize=(8, 5))
        sns.lineplot(data=df, x="Concentration", y="Response", marker='o')
        plt.xscale('log')
        plt.title("Dose-Response Curve")
        plot_path, plot_url = save_plot()
        plt.savefig(plot_path)
        plt.close()
        return "Dose-response curve plotted.", plot_url
    except Exception as e:
        return f"Error: {str(e)}", None

def synergy_heatmap_plot(df):
    try:
        plt.figure(figsize=(10, 8))
        sns.heatmap(df, cmap="coolwarm", annot=True)
        plt.title("Drug Synergy Heatmap (FICI)")
        plot_path, plot_url = save_plot()
        plt.savefig(plot_path)
        plt.close()
        return "Synergy heatmap generated.", plot_url
    except Exception as e:
        return f"Error: {str(e)}", None

def flow_cytometry_histogram_plot(df):
    try:
        if 'Fluorescence' not in df.columns:
            return "Missing 'Fluorescence' column", None
        plt.figure(figsize=(8, 5))
        sns.histplot(df['Fluorescence'], bins=50)
        plt.title("Flow Cytometry Histogram")
        plot_path, plot_url = save_plot()
        plt.savefig(plot_path)
        plt.close()
        return "Flow cytometry histogram plotted.", plot_url
    except Exception as e:
        return f"Error: {str(e)}", None

def qpcr_ct_plot(df):
    try:
        if 'Sample' not in df.columns or 'Ct' not in df.columns:
            return "Missing columns: 'Sample' and 'Ct'", None
        plt.figure(figsize=(8, 5))
        sns.barplot(data=df, x="Sample", y="Ct")
        plt.title("qPCR Cycle Threshold")
        plot_path, plot_url = save_plot()
        plt.savefig(plot_path)
        plt.close()
        return "qPCR Ct plot generated.", plot_url
    except Exception as e:
        return f"Error: {str(e)}", None

def western_blot_densitometry_plot(df):
    try:
        if 'Band' not in df.columns or 'Intensity' not in df.columns:
            return "Missing columns: 'Band' and 'Intensity'", None
        plt.figure(figsize=(8, 5))
        sns.barplot(data=df, x="Band", y="Intensity")
        plt.title("Western Blot Densitometry")
        plot_path, plot_url = save_plot()
        plt.savefig(plot_path)
        plt.close()
        return "Western blot densitometry plotted.", plot_url
    except Exception as e:
        return f"Error: {str(e)}", None

def elisa_standard_curve_plot(df):
    try:
        if 'Concentration' not in df.columns or 'Absorbance' not in df.columns:
            return "Missing columns: 'Concentration' and 'Absorbance'", None
        plt.figure(figsize=(8, 5))
        sns.lineplot(data=df, x="Concentration", y="Absorbance")
        plt.title("ELISA Standard Curve")
        plot_path, plot_url = save_plot()
        plt.savefig(plot_path)
        plt.close()
        return "ELISA standard curve plotted.", plot_url
    except Exception as e:
        return f"Error: {str(e)}", None

def hemagglutination_titer_plot(df):
    try:
        if 'Dilution' not in df.columns or 'Titer' not in df.columns:
            return "Missing columns: 'Dilution' and 'Titer'", None
        plt.figure(figsize=(8, 5))
        sns.barplot(data=df, x="Dilution", y="Titer")
        plt.title("Hemagglutination Titer")
        plot_path, plot_url = save_plot()
        plt.savefig(plot_path)
        plt.close()
        return "Hemagglutination titer plotted.", plot_url
    except Exception as e:
        return f"Error: {str(e)}", None

def plaque_assay_plot(df):
    try:
        if 'Sample' not in df.columns or 'PFU' not in df.columns:
            return "Missing columns: 'Sample' and 'PFU'", None
        plt.figure(figsize=(8, 5))
        sns.barplot(data=df, x="Sample", y="PFU")
        plt.title("Plaque Assay PFU")
        plot_path, plot_url = save_plot()
        plt.savefig(plot_path)
        plt.close()
        return "Plaque assay plot generated.", plot_url
    except Exception as e:
        return f"Error: {str(e)}", None

def beta_lactamase_activity_plot(df):
    try:
        if 'Time' not in df.columns or 'Activity' not in df.columns:
            return "Missing columns: 'Time' and 'Activity'", None
        plt.figure(figsize=(8, 5))
        sns.lineplot(data=df, x="Time", y="Activity")
        plt.title("Beta-Lactamase Activity")
        plot_path, plot_url = save_plot()
        plt.savefig(plot_path)
        plt.close()
        return "Beta-lactamase activity plot generated.", plot_url
    except Exception as e:
        return f"Error: {str(e)}", None

def efflux_pump_activity_plot(df):
    try:
        if 'Time' not in df.columns or 'Accumulation' not in df.columns:
            return "Missing columns: 'Time' and 'Accumulation'", None
        plt.figure(figsize=(8, 5))
        sns.lineplot(data=df, x="Time", y="Accumulation")
        plt.title("Efflux Pump Activity")
        plot_path, plot_url = save_plot()
        plt.savefig(plot_path)
        plt.close()
        return "Efflux pump activity plot generated.", plot_url
    except Exception as e:
        return f"Error: {str(e)}", None

def mdr_index_heatmap_plot(df):
    try:
        plt.figure(figsize=(10, 8))
        sns.heatmap(df, cmap="viridis", annot=True)
        plt.title("Multi-Drug Resistance Index Heatmap")
        plot_path, plot_url = save_plot()
        plt.savefig(plot_path)
        plt.close()
        return "MDR index heatmap generated.", plot_url
    except Exception as e:
        return f"Error: {str(e)}", None

def susceptibility_profile_plot(df):
    try:
        if 'Antibiotic' not in df.columns or 'MIC' not in df.columns:
            return "Missing columns: 'Antibiotic' and 'MIC'", None
        plt.figure(figsize=(8, 5))
        sns.barplot(data=df, x="Antibiotic", y="MIC")
        plt.title("Antibiotic Susceptibility Profile")
        plot_path, plot_url = save_plot()
        plt.savefig(plot_path)
        plt.close()
        return "Susceptibility profile plotted.", plot_url
    except Exception as e:
        return f"Error: {str(e)}", None

def gel_band_quant_plot(df):
    try:
        if 'Band' not in df.columns or 'Quantity' not in df.columns:
            return "Missing columns: 'Band' and 'Quantity'", None
        plt.figure(figsize=(8, 5))
        sns.barplot(data=df, x="Band", y="Quantity")
        plt.title("Gel Band Quantification")
        plot_path, plot_url = save_plot()
        plt.savefig(plot_path)
        plt.close()
        return "Gel band quantification plotted.", plot_url
    except Exception as e:
        return f"Error: {str(e)}", None

# Stats functions (existing and new to reach 35-50 total analyses)

def growth_stats(df):
    try:
        if 'Time' not in df.columns or 'OD' not in df.columns:
            return {"Error": "Missing 'Time' and 'OD' columns"}
        df = df.sort_values("Time")
        df["logOD"] = np.log(df["OD"] + 1e-10)  # Avoid log(0)
        slope, _, _, _, _ = stats.linregress(df["Time"], df["logOD"])
        doubling_time = np.log(2) / slope if slope > 0 else np.nan
        return {
            "Growth Rate (1/hr)": round(slope, 4),
            "Doubling Time (hr)": round(doubling_time, 2)
        }
    except Exception as e:
        return {"Error": str(e)}

def inhibition_stats(df):
    try:
        antibiotics = df.columns[1:]
        means = df[antibiotics].mean()
        stds = df[antibiotics].std()
        stats_dict = {f"{ab} Mean (mm)": round(means[ab], 2) for ab in antibiotics}
        if len(antibiotics) >= 2:
            t_stat, p_value = stats.ttest_ind(df[antibiotics[0]], df[antibiotics[1]], nan_policy='omit')
            ttest_dict = {
                "Antibiotic 1": antibiotics[0],
                "Antibiotic 2": antibiotics[1],
                "t-statistic": round(t_stat, 4),
                "p-value": round(p_value, 4)
            }
        else:
            ttest_dict = {"Note": "Need at least 2 antibiotic columns"}
        return stats_dict, ttest_dict
    except Exception as e:
        return {"Error": str(e)}, {}

def cfu_stats(df):
    try:
        return df.describe()
    except Exception as e:
        return pd.DataFrame({"Error": [str(e)]})

def compute_diversity(df, method="shannon"):
    counts = df["count"].values
    proportions = counts / counts.sum() if counts.sum() > 0 else np.zeros_like(counts)
    if method == "shannon":
        index = entropy(proportions)
        return pd.DataFrame({"Index": ["Shannon"], "Value": [round(index, 4)]})
    elif method == "simpson":
        simpson = 1 - sum(proportions**2)
        return pd.DataFrame({"Index": ["Simpson"], "Value": [round(simpson, 4)]})
    else:
        raise ValueError("Invalid method selected.")

def mic_determination(df):
    try:
        if 'Concentration' not in df.columns or 'Growth' not in df.columns:
            return {"Error": "Missing 'Concentration' and 'Growth' columns"}
        df = df.sort_values("Concentration")
        mic = df[df['Growth'] <= 0]['Concentration'].min()
        return {"MIC": mic if not pd.isna(mic) else "Not determined"}
    except Exception as e:
        return {"Error": str(e)}

def mbc_determination(df):
    try:
        if 'Concentration' not in df.columns or 'Viability' not in df.columns:
            return {"Error": "Missing 'Concentration' and 'Viability' columns"}
        df = df.sort_values("Concentration")
        mbc = df[df['Viability'] <= 0.01]['Concentration'].min()  # 99% kill
        return {"MBC": mbc if not pd.isna(mbc) else "Not determined"}
    except Exception as e:
        return {"Error": str(e)}

def d_value_calc(df):
    try:
        if 'Time' not in df.columns or 'LogCFU' not in df.columns:
            return {"Error": "Missing 'Time' and 'LogCFU' columns"}
        slope, _, _, _, _ = stats.linregress(df["Time"], df["LogCFU"])
        d_value = -1 / slope if slope < 0 else np.nan
        return {"D-Value (min)": round(d_value, 2)}
    except Exception as e:
        return {"Error": str(e)}

def z_value_calc(df):
    try:
        if 'Temp' not in df.columns or 'DValue' not in df.columns:
            return {"Error": "Missing 'Temp' and 'DValue' columns"}
        df['LogD'] = np.log10(df['DValue'] + 1e-10)
        slope, _, _, _, _ = stats.linregress(df["Temp"], df['LogD'])
        z_value = -1 / slope if slope < 0 else np.nan
        return {"Z-Value (°C)": round(z_value, 2)}
    except Exception as e:
        return {"Error": str(e)}

def f_value_calc(df):
    try:
        if 'Time' not in df.columns or 'Temp' not in df.columns:
            return {"Error": "Missing 'Time' and 'Temp' columns"}
        ref_temp = 121.1  # Standard for F-value
        z = 10  # Default Z-value
        f_value = sum(10**((df['Temp'] - ref_temp)/z) * np.diff(df['Time'].append(pd.Series([0])))[:-1])
        return {"F-Value": round(f_value, 2)}
    except Exception as e:
        return {"Error": str(e)}

def lag_phase_est(df):
    try:
        if 'Time' not in df.columns or 'OD' not in df.columns:
            return {"Error": "Missing 'Time' and 'OD' columns"}
        df['LogOD'] = np.log(df['OD'] + 1e-10)
        slope, intercept, _, _, _ = stats.linregress(df["Time"], df['LogOD'])
        lag = (np.log(df['OD'].iloc[0] + 1e-10) - intercept) / slope
        return {"Lag Phase (hr)": round(lag, 2) if lag > 0 else 0}
    except Exception as e:
        return {"Error": str(e)}

def mu_max_calc(df):
    try:
        if 'Time' not in df.columns or 'OD' not in df.columns:
            return {"Error": "Missing 'Time' and 'OD' columns"}
        df = df.sort_values("Time")
        spline = UnivariateSpline(df['Time'], np.log(df['OD'] + 1e-10), s=0)
        mu_max = np.max(spline.derivative()(df['Time']))
        return {"Mu Max (1/hr)": round(mu_max, 4)}
    except Exception as e:
        return {"Error": str(e)}

def baranyi_model_fit(df):
    def baranyi(t, y0, mu, lag, ymax):
        a = mu * (t + (1/mu) * np.log(np.exp(-mu * t) + np.exp(-mu * lag) - np.exp(-mu * (t + lag))))
        return y0 + a - np.log(1 + (np.exp(a) - 1) / np.exp(ymax - y0))
    try:
        if 'Time' not in df.columns or 'LogCFU' not in df.columns:
            return {"Error": "Missing 'Time' and 'LogCFU' columns"}
        popt, _ = curve_fit(baranyi, df['Time'], df['LogCFU'], p0=[df['LogCFU'].min(), 0.1, 1, df['LogCFU'].max()])
        return {"y0": round(popt[0], 2), "mu": round(popt[1], 4), "lag": round(popt[2], 2), "ymax": round(popt[3], 2)}
    except Exception as e:
        return {"Error": str(e)}

def gompertz_model_fit(df):
    def gompertz(t, a, b, c):
        return a * np.exp(-np.exp(b - c * t))
    try:
        if 'Time' not in df.columns or 'OD' not in df.columns:
            return {"Error": "Missing 'Time' and 'OD' columns"}
        popt, _ = curve_fit(gompertz, df['Time'], df['OD'], p0=[df['OD'].max(), 1, 0.1])
        return {"A": round(popt[0], 2), "B": round(popt[1], 2), "C": round(popt[2], 4)}
    except Exception as e:
        return {"Error": str(e)}

def logistic_model_fit(df):
    def logistic(t, k, n0, r):
        return k / (1 + ((k - n0)/n0) * np.exp(-r * t))
    try:
        if 'Time' not in df.columns or 'Population' not in df.columns:
            return {"Error": "Missing 'Time' and 'Population' columns"}
        popt, _ = curve_fit(logistic, df['Time'], df['Population'], p0=[df['Population'].max(), df['Population'].min(), 0.1])
        return {"K": round(popt[0], 2), "N0": round(popt[1], 2), "r": round(popt[2], 4)}
    except Exception as e:
        return {"Error": str(e)}

def alpha_diversity(df):
    try:
        counts = df["count"].values
        richness = len(counts[counts > 0])
        shannon = entropy(counts / counts.sum() if counts.sum() > 0 else 1)
        evenness = shannon / np.log(richness) if richness > 1 else 0
        return pd.DataFrame({"Metric": ["Richness", "Shannon", "Evenness"], "Value": [richness, round(shannon, 4), round(evenness, 4)]})
    except Exception as e:
        return pd.DataFrame({"Error": [str(e)]})

def beta_diversity(df):
    # Assume df has samples as rows, taxa as columns
    try:
        dist = pairwise_distances(df, metric='braycurtis')
        return pd.DataFrame(dist, index=df.index, columns=df.index)
    except Exception as e:
        return pd.DataFrame({"Error": [str(e)]})

def phylogenetic_diversity(df):
    # Assume df is distance matrix or tree
    try:
        # Simple sum of branch lengths if tree
        if 'Newick' in df.columns:
            tree = Phylo.read(df['Newick'].iloc[0], 'newick')
            pd = sum([branch.length for branch in tree.get_terminals() + tree.get_nonterminals() if branch.length])
            return {"Phylogenetic Diversity": round(pd, 2)}
        else:
            return {"Error": "No tree data"}
    except Exception as e:
        return {"Error": str(e)}

def mutation_frequency(df):
    try:
        if 'Mutants' not in df.columns or 'Total' not in df.columns:
            return {"Error": "Missing 'Mutants' and 'Total' columns"}
        freq = (df['Mutants'] / df['Total']).mean()
        return {"Mutation Frequency": round(freq, 6)}
    except Exception as e:
        return {"Error": str(e)}

def pfu_calc(df):
    try:
        if 'Plaques' not in df.columns or 'Dilution' not in df.columns or 'Volume' not in df.columns:
            return {"Error": "Missing required columns"}
        pfu = (df['Plaques'] * df['Dilution']) / df['Volume']
        df['PFU'] = pfu
        return df.describe()
    except Exception as e:
        return pd.DataFrame({"Error": [str(e)]})

def synergy_testing(df):
    # Assume columns: DrugA_MIC, DrugB_MIC, Combo_MIC_A, Combo_MIC_B
    try:
        fici = (df['Combo_MIC_A'] / df['DrugA_MIC']) + (df['Combo_MIC_B'] / df['DrugB_MIC'])
        return {"FICI Mean": round(fici.mean(), 4), "Interpretation": "Synergy" if fici.mean() < 0.5 else "Additive" if fici.mean() < 1 else "Antagonism"}
    except Exception as e:
        return {"Error": str(e)}

def post_antibiotic_effect(df):
    # Assume columns: Time, Control_Growth, Treated_Growth
    try:
        pae = df['Time'][df['Treated_Growth'] > df['Control_Growth'].iloc[0]].min() - df['Time'][df['Control_Growth'] > df['Control_Growth'].iloc[0]].min()
        return {"PAE (hr)": round(pae, 2)}
    except Exception as e:
        return {"Error": str(e)}

def biofilm_biomass_stats(df):
    try:
        if 'Biomass' not in df.columns:
            return {"Error": "Missing 'Biomass' column"}
        return df['Biomass'].describe()
    except Exception as e:
        return pd.DataFrame({"Error": [str(e)]})

def biofilm_viability(df):
    try:
        if 'Live' not in df.columns or 'Dead' not in df.columns:
            return {"Error": "Missing 'Live' and 'Dead' columns"}
        viability = df['Live'] / (df['Live'] + df['Dead'])
        return {"Viability Ratio": round(viability.mean(), 4)}
    except Exception as e:
        return {"Error": str(e)}

def motility_assay_stats(df):
    try:
        if 'Distance' not in df.columns:
            return {"Error": "Missing 'Distance' column"}
        return {"Mean Motility (mm)": round(df['Distance'].mean(), 2)}
    except Exception as e:
        return {"Error": str(e)}

def chemotaxis_index(df):
    try:
        if 'Towards' not in df.columns or 'Away' not in df.columns:
            return {"Error": "Missing 'Towards' and 'Away' columns"}
        ci = df['Towards'] / df['Away']
        return {"Chemotaxis Index": round(ci.mean(), 4)}
    except Exception as e:
        return {"Error": str(e)}

def spore_count(df):
    try:
        if 'Spores' not in df.columns or 'Total' not in df.columns:
            return {"Error": "Missing 'Spores' and 'Total' columns"}
        percent = (df['Spores'] / df['Total']) * 100
        return {"Spore Percentage": round(percent.mean(), 2)}
    except Exception as e:
        return {"Error": str(e)}

def viability_staining(df):
    try:
        if 'Live' not in df.columns or 'Dead' not in df.columns:
            return {"Error": "Missing 'Live' and 'Dead' columns"}
        ratio = df['Live'] / df['Dead']
        return {"Live/Dead Ratio": round(ratio.mean(), 4)}
    except Exception as e:
        return {"Error": str(e)}

def flow_cytometry_analysis(df):
    try:
        if 'Fluorescence' not in df.columns:
            return {"Error": "Missing 'Fluorescence' column"}
        return df['Fluorescence'].describe()
    except Exception as e:
        return pd.DataFrame({"Error": [str(e)]})

def qpcr_relative_quant(df):
    try:
        if 'Ct_Target' not in df.columns or 'Ct_Reference' not in df.columns:
            return {"Error": "Missing 'Ct_Target' and 'Ct_Reference' columns"}
        delta_ct = df['Ct_Target'] - df['Ct_Reference']
        rq = 2 ** -delta_ct
        return {"Relative Quantification Mean": round(rq.mean(), 4)}
    except Exception as e:
        return {"Error": str(e)}

def qpcr_absolute_quant(df):
    try:
        if 'Ct' not in df.columns or 'Standard_Curve_Slope' not in df.columns:
            return {"Error": "Missing required columns"}
        # Assume standard curve provided
        copies = 10 ** ((df['Ct'] - df['Intercept']) / df['Standard_Curve_Slope'])
        return {"Absolute Copies Mean": round(copies.mean(), 2)}
    except Exception as e:
        return {"Error": str(e)}

def enzyme_activity(df):
    try:
        if 'Time' not in df.columns or 'Product' not in df.columns:
            return {"Error": "Missing 'Time' and 'Product' columns"}
        rate = np.diff(df['Product']) / np.diff(df['Time'])
        return {"Enzyme Activity Rate": round(rate.mean(), 4)}
    except Exception as e:
        return {"Error": str(e)}

def protein_conc(df):
    try:
        if 'Absorbance' not in df.columns:
            return {"Error": "Missing 'Absorbance' column"}
        # Assume Bradford or BCA standard
        conc = df['Absorbance'] * 1.0  # Placeholder factor
        return {"Protein Concentration Mean (mg/mL)": round(conc.mean(), 2)}
    except Exception as e:
        return {"Error": str(e)}

def dna_rna_conc(df):
    try:
        if 'A260' not in df.columns:
            return {"Error": "Missing 'A260' column"}
        conc = df['A260'] * 50  # For DNA, ng/uL
        return {"DNA/RNA Concentration Mean (ng/uL)": round(conc.mean(), 2)}
    except Exception as e:
        return {"Error": str(e)}

def disk_diffusion_interpret(df):
    try:
        if 'Zone' not in df.columns or 'Antibiotic' not in df.columns:
            return {"Error": "Missing 'Zone' and 'Antibiotic' columns"}
        # Placeholder interpretation
        interpret = df.apply(lambda row: "Sensitive" if row['Zone'] > 15 else "Resistant", axis=1)
        return pd.DataFrame({"Antibiotic": df['Antibiotic'], "Interpretation": interpret})
    except Exception as e:
        return pd.DataFrame({"Error": [str(e)]})

def mdr_index(df):
    try:
        if 'Resistant' not in df.columns or 'Total_AB' not in df.columns:
            return {"Error": "Missing 'Resistant' and 'Total_AB' columns"}
        index = df['Resistant'] / df['Total_AB']
        return {"MDR Index Mean": round(index.mean(), 4)}
    except Exception as e:
        return {"Error": str(e)}

def efflux_pump_stats(df):
    try:
        if 'Accumulation' not in df.columns:
            return {"Error": "Missing 'Accumulation' column"}
        return {"Mean Accumulation": round(df['Accumulation'].mean(), 2)}
    except Exception as e:
        return {"Error": str(e)}

def beta_lactamase_stats(df):
    try:
        if 'Activity' not in df.columns:
            return {"Error": "Missing 'Activity' column"}
        return {"Mean Activity": round(df['Activity'].mean(), 4)}
    except Exception as e:
        return {"Error": str(e)}

def hemagglutination_titer(df):
    try:
        if 'Dilution' not in df.columns or 'Agglutination' not in df.columns:
            return {"Error": "Missing 'Dilution' and 'Agglutination' columns"}
        titer = df[df['Agglutination'] > 0]['Dilution'].max()
        return {"Titer": titer}
    except Exception as e:
        return {"Error": str(e)}

def elisa_quant(df):
    try:
        if 'Absorbance' not in df.columns or 'Standard_Conc' not in df.columns:
            return {"Error": "Missing required columns"}
        # Simple linear fit for standard
        slope, intercept, _, _, _ = stats.linregress(df['Standard_Conc'], df['Absorbance'])
        quant = (df['Absorbance'] - intercept) / slope
        return {"Quantification Mean": round(quant.mean(), 2)}
    except Exception as e:
        return {"Error": str(e)}

def western_blot_quant(df):
    try:
        if 'Intensity' not in df.columns:
            return {"Error": "Missing 'Intensity' column"}
        return {"Mean Intensity": round(df['Intensity'].mean(), 2)}
    except Exception as e:
        return {"Error": str(e)}

def gel_band_quant(df):
    try:
        if 'Quantity' not in df.columns:
            return {"Error": "Missing 'Quantity' column"}
        return {"Mean Quantity": round(df['Quantity'].mean(), 2)}
    except Exception as e:
        return {"Error": str(e)}

def susceptibility_profile_stats(df):
    try:
        if 'MIC' not in df.columns:
            return {"Error": "Missing 'MIC' column"}
        return df['MIC'].describe()
    except Exception as e:
        return pd.DataFrame({"Error": [str(e)]})



# Manual calculation functions

def perform_growth_rate(initial_pop, final_pop, time):
    try:
        growth_rate = (math.log(final_pop) - math.log(initial_pop)) / time if time > 0 else 0
        return f"Growth Rate (μ): {growth_rate:.4f} per hour"
    except Exception as e:
        return f"Error in growth rate calculation: {str(e)}"

def perform_cfu(colony_count, dilution_factor, volume_plated):
    try:
        cfu = (colony_count * dilution_factor) / volume_plated if volume_plated > 0 else 0
        return f"Colony Forming Units (CFU/mL): {cfu:.2e}"
    except Exception as e:
        return f"Error in CFU calculation: {str(e)}"

def perform_zone_of_inhibition(diameter):
    try:
        radius = diameter / 2
        area = math.pi * (radius ** 2)
        return f"Zone of Inhibition Area: {area:.2f} mm²"
    except Exception as e:
        return f"Error in zone of inhibition calculation: {str(e)}"

def perform_diversity_index(counts, method="shannon"):
    try:
        total = sum(counts)
        if total == 0:
            return "Total count is zero."
        if method == "shannon":
            proportions = [c / total for c in counts if c > 0]
            shannon_index = -sum(p * math.log(p) for p in proportions)
            return f"Shannon Diversity Index (H'): {shannon_index:.4f}"
        elif method == "simpson":
            proportions = [c / total for c in counts]
            simpson_index = 1 - sum(p ** 2 for p in proportions)
            return f"Simpson Diversity Index (1 - D): {simpson_index:.4f}"
        else:
            return "Unknown method."
    except Exception as e:
        return f"Error: {str(e)}"


def calculate_reed_muench(dilutions, infected, total):
    try:
        prop_infected = [i / t if t != 0 else 0 for i, t in zip(infected, total)]
        cum_above = [sum(prop_infected[j:] ) for j in range(len(prop_infected))]
        cum_below = [sum(prop_infected[:j+1]) for j in range(len(prop_infected))]
        diffs = [abs(a - 0.5) for a in cum_above]  # Simplified
        min_index = diffs.index(min(diffs))
        if min_index + 1 < len(dilutions):
            d1, d2 = dilutions[min_index], dilutions[min_index + 1]
            r1, r2 = prop_infected[min_index], prop_infected[min_index + 1]
            tcid50 = d1 - ( (r1 - 0.5) / (r1 - r2) ) * (d1 - d2)
            return round(tcid50, 4)
        else:
            return "Calculation not possible."
    except Exception as e:
        return f"Error: {str(e)}"

def perform_mic_manual(concentrations, growths):
    try:
        mic = min([conc for conc, growth in zip(concentrations, growths) if growth <= 0], default="Not determined")
        return f"MIC: {mic}"
    except Exception as e:
        return f"Error: {str(e)}"

def perform_d_value_manual(times, log_cfus):
    try:
        slope, _, _, _, _ = stats.linregress(times, log_cfus)
        d_value = -1 / slope if slope < 0 else "Invalid"
        return f"D-Value: {round(d_value, 2)} min"
    except Exception as e:
        return f"Error: {str(e)}"

def perform_z_value_manual(temps, d_values):
    try:
        log_d = np.log10(d_values)
        slope, _, _, _, _ = stats.linregress(temps, log_d)
        z_value = -1 / slope if slope < 0 else "Invalid"
        return f"Z-Value: {round(z_value, 2)} °C"
    except Exception as e:
        return f"Error: {str(e)}"

def perform_lag_phase_manual(times, ods):
    try:
        log_od = np.log(ods + 1e-10)
        slope, intercept, _, _, _ = stats.linregress(times, log_od)
        lag = (log_od[0] - intercept) / slope
        return f"Lag Phase: {round(lag, 2)} hr"
    except Exception as e:
        return f"Error: {str(e)}"

def perform_mu_max_manual(times, ods):
    try:
        spline = UnivariateSpline(times, np.log(ods + 1e-10), s=0)
        mu_max = np.max(spline.derivative()(times))
        return f"Mu Max: {round(mu_max, 4)} 1/hr"
    except Exception as e:
        return f"Error: {str(e)}"

def perform_mbc_manual(concentrations, viabilities):
    try:
        mbc = min([conc for conc, v in zip(concentrations, viabilities) if v <= 0.01], default="Not determined")
        return f"MBC: {mbc}"
    except Exception as e:
        return f"Error: {str(e)}"

def perform_f_value_manual(times, temps, ref_temp=121.1, z=10):
    try:
        dt = np.diff(times, prepend=0)
        f = sum(dt * 10**((temps - ref_temp)/z))
        return f"F-Value: {round(f, 2)}"
    except Exception as e:
        return f"Error: {str(e)}"

def calculate_fici(drug_a_mic, drug_b_mic, combo_mic_a, combo_mic_b):
    try:
        fica = combo_mic_a / drug_a_mic if drug_a_mic != 0 else float('inf')
        ficb = combo_mic_b / drug_b_mic if drug_b_mic != 0 else float('inf')
        fici = fica + ficb
        interpretation = (
            "Synergy" if fici <= 0.5 else
            "Additive" if 0.5 < fici <= 1 else
            "Indifferent" if 1 < fici <= 4 else
            "Antagonism"
        )
        return f"FICI: {round(fici, 2)} ({interpretation})"
    except Exception as e:
        return f"Error: {str(e)}"

def calculate_mutation_frequency(mutants, total):
    try:
        if total == 0:
            return "Error: Total cannot be zero"
        frequency = mutants / total
        return f"Mutation Frequency: {frequency:.2e}"
    except Exception as e:
        return f"Error: {str(e)}"

def calculate_pfu(plaques, dilution, volume):
    try:
        if volume == 0:
            return "Error: Volume cannot be zero"
        pfu = (plaques * dilution) / volume
        return f"PFU: {pfu:.2e} PFU/mL"
    except Exception as e:
        return f"Error: {str(e)}"

def calculate_spore_percentage(spores, total):
    try:
        if total == 0:
            return "Error: Total cannot be zero"
        percentage = (spores / total) * 100
        return f"Spore Percentage: {round(percentage, 2)}%"
    except Exception as e:
        return f"Error: {str(e)}"

def calculate_live_dead_ratio(live, dead):
    try:
        total = live + dead
        if total == 0:
            return "Error: Total cells cannot be zero"
        ratio = live / dead if dead != 0 else float('inf')
        percentage = (live / total) * 100
        return f"Live/Dead Ratio: {round(ratio, 2)}, Live Percentage: {round(percentage, 2)}%"
    except Exception as e:
        return f"Error: {str(e)}"

def calculate_enzyme_activity(times, products):
    try:
        times = [float(t) for t in times.split(',')]
        products = [float(p) for p in products.split(',')]
        if len(times) != len(products):
            return "Error: Times and products lists must have equal length"
        slope, _, _, _, _ = stats.linregress(times, products)
        return f"Enzyme Activity: {round(slope, 4)} units/min"
    except Exception as e:
        return f"Error: {str(e)}"

def calculate_protein_concentration(absorbance, factor):
    try:
        concentration = absorbance * factor
        return f"Protein Concentration: {round(concentration, 2)} µg/mL"
    except Exception as e:
        return f"Error: {str(e)}"

def calculate_dna_rna_concentration(a260, factor):
    try:
        concentration = a260 * factor
        return f"DNA/RNA Concentration: {round(concentration, 2)} µg/mL"
    except Exception as e:
        return f"Error: {str(e)}"

def calculate_hemagglutination_titer(dilutions, agglutinations):
    try:
        dilutions = [float(d) for d in dilutions.split(',')]
        agglutinations = [int(a) for a in agglutinations.split(',')]
        if len(dilutions) != len(agglutinations):
            return "Error: Dilutions and agglutinations lists must have equal length"
        for i, agg in enumerate(agglutinations):
            if agg == 0:
                return f"Titer: {dilutions[i-1]}" if i > 0 else "Titer: Not determined"
        return f"Titer: {max(dilutions)}"
    except Exception as e:
        return f"Error: {str(e)}"

def calculate_elisa_quantification(absorbances, standard_concs):
    try:
        absorbances = [float(a) for a in absorbances.split(',')]
        standard_concs = [float(c) for c in standard_concs.split(',')]
        if len(absorbances) != len(standard_concs):
            return "Error: Absorbances and standard concentrations lists must have equal length"
        slope, _, _, _, _ = stats.linregress(standard_concs, absorbances)
        return f"ELISA Quantification Slope: {round(slope, 4)} absorbance units/concentration"
    except Exception as e:
        return f"Error: {str(e)}"