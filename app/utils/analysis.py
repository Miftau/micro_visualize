import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from io import BytesIO

def plot_growth_curve(df):
    fig, ax = plt.subplots()
    if 'Time' in df.columns and 'OD' in df.columns:
        sns.lineplot(x='Time', y='OD', data=df, ax=ax)
        ax.set_title("Growth Curve")
        ax.set_xlabel("Time (hours)")
        ax.set_ylabel("Optical Density (OD600)")
    else:
        ax.text(0.5, 0.5, "Missing 'Time' or 'OD' columns", ha='center')
    return fig

def plot_microbiome_composition(df):
    fig, ax = plt.subplots()
    if 'Taxon' in df.columns and 'Abundance' in df.columns:
        sns.barplot(x='Taxon', y='Abundance', data=df, ax=ax)
        ax.set_title("Microbiome Composition")
        ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right')
    else:
        ax.text(0.5, 0.5, "Missing 'Taxon' or 'Abundance' columns", ha='center')
    return fig

def plot_amr_heatmap(df):
    fig, ax = plt.subplots()
    sns.heatmap(df.corr(), cmap='coolwarm', ax=ax)
    ax.set_title("AMR Heatmap")
    return fig

def plot_gene_expression(df):
    fig, ax = plt.subplots()
    if 'Gene' in df.columns and 'Expression' in df.columns:
        sns.violinplot(x='Gene', y='Expression', data=df, ax=ax)
        ax.set_title("Gene Expression Distribution")
    else:
        ax.text(0.5, 0.5, "Missing 'Gene' or 'Expression' columns", ha='center')
    return fig
