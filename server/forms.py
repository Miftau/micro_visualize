# forms.py (Fully Extended)

from flask_wtf import FlaskForm
from wtforms import (
    SubmitField, DecimalField, IntegerField, SelectField, FloatField, TextAreaField, FileField
)
from wtforms.validators import DataRequired

class UploadForm(FlaskForm):
    file = FileField("Upload Excel or CSV File", validators=[DataRequired()])
    analysis_type = SelectField("Analysis Type", choices=[
        ("growth_curve", "Bacterial Growth Curve"),
        ("amr_heatmap", "Antimicrobial Resistance Heatmap"),
        ("amr_prediction", "AMR Prediction (ML)"),  # Assuming this exists
        ("growth_stats", "Growth Rate & Doubling Time"),
        ("inhibition_stats", "Zone of Inhibition Stats"),
        ("cfu_stats", "CFU Count Statistics"),
        ("diversity_stats", "Microbial Diversity Indices"),
        ("time_kill_curve", "Time-Kill Curve"),
        ("survivorship_curve", "Survivorship Curve Plot"),
        ("ph_change_curve", "pH Change Curve"),
        ("pca_microbiome", "PCA for Microbiome"),
        ("rarefaction_curve", "Rarefaction Curve"),
        ("volcano_plot", "Volcano Plot"),
        ("gene_expression_heatmap", "Gene Expression Heatmap"),
        ("cooccurrence_network", "Co-occurrence Network"),
        ("biofilm_biomass_plot", "Biofilm Biomass Plot"),
        ("motility_plot", "Motility Plot"),
        ("fermentation_gas_plot", "Fermentation Gas Plot"),
        ("oxygen_consumption_plot", "Oxygen Consumption Plot"),
        ("quorum_sensing_plot", "Quorum Sensing Plot"),
        ("enzyme_kinetics_plot", "Enzyme Kinetics Plot"),
        ("phylogenetic_tree", "Phylogenetic Tree"),
        ("dose_response_curve", "Dose Response Curve"),
        ("synergy_heatmap", "Synergy Heatmap"),
        ("flow_cytometry_histogram", "Flow Cytometry Histogram"),
        ("qpcr_ct_plot", "qPCR Ct Plot"),
        ("western_blot_densitometry", "Western Blot Densitometry Plot"),
        ("elisa_standard_curve", "ELISA Standard Curve"),
        ("hemagglutination_titer_plot", "Hemagglutination Titer Plot"),
        ("plaque_assay_plot", "Plaque Assay Plot"),
        ("beta_lactamase_activity_plot", "Beta Lactamase Activity Plot"),
        ("efflux_pump_activity", "Efflux Pump Activity Plot"),
        ("mdr_index_heatmap", "MDR Index Heatmap"),
        ("susceptibility_profile", "Susceptibility Profile Plot"),
        ("gel_band_quant_plot", "Gel Band Quant Plot"),
        ("mic_determination", "MIC Determination"),
        ("mbc_determination", "MBC Determination"),
        ("d_value_calc", "D-Value Calculation"),
        ("z_value_calc", "Z-Value Calculation"),
        ("f_value_calc", "F-Value Calculation"),
        ("lag_phase_est", "Lag Phase Estimation"),
        ("mu_max_calc", "Maximum Specific Growth Rate"),
        ("baranyi_model_fit", "Baranyi Model Fit"),
        ("gompertz_model_fit", "Gompertz Model Fit"),
        ("logistic_model_fit", "Logistic Model Fit"),
        ("alpha_diversity", "Alpha Diversity Metrics"),
        ("beta_diversity", "Beta Diversity"),
        ("phylogenetic_diversity", "Phylogenetic Diversity"),
        ("mutation_frequency", "Mutation Frequency"),
        ("pfu_calc", "Plaque Forming Units Stats"),
        ("synergy_testing", "Synergy Testing"),
        ("post_antibiotic_effect", "Post-Antibiotic Effect"),
        ("biofilm_biomass_stats", "Biofilm Biomass Stats"),
        ("biofilm_viability", "Biofilm Viability"),
        ("motility_assay_stats", "Motility Assay Stats"),
        ("chemotaxis_index", "Chemotaxis Index"),
        ("spore_count", "Spore Count"),
        ("viability_staining", "Viability Staining"),
        ("flow_cytometry_analysis", "Flow Cytometry Analysis"),
        ("qpcr_relative_quant", "qPCR Relative Quant"),
        ("qpcr_absolute_quant", "qPCR Absolute Quant"),
        ("enzyme_activity", "Enzyme Activity"),
        ("protein_conc", "Protein Concentration"),
        ("dna_rna_conc", "DNA RNA Concentration"),
        ("disk_diffusion_interpret", "Disk Diffusion Interpret"),
        ("mdr_index", "MDR Index"),
        ("efflux_pump_stats", "Efflux Pump Stats"),
        ("beta_lactamase_stats", "Beta Lactamase Stats"),
        ("hemagglutination_titer", "Hemagglutination Titer"),
        ("elisa_quant", "ELISA Quant"),
        ("western_blot_quant", "Western Blot Quant"),
        ("gel_band_quant", "Gel Band Quant"),
        ("susceptibility_profile_stats", "Susceptibility Profile Stats")
    ], validators=[DataRequired()])
    submit = SubmitField("Analyze")

class AnalysisTypeForm(FlaskForm):
    manual_analysis_type = SelectField("Analysis Type", choices=[
        ("growth_stats", "Growth Rate & Doubling Time"),
        ("cfu_stats", "CFU Count Statistics"),
        ("zone", "Zone of Inhibition"),
        ("diversity", "Microbial Diversity Indices"),
        ("reed_muench", "Reed-Muench TCID50"),
        ("mic_manual", "MIC Manual"),
        ("d_value_manual", "D-Value Manual"),
        ("z_value_manual", "Z-Value Manual"),
        ("lag_phase_manual", "Lag Phase Manual"),
        ("mu_max_manual", "Mu Max Manual"),
        ("mbc_manual", "MBC Manual"),
        ("f_value_manual", "F-Value Manual"),
        ("synergy_testing_manual", "Synergy Testing Manual"),
        ("mutation_frequency_manual", "Mutation Frequency Manual"),
        ("pfu_manual", "PFU Manual"),
        ("spore_count_manual", "Spore Count Manual"),
        ("viability_staining_manual", "Viability Staining Manual"),
        ("enzyme_activity_manual", "Enzyme Activity Manual"),
        ("protein_conc_manual", "Protein Concentration Manual"),
        ("dna_rna_conc_manual", "DNA RNA Conc Manual"),
        ("hemagglutination_titer_manual", "Hemagglutination Titer Manual"),
        ("elisa_quant_manual", "ELISA Quant Manual")

    ], validators=[DataRequired()])
    submit = SubmitField("Select Analysis")

class GrowthRateForm(FlaskForm):
    initial_population = FloatField('Initial Population (N₀)', validators=[DataRequired()])
    final_population = FloatField('Final Population (N)', validators=[DataRequired()])
    time_start = FloatField('Time_Start (t1)', validators=[DataRequired()])
    time_end = FloatField('Time (t2)', validators=[DataRequired()])
    submit = SubmitField("Calculate Growth Rate")

class CFUForm(FlaskForm):
    colony_count = IntegerField("Colony Count", validators=[DataRequired()])
    dilution_factor = FloatField("Dilution Factor", validators=[DataRequired()])
    volume_plated = FloatField("Volume Plated (mL)", validators=[DataRequired()])
    submit = SubmitField("Calculate CFU/mL")

class ZoneOfInhibitionForm(FlaskForm):
    diameter = FloatField("Zone Diameter (mm)", validators=[DataRequired()])
    submit = SubmitField("Submit Zone of Inhibition")

class DiversityIndexForm(FlaskForm):
    species_counts = TextAreaField("Species Counts (comma-separated)", validators=[DataRequired()])
    method = SelectField("Method", choices=[("shannon", "Shannon Index"), ("simpson", "Simpson Index")])
    submit = SubmitField("Calculate Diversity Index")

class ReedMuenchForm(FlaskForm):
    dilutions = TextAreaField("Dilution Levels (comma-separated)", validators=[DataRequired()])
    infected = TextAreaField("Number Infected (comma-separated)", validators=[DataRequired()])
    total = TextAreaField("Total Tested (comma-separated)", validators=[DataRequired()])
    submit = SubmitField("Calculate TCID₅₀")

class MICManualForm(FlaskForm):
    concentrations = TextAreaField("Concentrations (comma-separated)", validators=[DataRequired()])
    growths = TextAreaField("Growth (0/1 comma-separated)", validators=[DataRequired()])
    submit = SubmitField("Calculate MIC")

class DValueManualForm(FlaskForm):
    times = TextAreaField("Times (comma-separated)", validators=[DataRequired()])
    log_cfus = TextAreaField("Log CFU (comma-separated)", validators=[DataRequired()])
    submit = SubmitField("Calculate D-Value")

class ZValueManualForm(FlaskForm):
    temps = TextAreaField("Temperatures (comma-separated)", validators=[DataRequired()])
    d_values = TextAreaField("D-Values (comma-separated)", validators=[DataRequired()])
    submit = SubmitField("Calculate Z-Value")

class LagPhaseManualForm(FlaskForm):
    times = TextAreaField("Times (comma-separated)", validators=[DataRequired()])
    ods = TextAreaField("ODs (comma-separated)", validators=[DataRequired()])
    submit = SubmitField("Estimate Lag Phase")

class MuMaxManualForm(FlaskForm):
    times = TextAreaField("Times (comma-separated)", validators=[DataRequired()])
    ods = TextAreaField("ODs (comma-separated)", validators=[DataRequired()])
    submit = SubmitField("Calculate Mu Max")

class MBCManualForm(FlaskForm):
    concentrations = TextAreaField("Concentrations (comma-separated)", validators=[DataRequired()])
    viabilities = TextAreaField("Viabilities (comma-separated)", validators=[DataRequired()])
    submit = SubmitField("Calculate MBC")

class FValueManualForm(FlaskForm):
    times = TextAreaField("Times (comma-separated)", validators=[DataRequired()])
    temps = TextAreaField("Temperatures (comma-separated)", validators=[DataRequired()])
    ref_temp = FloatField("Reference Temperature", default=121.1)
    z_value = FloatField("Z Value", default=10)
    submit = SubmitField("Calculate F-Value")

class SynergyTestingManualForm(FlaskForm):
    drug_a_mic = FloatField("Drug A MIC", validators=[DataRequired()])
    drug_b_mic = FloatField("Drug B MIC", validators=[DataRequired()])
    combo_mic_a = FloatField("Combo MIC A", validators=[DataRequired()])
    combo_mic_b = FloatField("Combo MIC B", validators=[DataRequired()])
    submit = SubmitField("Calculate FICI")

class MutationFrequencyManualForm(FlaskForm):
    mutants = IntegerField("Mutants", validators=[DataRequired()])
    total = IntegerField("Total", validators=[DataRequired()])
    submit = SubmitField("Calculate Mutation Frequency")

class PFUManualForm(FlaskForm):
    plaques = IntegerField("Plaques", validators=[DataRequired()])
    dilution = FloatField("Dilution", validators=[DataRequired()])
    volume = FloatField("Volume", validators=[DataRequired()])
    submit = SubmitField("Calculate PFU")

class SporeCountManualForm(FlaskForm):
    spores = IntegerField("Spores", validators=[DataRequired()])
    total = IntegerField("Total", validators=[DataRequired()])
    submit = SubmitField("Calculate Spore Percentage")

class ViabilityStainingManualForm(FlaskForm):
    live = IntegerField("Live", validators=[DataRequired()])
    dead = IntegerField("Dead", validators=[DataRequired()])
    submit = SubmitField("Calculate Live/Dead Ratio")

class EnzymeActivityManualForm(FlaskForm):
    times = TextAreaField("Times (comma-separated)", validators=[DataRequired()])
    products = TextAreaField("Products (comma-separated)", validators=[DataRequired()])
    submit = SubmitField("Calculate Enzyme Activity")

class ProteinConcManualForm(FlaskForm):
    absorbance = FloatField("Absorbance", validators=[DataRequired()])
    factor = FloatField("Factor", default=1.0)
    submit = SubmitField("Calculate Protein Concentration")

class DNARNAConcManualForm(FlaskForm):
    a260 = FloatField("A260", validators=[DataRequired()])
    factor = FloatField("Factor", default=50)
    submit = SubmitField("Calculate DNA/RNA Concentration")

class HemagglutinationTiterManualForm(FlaskForm):
    dilutions = TextAreaField("Dilutions (comma-separated)", validators=[DataRequired()])
    agglutinations = TextAreaField("Agglutinations (0/1 comma-separated)", validators=[DataRequired()])
    submit = SubmitField("Calculate Titer")

class ELISAQuantManualForm(FlaskForm):
    absorbances = TextAreaField("Absorbances (comma-separated)", validators=[DataRequired()])
    standard_concs = TextAreaField("Standard Concs (comma-separated)", validators=[DataRequired()])
    submit = SubmitField("Calculate Quantification")





