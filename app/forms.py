# microhub_flask/app/forms.py

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
        ("amr_prediction", "AMR Prediction (ML)"),
        ("growth_stats", "Growth Rate & Doubling Time"),
        ("inhibition_stats", "Zone of Inhibition Stats"),
        ("cfu_stats", "CFU Count Statistics"),
        ("diversity_stats", "Microbial Diversity Indices")
    ], validators=[DataRequired()])
    submit = SubmitField("Analyze")

class AnalysisTypeForm(FlaskForm):
    manual_analysis_type = SelectField("Analysis Type", choices=[
        ("growth_stats", "Growth Rate & Doubling Time"),
        ("cfu_stats", "CFU Count Statistics"),
        ("zone", "Zone of Inhibition"),
        ("diversity", "Microbial Diversity Indices")
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



