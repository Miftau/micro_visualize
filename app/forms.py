from flask_wtf import FlaskForm
from wtforms import FileField, SelectField, SubmitField
from wtforms.validators import DataRequired, ValidationError

class UploadForm(FlaskForm):
    file = FileField("Upload CSV or Excel File", validators=[DataRequired()])
    analysis_type = SelectField("Select Analysis Type", choices=[
        ("growth_curve", "Growth Curve"),
        ("microbiome", "Microbiome Composition"),
        ("amr", "AMR Heatmap"),
        ("gene_expression", "Gene Expression"),
        ("predict_resistance", "Predict Resistance")
    ], validators=[DataRequired()])
    submit = SubmitField("Analyze")

    def validate_file(self, field):
        filename = field.data.filename.lower()
        if not (filename.endswith(".csv") or filename.endswith(".xlsx")):
            raise ValidationError("Only .csv and .xlsx files are allowed.")
