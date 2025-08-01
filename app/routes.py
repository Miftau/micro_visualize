from flask import Blueprint, render_template, request, redirect, flash, url_for
from werkzeug.utils import secure_filename
import os
import pandas as pd
import matplotlib.pyplot as plt
from io import BytesIO
import base64

from .forms import UploadForm
from .utils import analysis as ana
from .utils import predict

main = Blueprint("main", __name__)

@main.route("/", methods=["GET", "POST"])
def index():
    form = UploadForm()
    plot_url, table = None, None

    if form.validate_on_submit():
        file = form.file.data
        analysis_type = form.analysis_type.data
        filename = secure_filename(file.filename)
        file_path = os.path.join(os.getcwd(), 'microhub_flask', 'uploads', filename)
        file.save(file_path)

        # Load DataFrame
        ext = filename.rsplit('.', 1)[1].lower()
        df = pd.read_csv(file_path) if ext == 'csv' else pd.read_excel(file_path)

        # Run analysis or prediction
        if analysis_type == "growth_curve":
            fig = ana.plot_growth_curve(df)
        elif analysis_type == "microbiome":
            fig = ana.plot_microbiome_composition(df)
        elif analysis_type == "amr":
            fig = ana.plot_amr_heatmap(df)
        elif analysis_type == "gene_expression":
            fig = ana.plot_gene_expression(df)
        elif analysis_type == "predict_resistance":
            result_df = predict.predict_resistance(df)
            table = result_df.to_html(classes="table table-bordered", index=False)
            return render_template("index.html", form=form, table=table)

        # Convert figure to base64 for inline HTML
        img = BytesIO()
        fig.savefig(img, format='png', bbox_inches='tight')
        img.seek(0)
        plot_url = 'data:image/png;base64,' + base64.b64encode(img.read()).decode()

    return render_template("index.html", form=form, plot_url=plot_url, table=table)
