from flask import Blueprint, render_template, request, flash
from werkzeug.utils import secure_filename
import os
import pandas as pd
from .forms import *
import json
from server.utils.analysis import *
from server.utils.predict import predict_resistance

bp = Blueprint('main', __name__)
UPLOAD_FOLDER = "uploads/"
# Load formulas dataset
def load_formulas():
    file_path = os.path.join(os.path.dirname(__file__), "formulas.json")
    with open(file_path, "r", encoding="utf-8") as f:
        return json.load(f)



@bp.route("/", methods=["GET", "POST"])
def index():
    form = UploadForm()
    result = plot_url = table = stats = ttest = None
    if form.validate_on_submit():
        if request.method == "POST":
            file = request.files["file"]
            analysis_type = request.form.get("analysis_type")
            if file and file.filename.endswith((".csv", ".xlsx")):
                filename = secure_filename(file.filename)
                filepath = os.path.join(UPLOAD_FOLDER, filename)
                file.save(filepath)
                try:
                    if filename.endswith(".xlsx"):
                        df = pd.read_excel(filepath)
                    else:
                        df = pd.read_csv(filepath)

                    if analysis_type == "growth_curve":
                        result, plot_url = perform_analysis(df, "growth_curve")
                    elif analysis_type == "amr_heatmap":
                        result, plot_url = perform_analysis(df, "amr_heatmap")
                    elif analysis_type == "amr_prediction":
                        table = predict_resistance(df)
                    elif analysis_type == "growth_stats":
                        stats = growth_stats(df)
                    elif analysis_type == "inhibition_stats":
                        stats, ttest = inhibition_stats(df)
                    elif analysis_type == "cfu_stats":
                        table = cfu_stats(df).to_html(classes="table table-bordered")
                    elif analysis_type == "diversity_stats":
                        table = compute_diversity(df).to_html(classes="table table-bordered")
                except Exception as e:
                    flash(f"Error processing file: {str(e)}", "danger")
            else:
                flash("Invalid file type. Please upload a CSV or Excel file.", "danger")
    else:
        for field, errors in form.errors.items():
            for error in errors:
                flash(f"{field}: {error}", "danger")

    return render_template("index.html", form=form, result=result, plot_url=plot_url, table=table, stats=stats,
                           ttest=ttest)


@bp.route("/manual-input", methods=["GET", "POST"])
def manual_input():
    from .forms import AnalysisTypeForm, GrowthRateForm, CFUForm, ZoneOfInhibitionForm, DiversityIndexForm, ReedMuenchForm

    manual_analysis_type_form = AnalysisTypeForm()
    growth_form = GrowthRateForm()
    cfu_form = CFUForm()
    zone_form = ZoneOfInhibitionForm()
    diversity_form = DiversityIndexForm()
    reed_muench_form = ReedMuenchForm()

    selected_analysis = request.args.get("analysis_type") or request.form.get("analysis_type")

    stats = table = error = None

    if request.method == "POST":
        if manual_analysis_type_form.validate_on_submit():
            selected_analysis = manual_analysis_type_form.manual_analysis_type.data
        try:
            if selected_analysis == "growth_stats" and growth_form.validate_on_submit():
                initial_pop = growth_form.initial_population.data
                final_pop = growth_form.final_population.data
                time_start = growth_form.time_start.data
                time_end = growth_form.time_end.data

                if time_end != time_start and initial_pop > 0 and final_pop > 0:
                    growth_rate = (np.log(final_pop / initial_pop)) / (time_end - time_start)
                    stats = {
                        "Initial Population": initial_pop,
                        "Final Population": final_pop,
                        "Time Interval (h)": time_end - time_start,
                        "Growth Rate (h⁻¹)": round(growth_rate, 4),
                        "Doubling Time (h)": round(np.log(2) / growth_rate, 4) if growth_rate > 0 else "N/A"
                    }
                else:
                    error = "Invalid input: Ensure time interval is non-zero and populations are positive."
                    flash(error, "danger")

            elif selected_analysis == "cfu_stats" and cfu_form.validate_on_submit():
                colony_count = cfu_form.colony_count.data
                dilution_factor = cfu_form.dilution_factor.data
                volume_plated = cfu_form.volume_plated.data

                if volume_plated > 0 and dilution_factor > 0:
                    cfu_ml = (colony_count * dilution_factor) / volume_plated
                    stats = {
                        "Colony Count": colony_count,
                        "Dilution Factor": dilution_factor,
                        "Volume Plated (mL)": volume_plated,
                        "CFU/mL": round(cfu_ml, 2)
                    }
                else:
                    error = "Invalid input: Dilution factor and volume plated must be positive."
                    flash(error, "danger")

            elif selected_analysis == "zone" and zone_form.validate_on_submit():
                diameter = zone_form.diameter.data
                stats = {"Zone Diameter (mm)": round(diameter, 2)}

            elif selected_analysis == "diversity" and diversity_form.validate_on_submit():
                counts = [int(x) for x in diversity_form.species_counts.data.split(",") if x.strip()]
                method = diversity_form.method.data
                if counts and all(c >= 0 for c in counts):
                    df = pd.DataFrame({"species": [f"S{i + 1}" for i in range(len(counts))], "count": counts})
                    table = compute_diversity(df, method=method).to_html(classes="table table-bordered")
                else:
                    error = "Invalid input: Species counts must be non-negative integers."
                    flash(error, "danger")

            elif selected_analysis == "reed_muench" and reed_muench_form.validate_on_submit():
                try:
                    dilutions = [float(x) for x in reed_muench_form.dilutions.data.split(",")]
                    infected = [int(x) for x in reed_muench_form.infected.data.split(",")]
                    total = [int(x) for x in reed_muench_form.total.data.split(",")]

                    if len(dilutions) == len(infected) == len(total):
                        tcid50 = calculate_reed_muench(dilutions, infected, total)
                        stats = {"Estimated TCID₅₀": tcid50}
                    else:
                        error = "All input lists must have the same length."
                        flash(error, "danger")
                except Exception as e:
                    flash(f"Error: {str(e)}", "danger")

            else:
                for form in [manual_analysis_type_form, growth_form, cfu_form, zone_form, diversity_form]:
                    for field, errors in form.errors.items():
                        for error in errors:
                            flash(f"{form[field].label.text}: {error}", "danger")

        except Exception as e:
            error = f"Error processing input: {str(e)}"
            flash(error, "danger")

    return render_template(
        "manual_input.html",
        analysis_type_form=manual_analysis_type_form,
        growth_form=growth_form,
        cfu_form=cfu_form,
        zone_form=zone_form,
        diversity_form=diversity_form,
        selected_analysis=selected_analysis,
        reed_muench_form=reed_muench_form,
        stats=stats,
        table=table,
        error=error
    )


@bp.route("/formulas")
def show_formulas():
    formulas = load_formulas()

    # Group formulas by category
    categories = {}
    for formula in formulas:
        category = formula.get("category", "Other")
        if category not in categories:
            categories[category] = []
        categories[category].append(formula)

    return render_template("formulas.html", categories=categories)



