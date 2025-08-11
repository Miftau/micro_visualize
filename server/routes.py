# routes.py (Fully Extended)

from flask import Blueprint, render_template, request, flash
from werkzeug.utils import secure_filename
import os
import pandas as pd
from .forms import *
import json
from server.utils.analysis import *
from server.utils.predict import predict_resistance  # Assume exists

bp = Blueprint('main', __name__)
UPLOAD_FOLDER = "uploads/"

def load_formulas():
    file_path = os.path.join(os.path.dirname(__file__), "formulas.json")
    with open(file_path, "r", encoding="utf-8") as f:
        return json.load(f)

@bp.route("/", methods=["GET", "POST"])
def index():
    form = UploadForm()
    result = plot_url = table = stats = ttest = None
    if form.validate_on_submit():
        file = request.files["file"]
        analysis_type = form.analysis_type.data
        if file and file.filename.endswith((".csv", ".xlsx")):
            filename = secure_filename(file.filename)
            filepath = os.path.join(UPLOAD_FOLDER, filename)
            file.save(filepath)
            try:
                if filename.endswith(".xlsx"):
                    df = pd.read_excel(filepath)
                else:
                    df = pd.read_csv(filepath)

                # Visual analyses
                if analysis_type in [
                    "growth_curve", "amr_heatmap", "time_kill_curve", "survivorship_curve", "ph_change_curve",
                    "pca_microbiome", "rarefaction_curve", "volcano_plot", "gene_expression_heatmap",
                    "cooccurrence_network", "biofilm_biomass_plot", "motility_plot", "fermentation_gas_plot",
                    "oxygen_consumption_plot", "quorum_sensing_plot", "enzyme_kinetics_plot", "phylogenetic_tree",
                    "dose_response_curve", "synergy_heatmap", "flow_cytometry_histogram", "qpcr_ct_plot",
                    "western_blot_densitometry", "elisa_standard_curve", "hemagglutination_titer_plot",
                    "plaque_assay_plot", "beta_lactamase_activity_plot", "efflux_pump_activity",
                    "mdr_index_heatmap", "susceptibility_profile", "gel_band_quant_plot"
                ]:
                    result, plot_url = perform_analysis(df, analysis_type)
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
                elif analysis_type == "mic_determination":
                    stats = mic_determination(df)
                elif analysis_type == "mbc_determination":
                    stats = mbc_determination(df)
                elif analysis_type == "d_value_calc":
                    stats = d_value_calc(df)
                elif analysis_type == "z_value_calc":
                    stats = z_value_calc(df)
                elif analysis_type == "f_value_calc":
                    stats = f_value_calc(df)
                elif analysis_type == "lag_phase_est":
                    stats = lag_phase_est(df)
                elif analysis_type == "mu_max_calc":
                    stats = mu_max_calc(df)
                elif analysis_type == "baranyi_model_fit":
                    stats = baranyi_model_fit(df)
                elif analysis_type == "gompertz_model_fit":
                    stats = gompertz_model_fit(df)
                elif analysis_type == "logistic_model_fit":
                    stats = logistic_model_fit(df)
                elif analysis_type == "alpha_diversity":
                    table = alpha_diversity(df).to_html(classes="table table-bordered")
                elif analysis_type == "beta_diversity":
                    table = beta_diversity(df).to_html(classes="table table-bordered")
                elif analysis_type == "phylogenetic_diversity":
                    stats = phylogenetic_diversity(df)
                elif analysis_type == "mutation_frequency":
                    stats = mutation_frequency(df)
                elif analysis_type == "pfu_calc":
                    table = pfu_calc(df).to_html(classes="table table-bordered")
                elif analysis_type == "synergy_testing":
                    stats = synergy_testing(df)
                elif analysis_type == "post_antibiotic_effect":
                    stats = post_antibiotic_effect(df)
                elif analysis_type == "biofilm_biomass_stats":
                    table = biofilm_biomass_stats(df).to_html(classes="table table-bordered")
                elif analysis_type == "biofilm_viability":
                    stats = biofilm_viability(df)
                elif analysis_type == "motility_assay_stats":
                    stats = motility_assay_stats(df)
                elif analysis_type == "chemotaxis_index":
                    stats = chemotaxis_index(df)
                elif analysis_type == "spore_count":
                    stats = spore_count(df)
                elif analysis_type == "viability_staining":
                    stats = viability_staining(df)
                elif analysis_type == "flow_cytometry_analysis":
                    table = flow_cytometry_analysis(df).to_html(classes="table table-bordered")
                elif analysis_type == "qpcr_relative_quant":
                    stats = qpcr_relative_quant(df)
                elif analysis_type == "qpcr_absolute_quant":
                    stats = qpcr_absolute_quant(df)
                elif analysis_type == "enzyme_activity":
                    stats = enzyme_activity(df)
                elif analysis_type == "protein_conc":
                    stats = protein_conc(df)
                elif analysis_type == "dna_rna_conc":
                    stats = dna_rna_conc(df)
                elif analysis_type == "disk_diffusion_interpret":
                    table = disk_diffusion_interpret(df).to_html(classes="table table-bordered")
                elif analysis_type == "mdr_index":
                    stats = mdr_index(df)
                elif analysis_type == "efflux_pump_stats":
                    stats = efflux_pump_stats(df)
                elif analysis_type == "beta_lactamase_stats":
                    stats = beta_lactamase_stats(df)
                elif analysis_type == "hemagglutination_titer":
                    stats = hemagglutination_titer(df)
                elif analysis_type == "elisa_quant":
                    stats = elisa_quant(df)
                elif analysis_type == "western_blot_quant":
                    stats = western_blot_quant(df)
                elif analysis_type == "gel_band_quant":
                    stats = gel_band_quant(df)
                elif analysis_type == "susceptibility_profile_stats":
                    table = susceptibility_profile_stats(df).to_html(classes="table table-bordered")
            except Exception as e:
                flash(f"Error processing file: {str(e)}", "danger")
        else:
            flash("Invalid file type. Please upload a CSV or Excel file.", "danger")
    else:
        for field, errors in form.errors.items():
            for error in errors:
                flash(f"{field}: {error}", "danger")

    return render_template("index.html", form=form, result=result, plot_url=plot_url, table=table, stats=stats, ttest=ttest)

@bp.route("/manual-input", methods=["GET", "POST"])
def manual_input():
    analysis_type_form = AnalysisTypeForm()
    growth_form = GrowthRateForm()
    cfu_form = CFUForm()
    zone_form = ZoneOfInhibitionForm()
    diversity_form = DiversityIndexForm()
    reed_muench_form = ReedMuenchForm()
    mic_manual_form = MICManualForm()
    d_value_manual_form = DValueManualForm()
    z_value_manual_form = ZValueManualForm()
    lag_phase_manual_form = LagPhaseManualForm()
    mu_max_manual_form = MuMaxManualForm()
    mbc_manual_form = MBCManualForm()
    f_value_manual_form = FValueManualForm()
    synergy_testing_manual_form = SynergyTestingManualForm()
    mutation_frequency_manual_form = MutationFrequencyManualForm()
    pfu_manual_form = PFUManualForm()
    spore_count_manual_form = SporeCountManualForm()
    viability_staining_manual_form = ViabilityStainingManualForm()
    enzyme_activity_manual_form = EnzymeActivityManualForm()
    protein_conc_manual_form = ProteinConcManualForm()
    dna_rna_conc_manual_form = DNARNAConcManualForm()
    hemagglutination_titer_manual_form = HemagglutinationTiterManualForm()
    elisa_quant_manual_form = ELISAQuantManualForm()

    selected_analysis = request.args.get("analysis_type") or request.form.get("manual_analysis_type")

    stats = table = error = None

    if request.method == "POST":
        if analysis_type_form.validate_on_submit():
            selected_analysis = analysis_type_form.manual_analysis_type.data
        try:
            if selected_analysis == "growth_stats" and growth_form.validate_on_submit():
                initial_pop = growth_form.initial_population.data
                final_pop = growth_form.final_population.data
                time_start = growth_form.time_start.data
                time_end = growth_form.time_end.data
                time = time_end - time_start
                if time > 0 and initial_pop > 0 and final_pop > 0:
                    result = perform_growth_rate(initial_pop, final_pop, time)
                    stats = {"Result": result}
                else:
                    error = "Invalid input."
                    flash(error, "danger")

            elif selected_analysis == "cfu_stats" and cfu_form.validate_on_submit():
                colony_count = cfu_form.colony_count.data
                dilution_factor = cfu_form.dilution_factor.data
                volume_plated = cfu_form.volume_plated.data
                result = perform_cfu(colony_count, dilution_factor, volume_plated)
                stats = {"Result": result}

            elif selected_analysis == "zone" and zone_form.validate_on_submit():
                diameter = zone_form.diameter.data
                result = perform_zone_of_inhibition(diameter)
                stats = {"Result": result}

            elif selected_analysis == "diversity" and diversity_form.validate_on_submit():
                counts = [int(x.strip()) for x in diversity_form.species_counts.data.split(",") if x.strip()]
                method = diversity_form.method.data
                result = perform_diversity_index(counts, method)
                stats = {"Result": result}

            elif selected_analysis == "reed_muench" and reed_muench_form.validate_on_submit():
                dilutions = [float(x) for x in reed_muench_form.dilutions.data.split(",")]
                infected = [int(x) for x in reed_muench_form.infected.data.split(",")]
                total = [int(x) for x in reed_muench_form.total.data.split(",")]
                tcid50 = calculate_reed_muench(dilutions, infected, total)
                stats = {"Estimated TCID₅₀": tcid50}

            elif selected_analysis == "mic_manual" and mic_manual_form.validate_on_submit():
                concentrations = [float(x) for x in mic_manual_form.concentrations.data.split(",")]
                growths = [float(x) for x in mic_manual_form.growths.data.split(",")]
                result = perform_mic_manual(concentrations, growths)
                stats = {"Result": result}

            elif selected_analysis == "d_value_manual" and d_value_manual_form.validate_on_submit():
                times = [float(x) for x in d_value_manual_form.times.data.split(",")]
                log_cfus = [float(x) for x in d_value_manual_form.log_cfus.data.split(",")]
                result = perform_d_value_manual(times, log_cfus)
                stats = {"Result": result}

            elif selected_analysis == "z_value_manual" and z_value_manual_form.validate_on_submit():
                temps = [float(x) for x in z_value_manual_form.temps.data.split(",")]
                d_values = [float(x) for x in z_value_manual_form.d_values.data.split(",")]
                result = perform_z_value_manual(temps, d_values)
                stats = {"Result": result}

            elif selected_analysis == "lag_phase_manual" and lag_phase_manual_form.validate_on_submit():
                times = [float(x) for x in lag_phase_manual_form.times.data.split(",")]
                ods = [float(x) for x in lag_phase_manual_form.ods.data.split(",")]
                df = pd.DataFrame({"Time": times, "OD": ods})
                stats = lag_phase_est(df)

            elif selected_analysis == "mu_max_manual" and mu_max_manual_form.validate_on_submit():
                times = [float(x) for x in mu_max_manual_form.times.data.split(",")]
                ods = [float(x) for x in mu_max_manual_form.ods.data.split(",")]
                result = perform_mu_max_manual(times, ods)
                stats = {"Result": result}

            elif selected_analysis == "mbc_manual" and mbc_manual_form.validate_on_submit():
                concentrations = [float(x) for x in mbc_manual_form.concentrations.data.split(",")]
                viabilities = [float(x) for x in mbc_manual_form.viabilities.data.split(",")]
                result = perform_mbc_manual(concentrations, viabilities)
                stats = {"Result": result}

            elif selected_analysis == "f_value_manual" and f_value_manual_form.validate_on_submit():
                times = [float(x) for x in f_value_manual_form.times.data.split(",")]
                temps = [float(x) for x in f_value_manual_form.temps.data.split(",")]
                ref_temp = f_value_manual_form.ref_temp.data
                z_value = f_value_manual_form.z_value.data
                result = perform_f_value_manual(times, temps, ref_temp, z_value)
                stats = {"Result": result}

            elif selected_analysis == "synergy_testing_manual" and synergy_testing_manual_form.validate_on_submit():
                drug_a_mic = synergy_testing_manual_form.drug_a_mic.data
                drug_b_mic = synergy_testing_manual_form.drug_b_mic.data
                combo_mic_a = synergy_testing_manual_form.combo_mic_a.data
                combo_mic_b = synergy_testing_manual_form.combo_mic_b.data
                df = pd.DataFrame({"DrugA_MIC": [drug_a_mic], "DrugB_MIC": [drug_b_mic], "Combo_MIC_A": [combo_mic_a], "Combo_MIC_B": [combo_mic_b]})
                stats = synergy_testing(df)

            elif selected_analysis == "mutation_frequency_manual" and mutation_frequency_manual_form.validate_on_submit():
                mutants = mutation_frequency_manual_form.mutants.data
                total = mutation_frequency_manual_form.total.data
                freq = mutants / total if total > 0 else 0
                stats = {"Mutation Frequency": round(freq, 6)}

            elif selected_analysis == "pfu_manual" and pfu_manual_form.validate_on_submit():
                plaques = pfu_manual_form.plaques.data
                dilution = pfu_manual_form.dilution.data
                volume = pfu_manual_form.volume.data
                pfu = (plaques * dilution) / volume if volume > 0 else 0
                stats = {"PFU": round(pfu, 2)}

            elif selected_analysis == "spore_count_manual" and spore_count_manual_form.validate_on_submit():
                spores = spore_count_manual_form.spores.data
                total = spore_count_manual_form.total.data
                percent = (spores / total * 100) if total > 0 else 0
                stats = {"Spore Percentage": round(percent, 2)}

            elif selected_analysis == "viability_staining_manual" and viability_staining_manual_form.validate_on_submit():
                live = viability_staining_manual_form.live.data
                dead = viability_staining_manual_form.dead.data
                ratio = live / dead if dead > 0 else "Infinite"
                stats = {"Live/Dead Ratio": ratio}

            elif selected_analysis == "enzyme_activity_manual" and enzyme_activity_manual_form.validate_on_submit():
                times = [float(x) for x in enzyme_activity_manual_form.times.data.split(",")]
                products = [float(x) for x in enzyme_activity_manual_form.products.data.split(",")]
                df = pd.DataFrame({"Time": times, "Product": products})
                stats = enzyme_activity(df)

            elif selected_analysis == "protein_conc_manual" and protein_conc_manual_form.validate_on_submit():
                absorbance = protein_conc_manual_form.absorbance.data
                factor = protein_conc_manual_form.factor.data
                conc = absorbance * factor
                stats = {"Protein Concentration": round(conc, 2)}

            elif selected_analysis == "dna_rna_conc_manual" and dna_rna_conc_manual_form.validate_on_submit():
                a260 = dna_rna_conc_manual_form.a260.data
                factor = dna_rna_conc_manual_form.factor.data
                conc = a260 * factor
                stats = {"DNA/RNA Concentration": round(conc, 2)}

            elif selected_analysis == "hemagglutination_titer_manual" and hemagglutination_titer_manual_form.validate_on_submit():
                dilutions = [float(x) for x in hemagglutination_titer_manual_form.dilutions.data.split(",")]
                agglutinations = [int(x) for x in hemagglutination_titer_manual_form.agglutinations.data.split(",")]
                df = pd.DataFrame({"Dilution": dilutions, "Agglutination": agglutinations})
                stats = hemagglutination_titer(df)

            elif selected_analysis == "elisa_quant_manual" and elisa_quant_manual_form.validate_on_submit():
                absorbances = [float(x) for x in elisa_quant_manual_form.absorbances.data.split(",")]
                standard_concs = [float(x) for x in elisa_quant_manual_form.standard_concs.data.split(",")]
                df = pd.DataFrame({"Absorbance": absorbances, "Standard_Conc": standard_concs})
                stats = elisa_quant(df)

            # Add similar handling for more manual forms

            else:
                forms_list = [analysis_type_form, growth_form, cfu_form, zone_form, diversity_form, reed_muench_form,
                              mic_manual_form, d_value_manual_form, z_value_manual_form, lag_phase_manual_form,
                              mu_max_manual_form, mbc_manual_form, f_value_manual_form, synergy_testing_manual_form,
                              mutation_frequency_manual_form, pfu_manual_form, spore_count_manual_form,
                              viability_staining_manual_form, enzyme_activity_manual_form, protein_conc_manual_form,
                              dna_rna_conc_manual_form, hemagglutination_titer_manual_form, elisa_quant_manual_form]
                for f in forms_list:
                    for field, errors in f.errors.items():
                        for error in errors:
                            flash(f"{f[field].label.text}: {error}", "danger")

        except Exception as e:
            error = f"Error processing input: {str(e)}"
            flash(error, "danger")

    return render_template(
        "manual_input.html",
        analysis_type_form=analysis_type_form,
        growth_form=growth_form,
        cfu_form=cfu_form,
        zone_form=zone_form,
        diversity_form=diversity_form,
        reed_muench_form=reed_muench_form,
        mic_manual_form=mic_manual_form,
        d_value_manual_form=d_value_manual_form,
        z_value_manual_form=z_value_manual_form,
        lag_phase_manual_form=lag_phase_manual_form,
        mu_max_manual_form=mu_max_manual_form,
        mbc_manual_form=mbc_manual_form,
        f_value_manual_form=f_value_manual_form,
        synergy_testing_manual_form=synergy_testing_manual_form,
        mutation_frequency_manual_form=mutation_frequency_manual_form,
        pfu_manual_form=pfu_manual_form,
        spore_count_manual_form=spore_count_manual_form,
        viability_staining_manual_form=viability_staining_manual_form,
        enzyme_activity_manual_form=enzyme_activity_manual_form,
        protein_conc_manual_form=protein_conc_manual_form,
        dna_rna_conc_manual_form=dna_rna_conc_manual_form,
        hemagglutination_titer_manual_form=hemagglutination_titer_manual_form,
        elisa_quant_manual_form=elisa_quant_manual_form,
        # Add more
        selected_analysis=selected_analysis,
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





