import streamlit as st
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import joblib
import io
from sklearn.metrics import classification_report, confusion_matrix, accuracy_score

# User login (simple authentication)
USER_CREDENTIALS = {"admin": "microhub123", "researcher": "lab2025"}

def login():
    st.sidebar.title("🔐 Login")
    username = st.sidebar.text_input("Username")
    password = st.sidebar.text_input("Password", type="password")
    if st.sidebar.button("Login"):
        if USER_CREDENTIALS.get(username) == password:
            st.session_state["authenticated"] = True
            st.session_state["username"] = username
            st.success(f"Welcome, {username}!")
        else:
            st.error("Invalid credentials")

def logout():
    if st.sidebar.button("Logout"):
        st.session_state["authenticated"] = False
        st.session_state["username"] = None

# Authenticate user
if "authenticated" not in st.session_state:
    st.session_state["authenticated"] = False

if not st.session_state["authenticated"]:
    login()
    st.stop()
else:
    logout()

# Page config
st.set_page_config(page_title="Microbiology Analyzer", layout="wide")
st.title("🧫 MicroHub - Extended Microbiology Data Analyzer")

# Sidebar navigation
page = st.sidebar.selectbox("Navigate", [
    "Upload & Explore", "Visualizations", "Predict Resistance", "Evaluate Model"
])

uploaded_file = st.sidebar.file_uploader("Upload CSV or Excel File", type=["csv", "xlsx", "xls"])

if uploaded_file:
    if uploaded_file.name.endswith(".csv"):
        df = pd.read_csv(uploaded_file)
    else:
        df = pd.read_excel(uploaded_file)

    st.sidebar.success("✅ File Uploaded")
    st.session_state["data"] = df

if "data" in st.session_state:
    df = st.session_state["data"]
    
    if page == "Upload & Explore":
        st.subheader("📄 Data Preview")
        st.dataframe(df.head())
        st.markdown(f"**Shape:** {df.shape[0]} rows × {df.shape[1]} columns")

    elif page == "Visualizations":
        st.subheader("📊 Visualize Microbiology Data")
        viz = st.selectbox("Select Visualization", [
            "Growth Curve (OD600)", 
            "Microbiome Composition", 
            "Antibiotic Resistance Heatmap", 
            "Gene Expression Heatmap"
        ])

        plot_buf = io.BytesIO()

        if viz == "Growth Curve (OD600)":
            time_col = st.selectbox("Time Column", df.columns)
            od_col = st.selectbox("OD600 Column", df.columns)
            sample_col = st.selectbox("Sample Column", df.columns)
            fig = plt.figure()
            sns.lineplot(data=df, x=time_col, y=od_col, hue=sample_col)
            plt.title("Growth Curve")
            st.pyplot(fig)
            fig.savefig(plot_buf, format="pdf")

        elif viz == "Microbiome Composition":
            genus_col = st.selectbox("Genus Column", df.columns)
            sample_col = st.selectbox("Sample Column", df.columns)
            ab_col = st.selectbox("Abundance Column", df.columns)
            fig = plt.figure()
            pivot = df.pivot_table(index=sample_col, columns=genus_col, values=ab_col, aggfunc="sum", fill_value=0)
            pivot.plot(kind="bar", stacked=True)
            plt.title("Microbiome Composition")
            st.pyplot(plt.gcf())
            plt.savefig(plot_buf, format="pdf")

        elif viz == "Antibiotic Resistance Heatmap":
            strain_col = st.selectbox("Strain Column", df.columns)
            abx_col = st.selectbox("Antibiotic Column", df.columns)
            mic_col = st.selectbox("MIC Column", df.columns)
            fig = plt.figure()
            pivot = df.pivot(index=strain_col, columns=abx_col, values=mic_col)
            sns.heatmap(pivot, annot=True, cmap="coolwarm")
            plt.title("AMR Heatmap")
            st.pyplot(fig)
            fig.savefig(plot_buf, format="pdf")

        elif viz == "Gene Expression Heatmap":
            gene_col = st.selectbox("Gene Column", df.columns)
            cond_col = st.selectbox("Condition Column", df.columns)
            exp_col = st.selectbox("Expression Column", df.columns)
            fig = plt.figure()
            pivot = df.pivot(index=gene_col, columns=cond_col, values=exp_col)
            sns.heatmap(pivot, cmap="viridis")
            plt.title("Gene Expression")
            st.pyplot(fig)
            fig.savefig(plot_buf, format="pdf")

        st.download_button("⬇️ Download Plot as PDF", plot_buf.getvalue(), file_name="plot.pdf", mime="application/pdf")

    elif page == "Predict Resistance":
        st.subheader("🧠 Predict Antibiotic Resistance")
        try:
            model = joblib.load("rf_resistance.pkl")
            features = st.multiselect("Select Feature Columns", df.columns)
            if st.button("Run Prediction"):
                preds = model.predict(df[features])
                df["Predicted_Resistance"] = preds
                st.dataframe(df[features + ["Predicted_Resistance"]].head())

                st.download_button(
                    "⬇️ Download Predictions as CSV",
                    df.to_csv(index=False).encode("utf-8"),
                    "predictions.csv",
                    "text/csv"
                )

                fig = plt.figure()
                sns.countplot(x="Predicted_Resistance", data=df)
                plt.title("Prediction Distribution")
                st.pyplot(fig)

        except Exception as e:
            st.error(f"Prediction failed: {e}")

    elif page == "Evaluate Model":
        st.subheader("📈 Evaluate Model Performance")
        true_col = st.selectbox("Select True Label Column", df.columns)
        pred_col = st.selectbox("Select Predicted Label Column", df.columns)

        if st.button("Evaluate"):
            report = classification_report(df[true_col], df[pred_col], output_dict=True)
            st.dataframe(pd.DataFrame(report).transpose())

            fig = plt.figure()
            cm = confusion_matrix(df[true_col], df[pred_col])
            sns.heatmap(cm, annot=True, fmt="d", cmap="Blues")
            plt.title("Confusion Matrix")
            st.pyplot(fig)
