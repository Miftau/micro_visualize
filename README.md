
```markdown
# 🧪 MicroAnalyze – Microbiology Data Analysis App

Microhub is a **Flask-based web application** for microbiologists and biotech professionals to perform:
- 📈 Growth rate calculations
- 🦠 Zone of inhibition analysis
- 📊 CFU and diversity indices
- 🧬 Virus quantification (Reed & Muench)
- 📂 CSV upload **and** manual data entry
- 📥 Export results and charts

---

## 🚀 Features
- Toggle between CSV upload and manual input
- Instant calculations with downloadable results
- Support for multiple microbiology-specific methods
- Ready for deployment on **Render** or other PaaS

---

## 🛠 Installation
```bash
git clone https://github.com/username/microhub.git
cd microhub
python -m venv venv
source venv/bin/activate   # Linux/macOS
venv\Scripts\activate      # Windows
pip install -r requirements.txt
flask run
