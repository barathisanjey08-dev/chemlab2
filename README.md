# ChemLab Studio — Molecular Discovery Workspace

**Developer:** Sanjey Barathi M I  
**Programme:** Artificial Intelligence and Machine Learning (AI & ML)

## Features
- Curated molecule library with plain-language explanations, uses, and cautions
- SMILES parsing and canonicalization with RDKit
- Browser-rendered 2D structure and interactive generated 3D conformer viewer
- Molecular descriptors, charts, and two-molecule comparisons
- Formal A4 student molecular-analysis PDF report styled after the supplied reference: student name and registration number, SMILES, molecular formula, 2D structure, generated 3D ball-and-stick model, descriptor/value/interpretation table, point-wise interpretation, learning note, and page numbering
- OpenAI-powered chemistry assistant using the Responses API
- Chemistry theory and mini quiz
- Display-only data tables convert mixed descriptor values to text to avoid Streamlit/PyArrow serialization errors

## Deploy on Streamlit Community Cloud
- Repository branch: `main`
- Main file path: `app.py`
- Dependencies file: `requirements.txt`

## Configure OpenAI
1. Create an API key at https://platform.openai.com/api-keys.
2. In Streamlit Community Cloud, open **Manage app → Settings → Secrets**.
3. Add this TOML line with your own key:

```toml
OPENAI_API_KEY = "your-key-here"
```

4. Save and wait for the app to restart. Never put API keys in `app.py`, commit them to GitHub, or share them in screenshots.
5. The app currently requests `gpt-4.1-mini`. If your API project does not have access to this model, change the `model` value in `openai_answer()` to a model enabled for your project.

**Billing note:** ChatGPT subscriptions and OpenAI API billing are separate. API use is usage-based and may require billing/credits to be configured.

## If the AI assistant fails
The app now displays a diagnostic next to the fallback answer. Common causes:
- **Authentication error:** the key is incorrect or revoked; create a new one and update Streamlit Secrets.
- **Quota/billing error:** check API billing, credits and project usage limits.
- **Model access error:** select a model available to your API project.
- **Connection error:** retry and check the OpenAI service status.

## Report notes
- The PDF is generated in memory; no API key is needed to create a report.
- The 3D image is a computationally generated conformer, not necessarily an experimentally observed geometry.
- Descriptor values are computational estimates, not experimental measurements, and do not independently establish safety, medical suitability, or biological activity.
- Confirm your institution's report requirements before submitting.

## Dependencies
All runtime packages are listed in `requirements.txt`, including RDKit, ReportLab, OpenAI, Matplotlib and Pillow.
