# ChemLab Studio — Molecular Discovery Workspace

**Developer:** Sanjey Barathi M I  
**Programme:** Artificial Intelligence and Machine Learning (AI & ML)

## Features
- Curated molecule library with plain-language explanations, uses, and cautions
- SMILES parsing and canonicalization with RDKit
- Browser-rendered 2D structure and interactive generated 3D conformer viewer
- Molecular descriptors, charts, and two-molecule comparisons
- Professional A4 PDF practical reports with student name, registration number, course, date, aims, theory, procedure, molecular input, computed results, interpretation, conclusion, and page footer
- OpenAI-powered chemistry assistant using the Responses API
- Theory and mini quiz

## Deploy on Streamlit Community Cloud
- Repository branch: `main`
- Main file path: `app.py`
- Dependencies file: `requirements.txt`

## Configure OpenAI
1. Create an API key at https://platform.openai.com/api-keys.
2. In Streamlit Community Cloud, open **Manage app → Settings → Secrets**.
3. Add this TOML line, using your own key:

```toml
OPENAI_API_KEY = "your-key-here"
```

4. Save and wait for the app to restart. Never put API keys in `app.py`, commit them to GitHub, or share them in screenshots. ChatGPT subscriptions and API billing are separate; API use may require billing and is usage-based. The app currently requests `gpt-4.1-mini`; if that model is not enabled for your API project, change the model identifier in `openai_answer()` to a model available to your account.

## Notes
The 2D and 3D viewers load JavaScript libraries from public CDNs and need internet access. RDKit descriptor values and generated conformers are computational results, not experimental measurements. They do not by themselves establish safety, medical suitability, or biological activity.
