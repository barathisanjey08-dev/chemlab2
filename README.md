# ChemLab Virtual Laboratory

A responsive chemoinformatics website built with Flask, RDKit, and 3Dmol.js.

**Developer:** Sanjey Barathi M I  
**Program:** Artificial Intelligence and Machine Learning (AI & ML)

## Features
- SMILES input and example molecule library
- RDKit-powered 2D structure and molecular descriptors
- Interactive 3D molecular viewer
- Molecule comparison
- CSV export and PDF report
- Theory and knowledge-check sections

## Run locally
Use Python 3.11:
```bash
python -m venv .venv
# Windows: .venv\Scripts\activate
pip install -r requirements.txt
python app.py
```
Open http://127.0.0.1:5000

## Render deployment
Build command: `pip install -r requirements.txt`  
Start command: `gunicorn app:app`

The 3D viewer loads 3Dmol.js from its public CDN and requires internet access.
