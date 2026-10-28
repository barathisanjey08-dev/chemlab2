from flask import Flask, render_template, request, jsonify, send_file
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors
from rdkit.Chem.Draw import rdMolDraw2D
from io import BytesIO
from reportlab.lib.pagesizes import letter
from reportlab.pdfgen import canvas

app = Flask(__name__)
MOLECULES = {
    "Water": "O", "Ethanol": "CCO", "Caffeine": "Cn1c(=O)c2c(ncn2C)n(C)c1=O",
    "Aspirin": "CC(=O)Oc1ccccc1C(=O)O", "Glucose": "C(C1C(C(C(C(O1)O)O)O)O)O",
    "Benzene": "c1ccccc1", "Acetone": "CC(=O)C",
    "Ibuprofen": "CC(C)Cc1ccc(cc1)[C@@H](C)C(=O)O",
    "Paracetamol": "CC(=O)NC1=CC=C(O)C=C1", "Nicotine": "CN1CCC[C@H]1c2cccnc2"
}

def parse_mol(smiles):
    if not isinstance(smiles, str) or not smiles.strip():
        raise ValueError("Enter a SMILES string first.")
    mol = Chem.MolFromSmiles(smiles.strip())
    if mol is None:
        raise ValueError("Invalid SMILES. Check the notation and try again.")
    return mol

def get_descriptors(mol):
    return {
        "Molecular formula": rdMolDescriptors.CalcMolFormula(mol),
        "Molecular weight (g/mol)": round(Descriptors.MolWt(mol), 3),
        "LogP": round(Descriptors.MolLogP(mol), 3),
        "TPSA (Å²)": round(rdMolDescriptors.CalcTPSA(mol), 3),
        "H-bond donors": int(rdMolDescriptors.CalcNumHBD(mol)),
        "H-bond acceptors": int(rdMolDescriptors.CalcNumHBA(mol)),
        "Rotatable bonds": int(rdMolDescriptors.CalcNumRotatableBonds(mol)),
        "Heavy atoms": int(mol.GetNumHeavyAtoms()),
        "Ring count": int(rdMolDescriptors.CalcNumRings(mol))
    }

def structure_svg(mol):
    copy = Chem.Mol(mol)
    Chem.rdDepictor.Compute2DCoords(copy)
    drawer = rdMolDraw2D.MolDraw2DSVG(520, 340)
    drawer.DrawMolecule(copy)
    drawer.FinishDrawing()
    return drawer.GetDrawingText()

def structure_3d(mol):
    mol3d = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = 42
    if AllChem.EmbedMolecule(mol3d, params) == 0:
        try:
            AllChem.MMFFOptimizeMolecule(mol3d, maxIters=150)
        except Exception:
            pass
    return Chem.MolToMolBlock(mol3d)

@app.get("/")
def index():
    return render_template("index.html", molecules=MOLECULES)

@app.get("/health")
def health():
    return {"status": "ok"}

@app.post("/api/analyze")
def analyze():
    data = request.get_json(silent=True) or {}
    try:
        mol = parse_mol(data.get("smiles", ""))
        return jsonify(ok=True, name=data.get("name") or "Custom molecule",
                       smiles=Chem.MolToSmiles(mol), descriptors=get_descriptors(mol),
                       svg=structure_svg(mol), molblock=structure_3d(mol),
                       atoms=mol.GetNumAtoms(), bonds=mol.GetNumBonds())
    except Exception as exc:
        return jsonify(ok=False, error=str(exc)), 400

@app.post("/api/compare")
def compare():
    data = request.get_json(silent=True) or {}
    try:
        a, b = parse_mol(data.get("smiles_a", "")), parse_mol(data.get("smiles_b", ""))
        return jsonify(ok=True, a={"name": data.get("name_a", "Molecule A"), "descriptors": get_descriptors(a)},
                       b={"name": data.get("name_b", "Molecule B"), "descriptors": get_descriptors(b)})
    except Exception as exc:
        return jsonify(ok=False, error=str(exc)), 400

@app.post("/api/report.pdf")
def report_pdf():
    data = request.get_json(silent=True) or {}
    try:
        mol = parse_mol(data.get("smiles", ""))
        buff = BytesIO()
        pdf = canvas.Canvas(buff, pagesize=letter)
        y = 740
        pdf.setFont("Helvetica-Bold", 18); pdf.drawString(50, y, "ChemLab Molecular Report"); y -= 32
        pdf.setFont("Helvetica", 11)
        for line in [f"Developer: Sanjey Barathi M I", f"Molecule: {data.get('name', 'Custom molecule')}",
                     f"Canonical SMILES: {Chem.MolToSmiles(mol)}", "Calculated descriptors:"]:
            pdf.drawString(50, y, line); y -= 22
        for key, value in get_descriptors(mol).items():
            pdf.drawString(60, y, f"{key}: {value}"); y -= 18
        pdf.save(); buff.seek(0)
        return send_file(buff, mimetype="application/pdf", as_attachment=True, download_name="chemlab-report.pdf")
    except Exception as exc:
        return jsonify(ok=False, error=str(exc)), 400

if __name__ == "__main__":
    app.run(host="0.0.0.0", port=5000, debug=True)
