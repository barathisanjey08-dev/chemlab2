import streamlit as st
import streamlit.components.v1 as components
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors
from io import BytesIO
from reportlab.lib.pagesizes import letter
from reportlab.pdfgen import canvas

st.set_page_config(page_title="ChemLab Virtual Laboratory", page_icon="🧪", layout="wide")

MOLECULES = {
    "Water": "O",
    "Ethanol": "CCO",
    "Caffeine": "Cn1c(=O)c2c(ncn2C)n(C)c1=O",
    "Aspirin": "CC(=O)Oc1ccccc1C(=O)O",
    "Glucose": "C(C1C(C(C(C(O1)O)O)O)O)O",
    "Benzene": "c1ccccc1",
    "Acetone": "CC(=O)C",
    "Ibuprofen": "CC(C)Cc1ccc(cc1)[C@@H](C)C(=O)O",
    "Paracetamol": "CC(=O)NC1=CC=C(O)C=C1",
    "Nicotine": "CN1CCC[C@H]1c2cccnc2",
}

def parse_mol(smiles):
    if not isinstance(smiles, str) or not smiles.strip():
        raise ValueError("Enter a SMILES string first.")
    mol = Chem.MolFromSmiles(smiles.strip())
    if mol is None:
        raise ValueError("That SMILES string is invalid. Please check it and try again.")
    return mol

def descriptors(mol):
    return {
        "Molecular formula": rdMolDescriptors.CalcMolFormula(mol),
        "Molecular weight (g/mol)": round(Descriptors.MolWt(mol), 3),
        "LogP": round(Descriptors.MolLogP(mol), 3),
        "TPSA (Å²)": round(rdMolDescriptors.CalcTPSA(mol), 3),
        "H-bond donors": int(rdMolDescriptors.CalcNumHBD(mol)),
        "H-bond acceptors": int(rdMolDescriptors.CalcNumHBA(mol)),
        "Rotatable bonds": int(rdMolDescriptors.CalcNumRotatableBonds(mol)),
        "Heavy atoms": int(mol.GetNumHeavyAtoms()),
        "Ring count": int(rdMolDescriptors.CalcNumRings(mol)),
    }

def molblock_3d(mol):
    m = Chem.AddHs(Chem.Mol(mol))
    params = AllChem.ETKDGv3()
    params.randomSeed = 42
    status = AllChem.EmbedMolecule(m, params)
    if status == 0:
        try:
            AllChem.MMFFOptimizeMolecule(m, maxIters=150)
        except Exception:
            pass
    return Chem.MolToMolBlock(m)

def smiles_2d_html(smiles):
    # SmilesDrawer draws the 2D structure in the browser; avoids RDKit's native drawing dependency.
    safe_smiles = smiles.replace("\\", "\\\\").replace("'", "\\'").replace("</", "<\\/")
    return f"""
    <div id="mol2d" style="width:100%;min-height:300px;display:flex;align-items:center;justify-content:center;background:#f8fbfc;border-radius:10px">
      <svg id="smiles-svg" width="520" height="330" style="max-width:100%;height:auto"></svg>
    </div>
    <script src="https://unpkg.com/smiles-drawer@2.0.3/dist/smiles-drawer.min.js"></script>
    <script>
    const s = '{safe_smiles}';
    const drawer = new SmilesDrawer.SvgDrawer({{width:520,height:330, bondThickness:1.5, fontSize:14, padding:20}});
    SmilesDrawer.parse(s, function(tree) {{
      drawer.draw(tree, 'smiles-svg', 'light', false);
    }}, function(err) {{
      document.getElementById('mol2d').innerHTML = '<p style="color:#a33">Could not draw 2D structure in this browser.</p>';
    }});
    </script>
    """

def viewer_3d_html(block, style):
    safe_block = block.replace("\\", "\\\\").replace("`", "\\`").replace("</", "<\\/")
    style_js = {
        "Ball and stick": "{stick:{radius:0.18},sphere:{scale:0.28}}",
        "Stick": "{stick:{radius:0.2}}",
        "Space-filling": "{sphere:{scale:0.9}}",
        "Wireframe": "{line:{linewidth:2}}",
    }[style]
    return f"""
    <div id="viewer3d" style="width:100%;height:390px;background:#f8fbfc;border-radius:10px;position:relative"></div>
    <script src="https://3Dmol.org/build/3Dmol-min.js"></script>
    <script>
    (function(){{
      const el=document.getElementById('viewer3d');
      if(!window.$3Dmol){{el.innerHTML='<p style="padding:20px;color:#333">3D viewer needs an internet connection. Reload this page to retry.</p>';return;}}
      const v=$3Dmol.createViewer(el,{{backgroundColor:'#f8fbfc',antialias:true}});
      v.addModel(`{safe_block}`,'mol');
      v.setStyle({{}},{style_js});
      v.zoomTo();v.render();
      window.addEventListener('resize',()=>{{v.resize();v.render();}});
    }})();
    </script>
    """

def pdf_bytes(name, smiles, props):
    buff = BytesIO()
    pdf = canvas.Canvas(buff, pagesize=letter)
    y = 745
    pdf.setTitle("ChemLab Molecular Report")
    pdf.setFont("Helvetica-Bold", 18)
    pdf.drawString(50, y, "ChemLab Molecular Report")
    y -= 32
    pdf.setFont("Helvetica", 11)
    for line in [f"Developer: Sanjey Barathi M I", f"Molecule: {name}", f"Canonical SMILES: {smiles}", "Calculated descriptors:"]:
        pdf.drawString(50, y, line)
        y -= 22
    for key, value in props.items():
        pdf.drawString(60, y, f"{key}: {value}")
        y -= 18
        if y < 55:
            pdf.showPage()
            y = 745
            pdf.setFont("Helvetica", 11)
    pdf.save()
    buff.seek(0)
    return buff.getvalue()

st.markdown("""
<style>
.stApp {background:#08141f;color:#eff7fa}
.block-container {max-width:1300px;padding-top:2rem}
h1,h2,h3 {letter-spacing:-.03em}
.hero {padding:2rem;border:1px solid #265052;border-radius:18px;background:linear-gradient(120deg,#12313d,#102130);margin-bottom:1rem}
.eyebrow {letter-spacing:.16em;color:#4ce0c1;font-size:.78rem;font-weight:700}
.small-muted {color:#9eb3bf;font-size:.9rem}
div[data-testid="stMetric"] {background:#101f2c;border:1px solid #233c4b;padding:1rem;border-radius:12px}
</style>
""", unsafe_allow_html=True)

st.sidebar.markdown("## 🧪 ChemLab")
st.sidebar.caption("VIRTUAL LABORATORY")
page = st.sidebar.radio("WORKSPACE", ["Molecular Lab", "Compare Molecules", "Theory & Learning", "Assessment"])
st.sidebar.divider()
st.sidebar.caption("DEVELOPED BY")
st.sidebar.markdown("**Sanjey Barathi M I**")
st.sidebar.caption("Artificial Intelligence & Machine Learning")

st.markdown('<div class="hero"><div class="eyebrow">INTERACTIVE COMPUTATIONAL CHEMISTRY</div><h1>🧪 ChemLab Virtual Laboratory</h1><p>Explore molecular structures, calculate chemical descriptors, and inspect molecules in an interactive 3D viewer.</p></div>', unsafe_allow_html=True)

if page == "Molecular Lab":
    st.subheader("01 · Define your molecule")
    c1, c2 = st.columns([1, 1])
    with c1:
        example = st.selectbox("Example molecules", ["Custom molecule"] + list(MOLECULES.keys()))
    default_smiles = MOLECULES.get(example, "CCO")
    with c2:
        name = st.text_input("Molecule name", value=example if example != "Custom molecule" else "Custom molecule")
    smiles = st.text_input("SMILES notation", value=default_smiles, help="Example: CCO for ethanol")
    style = st.selectbox("3D model style", ["Ball and stick", "Stick", "Space-filling", "Wireframe"])
    if st.button("Analyze molecule ↗", type="primary", use_container_width=True):
        try:
            mol = parse_mol(smiles)
            st.session_state["analysis"] = {
                "name": name or "Custom molecule",
                "smiles": Chem.MolToSmiles(mol),
                "descriptors": descriptors(mol),
                "molblock": molblock_3d(mol),
            }
        except Exception as e:
            st.error(str(e))
    result = st.session_state.get("analysis")
    if result:
        st.divider()
        st.markdown(f"### Analysis results · {result['name']}")
        st.code(result["smiles"], language=None)
        left, right = st.columns(2)
        with left:
            st.markdown("#### 2D molecular structure")
            components.html(smiles_2d_html(result["smiles"]), height=340, scrolling=False)
        with right:
            st.markdown("#### Interactive 3D model")
            components.html(viewer_3d_html(result["molblock"], style), height=410, scrolling=False)
        st.markdown("#### Molecular descriptors")
        props = result["descriptors"]
        cols = st.columns(3)
        for i, (key, value) in enumerate(props.items()):
            cols[i % 3].metric(key, value)
        df = pd.DataFrame([{"Property": k, "Value": v} for k, v in props.items()])
        st.download_button("⬇ Export CSV", df.to_csv(index=False).encode("utf-8"), "chemlab-descriptors.csv", "text/csv")
        st.download_button("⬇ Download PDF report", pdf_bytes(result["name"], result["smiles"], props),
                           "chemlab-report.pdf", "application/pdf")

elif page == "Compare Molecules":
    st.subheader("Side-by-side molecular comparison")
    col1, col2 = st.columns(2)
    with col1:
        a_name = st.selectbox("Molecule A", list(MOLECULES.keys()), index=1)
    with col2:
        b_name = st.selectbox("Molecule B", list(MOLECULES.keys()), index=2)
    if st.button("Compare properties ↗", type="primary"):
        try:
            da, db = descriptors(parse_mol(MOLECULES[a_name])), descriptors(parse_mol(MOLECULES[b_name]))
            st.dataframe(pd.DataFrame({"Descriptor": list(da.keys()), a_name: list(da.values()), b_name: list(db.values())}), use_container_width=True, hide_index=True)
        except Exception as e:
            st.error(str(e))

elif page == "Theory & Learning":
    st.subheader("Quick reference")
    with st.expander("SMILES notation", expanded=True):
        st.write("SMILES is a text notation for molecular structures. Atoms, bonds, branches, and rings are encoded as characters. For example, `CCO` represents ethanol.")
    with st.expander("Molecular weight"):
        st.write("Calculated from the atomic masses of the atoms in a molecule and commonly expressed in g/mol.")
    with st.expander("LogP"):
        st.write("An estimate of partitioning between a nonpolar phase and water; it is often used as a model descriptor of hydrophobicity.")
    with st.expander("TPSA"):
        st.write("Topological polar surface area estimates contributions from polar atoms and can support molecular comparisons.")
    with st.expander("3D conformers"):
        st.write("A computed conformer is a modelled shape, not necessarily the only shape a flexible molecule can adopt.")

else:
    st.subheader("Knowledge check")
    questions = [
        ("What does SMILES represent?", ["A text notation for molecular structure", "A unit of molecular weight", "A 3D image"], 0),
        ("Which descriptor estimates hydrophobicity?", ["TPSA", "LogP", "Ring count"], 1),
        ("A computed 3D conformer is…", ["The only possible shape", "A modelled molecular conformation", "An experimental photo"], 1),
        ("Molecular weight is commonly expressed in…", ["g/mol", "Å²", "pH"], 0),
    ]
    answers = []
    for i, (q, options, correct) in enumerate(questions):
        ans = st.radio(f"{i+1}. {q}", options, key=f"quiz_{i}")
        answers.append(options.index(ans))
    if st.button("Check answers", type="primary"):
        score = sum(answers[i] == q[2] for i, q in enumerate(questions))
        st.success(f"Score: {score}/{len(questions)}" + (" — Excellent work!" if score == len(questions) else " — Review the theory section and try again."))

st.divider()
st.caption("CHEMLAB VIRTUAL LABORATORY · Educational project by Sanjey Barathi M I")
