import streamlit as st
import streamlit.components.v1 as components
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors
from io import BytesIO
from datetime import date
from reportlab.lib import colors
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.enums import TA_CENTER
from reportlab.lib.units import cm
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle, PageBreak

st.set_page_config(page_title="ChemLab | Virtual Chemistry Lab", page_icon="🧪", layout="wide")

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

EXPERIMENTS = {
    "Molecular structure visualization": {
        "aim": "To represent a molecule from its SMILES notation and inspect its 2D structure and computed 3D conformer.",
        "theory": "SMILES is a line notation for describing molecular connectivity. A 2D depiction communicates atom connectivity, while a computed 3D conformer provides one plausible spatial arrangement. A computed conformer is a model and is not experimental structural evidence.",
        "procedure": "Enter or select a molecule, validate its SMILES string, generate a 2D depiction in the browser, and inspect the computed 3D model using the available display styles.",
    },
    "Molecular descriptor analysis": {
        "aim": "To calculate selected physicochemical and structural descriptors for a molecule using RDKit.",
        "theory": "Molecular descriptors summarize features of a chemical structure. Examples include molecular weight, calculated LogP, topological polar surface area (TPSA), hydrogen-bond donor/acceptor counts, rotatable bonds, heavy atoms, and ring count. These are computational estimates and should be interpreted in context.",
        "procedure": "Choose a sample molecule or enter a valid SMILES string. Run the analysis and review the descriptor table. Export the values for inclusion in a laboratory report.",
    },
    "Comparison of molecular properties": {
        "aim": "To compare calculated descriptors for two selected molecules and identify similarities and differences.",
        "theory": "Comparing the same descriptors across molecules can help explain differences in size, polarity, flexibility, and estimated hydrophobicity. Descriptor values alone do not establish biological activity or safety.",
        "procedure": "Select two molecules, calculate their descriptors with the same software method, and compare their values in the table and chart.",
    },
    "SMILES parsing and validation": {
        "aim": "To validate a molecular SMILES string and convert it into a canonical representation.",
        "theory": "SMILES encodes atoms, bonds, branches, and rings as text. A cheminformatics parser checks whether the string can be interpreted as a molecular graph. Canonical SMILES is a software-generated representation and may differ between toolkits or settings.",
        "procedure": "Enter a SMILES string, run analysis, and inspect the canonical SMILES returned by RDKit along with the calculated molecular formula and descriptors.",
    },
}

def parse_mol(smiles):
    if not isinstance(smiles, str) or not smiles.strip():
        raise ValueError("Enter a SMILES string first.")
    mol = Chem.MolFromSmiles(smiles.strip())
    if mol is None:
        raise ValueError("This SMILES string could not be parsed. Check the notation and try again.")
    return mol

def get_descriptors(mol):
    return {
        "Molecular formula": rdMolDescriptors.CalcMolFormula(mol),
        "Molecular weight (g/mol)": round(Descriptors.MolWt(mol), 3),
        "Calculated LogP": round(Descriptors.MolLogP(mol), 3),
        "TPSA (Å²)": round(rdMolDescriptors.CalcTPSA(mol), 3),
        "H-bond donors": int(rdMolDescriptors.CalcNumHBD(mol)),
        "H-bond acceptors": int(rdMolDescriptors.CalcNumHBA(mol)),
        "Rotatable bonds": int(rdMolDescriptors.CalcNumRotatableBonds(mol)),
        "Heavy atoms": int(mol.GetNumHeavyAtoms()),
        "Ring count": int(rdMolDescriptors.CalcNumRings(mol)),
    }

def make_3d_molblock(mol):
    m = Chem.AddHs(Chem.Mol(mol))
    params = AllChem.ETKDGv3()
    params.randomSeed = 42
    status = AllChem.EmbedMolecule(m, params)
    if status == 0:
        try:
            AllChem.MMFFOptimizeMolecule(m, maxIters=200)
        except Exception:
            try:
                AllChem.UFFOptimizeMolecule(m, maxIters=200)
            except Exception:
                pass
    return Chem.MolToMolBlock(m)

def html_2d(smiles):
    safe = smiles.replace("\\", "\\\\").replace("'", "\\'").replace("</", "<\\/")
    return f"""
    <div id="mol2d" style="width:100%;height:330px;display:flex;align-items:center;justify-content:center;background:#f8fbfc;border-radius:12px">
      <svg id="smiles-svg" width="540" height="320" style="max-width:100%;height:auto"></svg>
    </div>
    <script src="https://unpkg.com/smiles-drawer@2.0.3/dist/smiles-drawer.min.js"></script>
    <script>
    const smilesText = '{safe}';
    const drawer = new SmilesDrawer.SvgDrawer({{width:540,height:320,bondThickness:1.5,fontSize:14,padding:22}});
    SmilesDrawer.parse(smilesText, function(tree) {{
      drawer.draw(tree, 'smiles-svg', 'light', false);
    }}, function() {{
      document.getElementById('mol2d').innerHTML = '<p style="color:#333;padding:20px">The 2D renderer could not draw this structure. Check the SMILES string or internet connection.</p>';
    }});
    </script>
    """

def html_3d(molblock, style):
    safe = molblock.replace("\\", "\\\\").replace("`", "\\`").replace("</", "<\\/")
    styles = {
        "Ball and stick": "{stick:{radius:0.18},sphere:{scale:0.28}}",
        "Stick": "{stick:{radius:0.20}}",
        "Space-filling": "{sphere:{scale:0.85}}",
        "Wireframe": "{line:{linewidth:2}}",
    }
    style_js = styles[style]
    return f"""
    <div id="viewer3d" style="width:100%;height:390px;background:#f8fbfc;border-radius:12px"></div>
    <script src="https://3Dmol.org/build/3Dmol-min.js"></script>
    <script>
    (function(){{
      const el=document.getElementById('viewer3d');
      if(!window.$3Dmol){{el.innerHTML='<p style="padding:20px;color:#333">The 3D viewer needs an internet connection. Reload and try again.</p>';return;}}
      const viewer=$3Dmol.createViewer(el,{{backgroundColor:'#f8fbfc',antialias:true}});
      viewer.addModel(`{safe}`,'mol');
      viewer.setStyle({{}},{style_js});
      viewer.zoomTo();viewer.render();
      window.addEventListener('resize',()=>{{viewer.resize();viewer.render();}});
    }})();
    </script>
    """

def create_pdf(student_name, reg_no, course, experiment, molecule_name, smiles, canonical, props, comparison_text=""):
    buf = BytesIO()
    doc = SimpleDocTemplate(buf, pagesize=A4, rightMargin=1.8*cm, leftMargin=1.8*cm, topMargin=1.6*cm, bottomMargin=1.6*cm)
    styles = getSampleStyleSheet()
    styles.add(ParagraphStyle(name="ReportTitle", parent=styles["Title"], alignment=TA_CENTER, textColor=colors.HexColor("#123d4a"), spaceAfter=8))
    styles.add(ParagraphStyle(name="SectionHead", parent=styles["Heading2"], textColor=colors.HexColor("#087f8c"), spaceBefore=12, spaceAfter=5))
    styles.add(ParagraphStyle(name="BodySmall", parent=styles["BodyText"], leading=15, spaceAfter=6))
    story = []
    exp = EXPERIMENTS[experiment]
    story.append(Paragraph("CHEMLAB VIRTUAL LABORATORY", styles["ReportTitle"]))
    story.append(Paragraph("Computational Chemistry Practical Report", styles["Heading2"]))
    story.append(Spacer(1, 8))
    info = [
        ["Student name", student_name or "Not provided", "Registration number", reg_no or "Not provided"],
        ["Course / programme", course or "Not provided", "Date", date.today().strftime("%d %B %Y")],
        ["Experiment", experiment, "Molecule", molecule_name],
    ]
    tbl = Table(info, colWidths=[3.0*cm, 5.4*cm, 3.2*cm, 4.0*cm])
    tbl.setStyle(TableStyle([
        ("BACKGROUND",(0,0),(-1,-1),colors.HexColor("#edf5f7")),
        ("TEXTCOLOR",(0,0),(0,-1),colors.HexColor("#123d4a")),
        ("TEXTCOLOR",(2,0),(2,-1),colors.HexColor("#123d4a")),
        ("FONTNAME",(0,0),(-1,-1),"Helvetica"),
        ("FONTSIZE",(0,0),(-1,-1),8),
        ("GRID",(0,0),(-1,-1),0.4,colors.HexColor("#c5d7dc")),
        ("VALIGN",(0,0),(-1,-1),"TOP"),
        ("LEFTPADDING",(0,0),(-1,-1),6),
        ("RIGHTPADDING",(0,0),(-1,-1),6),
        ("TOPPADDING",(0,0),(-1,-1),7),
        ("BOTTOMPADDING",(0,0),(-1,-1),7),
    ]))
    story.append(tbl)
    story.append(Paragraph("1. Aim", styles["SectionHead"]))
    story.append(Paragraph(exp["aim"], styles["BodySmall"]))
    story.append(Paragraph("2. Theory", styles["SectionHead"]))
    story.append(Paragraph(exp["theory"], styles["BodySmall"]))
    story.append(Paragraph("3. Software and method", styles["SectionHead"]))
    story.append(Paragraph("The analysis was performed in ChemLab Virtual Laboratory using RDKit. The SMILES string was parsed into a molecular graph and the listed descriptors were calculated from that structure. The 3D conformer is computationally generated.", styles["BodySmall"]))
    story.append(Paragraph("4. Input structure", styles["SectionHead"]))
    story.append(Paragraph(f"<b>Input SMILES:</b> {smiles}", styles["BodySmall"]))
    story.append(Paragraph(f"<b>Canonical SMILES:</b> {canonical}", styles["BodySmall"]))
    story.append(Paragraph("5. Results", styles["SectionHead"]))
    rows = [["Descriptor", "Calculated value"]] + [[str(k), str(v)] for k,v in props.items()]
    t = Table(rows, colWidths=[8.5*cm, 7.1*cm], repeatRows=1)
    t.setStyle(TableStyle([
        ("BACKGROUND",(0,0),(-1,0),colors.HexColor("#123d4a")),
        ("TEXTCOLOR",(0,0),(-1,0),colors.white),
        ("GRID",(0,0),(-1,-1),0.4,colors.HexColor("#c5d7dc")),
        ("ROWBACKGROUNDS",(0,1),(-1,-1),[colors.white,colors.HexColor("#f1f7f8")]),
        ("FONTSIZE",(0,0),(-1,-1),9),
        ("TOPPADDING",(0,0),(-1,-1),6),
        ("BOTTOMPADDING",(0,0),(-1,-1),6),
    ]))
    story.append(t)
    if comparison_text:
        story.append(Paragraph("6. Comparison", styles["SectionHead"]))
        story.append(Paragraph(comparison_text, styles["BodySmall"]))
    story.append(Paragraph("Conclusion", styles["SectionHead"]))
    story.append(Paragraph("The selected molecular structure was processed successfully and its calculated descriptors are reported above. These computational outputs are useful for learning and comparison; they should not be treated as experimental measurements or as a prediction of biological effect.", styles["BodySmall"]))
    story.append(Spacer(1, 14))
    story.append(Paragraph("Prepared using ChemLab Virtual Laboratory · Educational use", styles["Italic"]))
    doc.build(story)
    buf.seek(0)
    return buf.getvalue()

st.markdown("""
<style>
.stApp {background:#08141f;color:#eff7fa}
.block-container {max-width:1400px;padding-top:1.6rem}
.hero {padding:1.8rem;border:1px solid #265052;border-radius:18px;background:linear-gradient(120deg,#12313d,#102130);margin-bottom:1rem}
.eyebrow {letter-spacing:.16em;color:#4ce0c1;font-size:.78rem;font-weight:700}
div[data-testid="stMetric"] {background:#101f2c;border:1px solid #233c4b;padding:1rem;border-radius:12px}
</style>
""", unsafe_allow_html=True)

st.sidebar.markdown("## 🧪 ChemLab")
st.sidebar.caption("VIRTUAL LABORATORY")
page = st.sidebar.radio("WORKSPACE", ["Practical Report & Molecular Lab", "Compare Molecules", "Theory & Learning", "Assessment"])
st.sidebar.divider()
st.sidebar.caption("DEVELOPED BY")
st.sidebar.markdown("**Sanjey Barathi M I**")
st.sidebar.caption("Artificial Intelligence and Machine Learning (AI & ML)")

st.markdown('<div class="hero"><div class="eyebrow">INTERACTIVE COMPUTATIONAL CHEMISTRY</div><h1>🧪 ChemLab Virtual Laboratory</h1><p>Explore molecular structures, visualize 2D and 3D models, calculate descriptors, and generate a student-ready practical report.</p></div>', unsafe_allow_html=True)

if page == "Practical Report & Molecular Lab":
    st.header("1 · Student and practical details")
    with st.container(border=True):
        a,b = st.columns(2)
        with a:
            student_name = st.text_input("Student name *", placeholder="Enter your full name")
            reg_no = st.text_input("Registration number *", placeholder="e.g. RA2611026050XXX")
        with b:
            course = st.text_input("Course / programme", placeholder="e.g. B.Tech Artificial Intelligence and Machine Learning")
            experiment = st.selectbox("Practical / visualization topic", list(EXPERIMENTS.keys()))
        st.caption("Enter your actual registration number carefully. The report uses the details you provide here.")
    st.header("2 · Select or enter a molecule")
    a,b = st.columns([1,1])
    with a:
        example = st.selectbox("Example molecule", ["Custom molecule"] + list(MOLECULES.keys()))
    with b:
        molecule_name = st.text_input("Molecule name", value=example if example != "Custom molecule" else "Custom molecule")
    smiles = st.text_input("SMILES notation", value=MOLECULES.get(example, "CCO"), help="Example: CCO represents ethanol.")
    model_style = st.selectbox("3D visualization style", ["Ball and stick", "Stick", "Space-filling", "Wireframe"])
    run = st.button("Generate visualization & results", type="primary", use_container_width=True)
    if run:
        try:
            mol = parse_mol(smiles)
            molblock = make_3d_molblock(mol)
            st.session_state["lab_result"] = {
                "molblock": molblock, "smiles": smiles.strip(),
                "canonical": Chem.MolToSmiles(mol),
                "props": get_descriptors(mol),
                "name": molecule_name.strip() or "Molecule",
                "experiment": experiment,
            }
        except Exception as e:
            st.error(str(e))
    result = st.session_state.get("lab_result")
    if result:
        st.divider()
        st.header("3 · Molecular visualization")
        st.markdown(f"**Experiment:** {result['experiment']}  ·  **Molecule:** {result['name']}")
        tab2d, tab3d, tabdata = st.tabs(["2D structure", "Interactive 3D model", "Descriptor results"])
        with tab2d:
            st.caption("2D structure rendered in your browser. An internet connection is needed for the drawing library.")
            components.html(html_2d(result["canonical"]), height=350, scrolling=False)
        with tab3d:
            st.caption("Rotate the model by dragging; zoom with the mouse wheel. The conformer is computationally generated.")
            components.html(html_3d(result["molblock"], model_style), height=410, scrolling=False)
        with tabdata:
            st.code(result["canonical"], language=None)
            props = result["props"]
            st.dataframe(pd.DataFrame([{"Descriptor":k,"Value":v} for k,v in props.items()]), use_container_width=True, hide_index=True)
            numeric = [(k,v) for k,v in props.items() if isinstance(v,(int,float))]
            if numeric:
                st.bar_chart(pd.DataFrame(numeric, columns=["Descriptor","Value"]).set_index("Descriptor"))
        st.divider()
        st.header("4 · Export practical report")
        if not student_name.strip() or not reg_no.strip():
            st.warning("Fill in Student name and Registration number above before downloading your final report.")
        else:
            pdf = create_pdf(student_name.strip(), reg_no.strip(), course.strip(), result["experiment"], result["name"], result["smiles"], result["canonical"], result["props"])
            filename_reg = "".join(ch for ch in reg_no.strip() if ch.isalnum() or ch in "-_")
            st.download_button("📄 Download printable practical report (PDF)", data=pdf,
                file_name=f"ChemLab_Report_{filename_reg}.pdf", mime="application/pdf", type="primary", use_container_width=True)
            csv = pd.DataFrame([{"Descriptor":k,"Value":v} for k,v in result["props"].items()]).to_csv(index=False).encode("utf-8")
            st.download_button("⬇ Download descriptor data (CSV)", data=csv, file_name="chemlab_descriptors.csv", mime="text/csv")
            st.info("Open the downloaded PDF and print it, or use your browser/PDF reader's Print command. Review your institution's report format before submission.")

elif page == "Compare Molecules":
    st.header("Molecular comparison")
    st.write("Compare calculated descriptors for two sample molecules.")
    c1,c2 = st.columns(2)
    with c1: a_name = st.selectbox("Molecule A", list(MOLECULES.keys()), index=1)
    with c2: b_name = st.selectbox("Molecule B", list(MOLECULES.keys()), index=2)
    if st.button("Compare molecules", type="primary"):
        try:
            da,db = get_descriptors(parse_mol(MOLECULES[a_name])), get_descriptors(parse_mol(MOLECULES[b_name]))
            df = pd.DataFrame({"Descriptor":list(da.keys()),a_name:list(da.values()),b_name:list(db.values())})
            st.dataframe(df, use_container_width=True, hide_index=True)
            numeric = [k for k in da if isinstance(da[k],(int,float)) and isinstance(db[k],(int,float))]
            chart = pd.DataFrame({"Descriptor":numeric,a_name:[da[k] for k in numeric],b_name:[db[k] for k in numeric]}).set_index("Descriptor")
            st.bar_chart(chart)
            st.download_button("Download comparison CSV", df.to_csv(index=False).encode(), "molecule_comparison.csv", "text/csv")
        except Exception as e: st.error(str(e))

elif page == "Theory & Learning":
    st.header("Theory and visualization guide")
    with st.expander("1. SMILES notation", expanded=True):
        st.write("SMILES is a text representation of molecular connectivity. For example, `CCO` represents ethanol. A valid SMILES string can be parsed into a molecular graph.")
    with st.expander("2. 2D molecular visualization"):
        st.write("A 2D depiction is a readable layout of atoms and bonds. It communicates connectivity, but it does not directly show the molecule's full three-dimensional geometry.")
    with st.expander("3. Interactive 3D visualization"):
        st.write("The 3D viewer displays a generated conformer. Ball-and-stick emphasizes bonds and atoms; space-filling emphasizes approximate atomic volume. Computed conformers are models, not experimental measurements.")
    with st.expander("4. Molecular descriptors"):
        st.write("Molecular weight describes molecular mass; calculated LogP estimates hydrophobicity; TPSA estimates polar surface area; H-bond donor/acceptor counts and rotatable bonds summarize selected structural features.")
    with st.expander("5. Limitations"):
        st.write("Descriptor values depend on the input structure and software implementation. They do not alone prove solubility, biological activity, toxicity, or clinical effectiveness.")

else:
    st.header("Assessment")
    qs = [
        ("What does SMILES encode?", ["Molecular connectivity as text", "Only the molecular weight", "An experimental photograph"], 0),
        ("Which descriptor estimates hydrophobicity?", ["TPSA", "Calculated LogP", "Ring count"], 1),
        ("A generated 3D conformer is…", ["A computational model", "Always the only stable shape", "A direct experimental measurement"], 0),
        ("TPSA is commonly reported in…", ["Å²", "g/mol", "degrees Celsius"], 0),
        ("Do descriptors alone prove biological activity?", ["Yes", "No"], 1),
    ]
    answers = []
    for i,(q,opts,correct) in enumerate(qs):
        choice = st.radio(f"{i+1}. {q}", opts, key=f"assessment_{i}")
        answers.append(opts.index(choice))
    if st.button("Mark assessment", type="primary"):
        score = sum(answers[i] == q[2] for i,q in enumerate(qs))
        st.success(f"Your score: {score}/{len(qs)}")

st.divider()
st.caption("CHEMLAB VIRTUAL LABORATORY · Educational project by Sanjey Barathi M I")
