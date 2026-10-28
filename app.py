import os
import json
from xml.sax.saxutils import escape
from io import BytesIO
from datetime import date

import pandas as pd
import requests
import streamlit as st
import streamlit.components.v1 as components
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors

from reportlab.lib import colors
from reportlab.lib.enums import TA_CENTER
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.units import cm
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle

# ---------------------------------------------------------------------
# ChemLab Studio — educational cheminformatics workspace
# ---------------------------------------------------------------------
st.set_page_config(
    page_title="ChemLab Studio | Molecular Discovery",
    page_icon="🧪",
    layout="wide",
    initial_sidebar_state="expanded",
)

MOLECULES = {
    "Water": {"smiles": "O", "formula": "H2O", "category": "Everyday chemistry", "uses": "Drinking water, cleaning, cooling, and as a solvent in many reactions.", "plain": "Water is a small molecule made from two hydrogen atoms and one oxygen atom. Its uneven charge distribution helps it dissolve many salts and polar substances.", "caution": "This entry describes the molecule, not the safety of any particular water source."},
    "Ethanol": {"smiles": "CCO", "formula": "C2H6O", "category": "Solvent / consumer products", "uses": "Used as a solvent, in some sanitizers, laboratory work, and fuel blends.", "plain": "Ethanol is the alcohol found in alcoholic drinks. It mixes well with water and can dissolve many other substances. It is flammable and ingestion has health risks.", "caution": "Flammable. Do not ingest laboratory or industrial ethanol."},
    "Caffeine": {"smiles": "Cn1c(=O)c2c(ncn2C)n(C)c1=O", "formula": "C8H10N4O2", "category": "Food chemistry", "uses": "Present in coffee, tea, and some medicines; acts as a central nervous system stimulant.", "plain": "Caffeine can temporarily increase alertness by blocking adenosine receptors in the brain. People differ in sensitivity, and too much may cause sleep problems, anxiety, or a racing heartbeat.", "caution": "This is general educational information, not personal medical advice."},
    "Aspirin": {"smiles": "CC(=O)Oc1ccccc1C(=O)O", "formula": "C9H8O4", "category": "Pharmaceutical chemistry", "uses": "A medicine used for pain, fever, inflammation, and in some cases under medical direction to reduce clotting risk.", "plain": "Aspirin contains an aromatic ring and two oxygen-rich functional groups. Its chemical structure relates to how it is absorbed and how it interacts with enzymes in the body.", "caution": "Aspirin is not appropriate for everyone, including some children and people with bleeding risks. Follow a qualified clinician's advice."},
    "Glucose": {"smiles": "C(C1C(C(C(C(O1)O)O)O)O)O", "formula": "C6H12O6", "category": "Biochemistry", "uses": "A major energy source for cells; used in food science and laboratory solutions.", "plain": "Glucose is a simple sugar. Cells can break it down to release energy. In water, glucose mostly exists in ring-shaped forms rather than only the straight-chain form often drawn in textbooks.", "caution": "Blood-glucose interpretation is a medical matter and cannot be determined from this molecular page."},
    "Benzene": {"smiles": "c1ccccc1", "formula": "C6H6", "category": "Industrial organic chemistry", "uses": "An industrial starting material for making many chemicals and polymers.", "plain": "Benzene is a flat six-carbon ring with delocalized electrons. This electron sharing makes its bonding different from a simple ring with alternating fixed single and double bonds.", "caution": "Benzene is toxic and carcinogenic. It is not a substance for casual handling."},
    "Acetic acid": {"smiles": "CC(=O)O", "formula": "C2H4O2", "category": "Food / organic chemistry", "uses": "The acid responsible for vinegar's characteristic sourness; also used in chemical manufacturing.", "plain": "Acetic acid has a carboxyl group that can donate a proton. Household vinegar is a dilute water solution, not pure concentrated acetic acid.", "caution": "Concentrated acetic acid can cause burns; concentration matters."},
    "Acetone": {"smiles": "CC(=O)C", "formula": "C3H6O", "category": "Solvent chemistry", "uses": "A solvent used in some cleaning, manufacturing, and laboratory processes.", "plain": "Acetone contains a carbonyl group between two carbon groups. It evaporates readily and dissolves many organic substances.", "caution": "Highly flammable; avoid breathing concentrated vapour."},
    "Paracetamol (acetaminophen)": {"smiles": "CC(=O)NC1=CC=C(O)C=C1", "formula": "C8H9NO2", "category": "Pharmaceutical chemistry", "uses": "A medicine used to relieve pain and reduce fever.", "plain": "Paracetamol has an aromatic ring, a hydroxyl group, and an amide group. Its molecular structure helps determine its interactions and properties.", "caution": "Excess doses can cause severe liver damage. Follow product directions and professional advice."},
    "Ibuprofen": {"smiles": "CC(C)Cc1ccc(cc1)[C@@H](C)C(=O)O", "formula": "C13H18O2", "category": "Pharmaceutical chemistry", "uses": "A non-steroidal anti-inflammatory medicine used for pain, fever, and inflammation.", "plain": "Ibuprofen has a hydrophobic aromatic region and an acidic carboxyl group. The combination affects how it interacts with water, membranes, and biological targets.", "caution": "It can cause side effects and is unsuitable for some people. Follow professional advice."},
    "Nicotine": {"smiles": "CN1CCC[C@H]1c2cccnc2", "formula": "C10H14N2", "category": "Natural products", "uses": "A naturally occurring alkaloid found in tobacco; studied in neuroscience and pharmacology.", "plain": "Nicotine contains two nitrogen-containing ring systems. It acts on nicotinic acetylcholine receptors and can produce dependence.", "caution": "Nicotine is addictive and toxic at high doses."},
    "Lactic acid": {"smiles": "CC(O)C(=O)O", "formula": "C3H6O3", "category": "Biochemistry / food chemistry", "uses": "Produced during metabolism and used in food processing and biodegradable polymer manufacture.", "plain": "Lactic acid has both an alcohol group and a carboxylic acid group. It is associated with fermentation and is a building block for polylactic acid plastics.", "caution": "The effects of a substance depend on concentration, route, and context."},
    "Urea": {"smiles": "NC(=O)N", "formula": "CH4N2O", "category": "Biochemistry / agriculture", "uses": "The body uses urea to remove nitrogen waste; it is also widely used in fertilizers.", "plain": "Urea is a small molecule with a carbonyl group attached to two nitrogen atoms. It is highly water-soluble and important in nitrogen chemistry.", "caution": "Do not confuse a chemical description with advice on fertilizer handling or medical testing."},
    "Citric acid": {"smiles": "O=C(O)CC(O)(CC(=O)O)C(=O)O", "formula": "C6H8O7", "category": "Food chemistry", "uses": "Provides tartness in foods and drinks and is used in descaling and cleaning products.", "plain": "Citric acid has three carboxylic acid groups and one alcohol group. It can bind some metal ions and is part of the citric acid cycle in cells.", "caution": "Concentrated acids may irritate skin and eyes."},
    "Limonene": {"smiles": "CC1=CCCC(C1)C(=C)C", "formula": "C10H16", "category": "Natural products", "uses": "A fragrant compound in citrus peel oils; used in fragrances and some cleaning products.", "plain": "Limonene is a terpene with a mostly carbon-and-hydrogen structure, so it is relatively hydrophobic and has a citrus-like smell.", "caution": "Oxidized limonene may cause skin sensitization in some people."},
    "Benzaldehyde": {"smiles": "O=Cc1ccccc1", "formula": "C7H6O", "category": "Organic chemistry", "uses": "Used in fragrance, flavour chemistry, and as an organic synthesis intermediate.", "plain": "Benzaldehyde is a benzene ring with an aldehyde group attached. The carbonyl group gives it different reactivity from benzene itself.", "caution": "Use appropriate lab safety procedures for chemical reagents."},
    "Glycine": {"smiles": "NCC(=O)O", "formula": "C2H5NO2", "category": "Biochemistry", "uses": "The simplest amino acid; a building block of proteins and a component in biochemical systems.", "plain": "Glycine has an amino group and a carboxyl group. In water, amino acids often exist as zwitterions, with both positive and negative charges.", "caution": "The displayed neutral SMILES is a simplified representation of its protonation state."},
    "ATP (simplified)": {"smiles": "Nc1ncnc2n(cnc12)[C@@H]1O[C@H](COP(=O)(O)OP(=O)(O)OP(=O)(O)O)[C@@H](O)[C@H]1O", "formula": "C10H16N5O13P3", "category": "Biochemistry", "uses": "A central energy-transfer molecule in cells.", "plain": "ATP contains adenine, ribose, and a chain of phosphate groups. Cells couple ATP breakdown to many processes that require energy.", "caution": "The simplified representation does not capture all biological ionization states."},
}

EXPERIMENTS = {
    "Molecular structure and 3D visualization": {
        "aim": "To represent a selected molecule from its SMILES notation and inspect its 2D connectivity and generated 3D conformer.",
        "theory": "SMILES is a text notation for molecular connectivity. A 2D depiction emphasizes atom connectivity; a 3D conformer is a computationally generated spatial arrangement. A generated conformer is a model, not an experimental structure.",
        "procedure": "Select or enter a molecule, validate its SMILES, generate its structure, and inspect the 2D depiction and interactive 3D model.",
    },
    "Molecular descriptor analysis": {
        "aim": "To calculate selected molecular descriptors from a valid molecular structure.",
        "theory": "Descriptors summarize properties such as molecular mass, calculated hydrophobicity, polar surface area, hydrogen-bond counts, flexibility, and ring content. They help compare molecules but do not alone establish biological activity or safety.",
        "procedure": "Choose a molecule, calculate its descriptors, review the table and graph, and export the results.",
    },
    "Structure–property comparison": {
        "aim": "To compare calculated descriptors for two molecular structures using a consistent method.",
        "theory": "Differences in size, polarity, flexibility, and hydrophobicity can help explain different behaviour in chemical systems. Descriptors are estimates and should be interpreted with context.",
        "procedure": "Choose two molecules, compute their descriptors, inspect the comparison table and chart, and record the major differences.",
    },
    "SMILES validation and canonicalization": {
        "aim": "To validate a SMILES string and obtain its canonical representation and molecular formula.",
        "theory": "A SMILES parser converts textual atom and bond notation into a molecular graph. Canonical SMILES is a toolkit-generated representation and may vary between software implementations.",
        "procedure": "Enter a SMILES string, validate it, and record the canonical representation and molecular descriptors.",
    },
}

def parse_molecule(smiles):
    if not smiles or not str(smiles).strip():
        raise ValueError("Enter a SMILES string first.")
    mol = Chem.MolFromSmiles(str(smiles).strip())
    if mol is None:
        raise ValueError("RDKit could not parse this SMILES. Please check the notation.")
    return mol

def descriptors(mol):
    return {
        "Molecular formula": rdMolDescriptors.CalcMolFormula(mol),
        "Molecular weight (g/mol)": round(Descriptors.MolWt(mol), 3),
        "Calculated LogP": round(Descriptors.MolLogP(mol), 3),
        "Topological polar surface area (Å²)": round(rdMolDescriptors.CalcTPSA(mol), 3),
        "Hydrogen-bond donors": int(rdMolDescriptors.CalcNumHBD(mol)),
        "Hydrogen-bond acceptors": int(rdMolDescriptors.CalcNumHBA(mol)),
        "Rotatable bonds": int(rdMolDescriptors.CalcNumRotatableBonds(mol)),
        "Heavy atoms": int(mol.GetNumHeavyAtoms()),
        "Ring count": int(rdMolDescriptors.CalcNumRings(mol)),
    }

def molblock_3d(mol):
    m = Chem.AddHs(Chem.Mol(mol))
    params = AllChem.ETKDGv3()
    params.randomSeed = 17
    status = AllChem.EmbedMolecule(m, params)
    if status == 0:
        try:
            AllChem.MMFFOptimizeMolecule(m, maxIters=250)
        except Exception:
            try:
                AllChem.UFFOptimizeMolecule(m, maxIters=250)
            except Exception:
                pass
    return Chem.MolToMolBlock(m)

def viewer_2d(smiles):
    safe = smiles.replace("\\", "\\\\").replace("'", "\\'").replace("</", "<\\/")
    return f"""
    <div id="sdraw" style="height:360px;width:100%;display:flex;align-items:center;justify-content:center;background:#fff;border-radius:16px">
      <svg id="chem-svg" width="620" height="340" style="max-width:100%;height:auto"></svg>
    </div>
    <script src="https://unpkg.com/smiles-drawer@2.0.3/dist/smiles-drawer.min.js"></script>
    <script>
    const sd = new SmilesDrawer.SvgDrawer({{width:620,height:340,bondThickness:1.6,fontSize:15,padding:28}});
    SmilesDrawer.parse('{safe}', t => sd.draw(t,'chem-svg','light',false),
      () => document.getElementById('sdraw').innerHTML='<p style="color:#333;padding:20px">2D drawing unavailable. Check your connection.</p>');
    </script>
    """

def viewer_3d(block, style, surface):
    safe = block.replace("\\", "\\\\").replace("`", "\\`").replace("</", "<\\/")
    styles = {
        "Ball and stick": "{stick:{radius:0.16},sphere:{scale:0.27}}",
        "Stick": "{stick:{radius:0.19}}",
        "Space filling": "{sphere:{scale:0.8}}",
        "Wireframe": "{line:{linewidth:2}}",
    }
    surface_js = "viewer.addSurface($3Dmol.SurfaceType.VDW,{opacity:0.55,colorscheme:'Jmol'});" if surface else ""
    return f"""
    <div id="view3d" style="height:430px;width:100%;background:#fff;border-radius:16px"></div>
    <script src="https://3Dmol.org/build/3Dmol-min.js"></script>
    <script>
    (function(){{
      const el=document.getElementById('view3d');
      if(!window.$3Dmol){{el.innerHTML='<p style="color:#333;padding:20px">3D viewer needs internet access. Reload and try again.</p>';return;}}
      const viewer=$3Dmol.createViewer(el,{{backgroundColor:'#ffffff',antialias:true}});
      viewer.addModel(`{safe}`,'mol');
      viewer.setStyle({{}},{styles[style]});
      {surface_js}
      viewer.zoomTo();viewer.render();
      window.addEventListener('resize',()=>{{viewer.resize();viewer.render();}});
    }})();
    </script>
    """

def pdf_report(student, reg, course, experiment, mol_name, smiles, canonical, props, note=""):
    """Create a polished, printable A4 practical report with student details."""
    buffer = BytesIO()
    doc = SimpleDocTemplate(
        buffer, pagesize=A4,
        leftMargin=1.65*cm, rightMargin=1.65*cm,
        topMargin=2.15*cm, bottomMargin=1.8*cm,
        title=f"ChemLab Studio Practical Report - {mol_name}",
        author="ChemLab Studio",
        subject="Educational cheminformatics practical report",
    )
    styles = getSampleStyleSheet()
    styles.add(ParagraphStyle(
        name="ChemTitle", parent=styles["Title"], alignment=TA_CENTER,
        textColor=colors.HexColor("#0B6670"), fontName="Helvetica-Bold",
        fontSize=20, leading=24, spaceAfter=4,
    ))
    styles.add(ParagraphStyle(
        name="ChemSubtitle", parent=styles["Normal"], alignment=TA_CENTER,
        textColor=colors.HexColor("#526777"), fontSize=9, leading=12,
        spaceAfter=13,
    ))
    styles.add(ParagraphStyle(
        name="ChemHead", parent=styles["Heading2"],
        textColor=colors.HexColor("#087E8B"), fontName="Helvetica-Bold",
        fontSize=11, leading=14, spaceBefore=11, spaceAfter=5,
        keepWithNext=True,
    ))
    styles.add(ParagraphStyle(
        name="ChemBody", parent=styles["BodyText"],
        textColor=colors.HexColor("#263746"), fontSize=9.2,
        leading=13.2, spaceAfter=5,
    ))
    styles.add(ParagraphStyle(
        name="ChemSmall", parent=styles["BodyText"],
        textColor=colors.HexColor("#526777"), fontSize=7.7, leading=10,
    ))
    exp = EXPERIMENTS[experiment]
    # Escape all user-controlled text before placing it in ReportLab Paragraph markup.
    safe = lambda value: escape(str(value or "Not provided"))
    story = [
        Paragraph("CHEMLAB STUDIO", styles["ChemTitle"]),
        Paragraph("MOLECULAR DISCOVERY LABORATORY  |  PRACTICAL REPORT", styles["ChemSubtitle"]),
    ]
    info = [
        [Paragraph("<b>STUDENT NAME</b>", styles["ChemSmall"]), Paragraph(safe(student), styles["ChemBody"]),
         Paragraph("<b>REGISTRATION NO.</b>", styles["ChemSmall"]), Paragraph(safe(reg), styles["ChemBody"])],
        [Paragraph("<b>PROGRAMME / COURSE</b>", styles["ChemSmall"]), Paragraph(safe(course), styles["ChemBody"]),
         Paragraph("<b>REPORT DATE</b>", styles["ChemSmall"]), Paragraph(date.today().strftime("%d %B %Y"), styles["ChemBody"])],
        [Paragraph("<b>PRACTICAL</b>", styles["ChemSmall"]), Paragraph(safe(experiment), styles["ChemBody"]),
         Paragraph("<b>MOLECULE</b>", styles["ChemSmall"]), Paragraph(safe(mol_name), styles["ChemBody"])],
    ]
    info_tbl = Table(info, colWidths=[2.7*cm, 5.1*cm, 2.8*cm, 5.0*cm], hAlign="LEFT")
    info_tbl.setStyle(TableStyle([
        ("BACKGROUND", (0,0), (-1,-1), colors.HexColor("#F1F7F8")),
        ("GRID", (0,0), (-1,-1), 0.45, colors.HexColor("#C7D9DD")),
        ("VALIGN", (0,0), (-1,-1), "MIDDLE"),
        ("LEFTPADDING", (0,0), (-1,-1), 7), ("RIGHTPADDING", (0,0), (-1,-1), 7),
        ("TOPPADDING", (0,0), (-1,-1), 7), ("BOTTOMPADDING", (0,0), (-1,-1), 4),
    ]))
    story.extend([
        info_tbl,
        Paragraph("1. AIM", styles["ChemHead"]), Paragraph(safe(exp["aim"]), styles["ChemBody"]),
        Paragraph("2. THEORY", styles["ChemHead"]), Paragraph(safe(exp["theory"]), styles["ChemBody"]),
        Paragraph("3. PROCEDURE", styles["ChemHead"]), Paragraph(safe(exp["procedure"]), styles["ChemBody"]),
        Paragraph("4. MOLECULAR INPUT", styles["ChemHead"]),
        Paragraph(f"<b>Molecule:</b> {safe(mol_name)}", styles["ChemBody"]),
        Paragraph(f"<b>Input SMILES:</b> <font name='Courier'>{safe(smiles)}</font>", styles["ChemBody"]),
        Paragraph(f"<b>Canonical SMILES:</b> <font name='Courier'>{safe(canonical)}</font>", styles["ChemBody"]),
        Paragraph("5. COMPUTED RESULTS", styles["ChemHead"]),
    ])
    rows = [[Paragraph("<b>MOLECULAR DESCRIPTOR</b>", styles["ChemSmall"]), Paragraph("<b>CALCULATED VALUE</b>", styles["ChemSmall"])]]
    rows += [[Paragraph(safe(k), styles["ChemBody"]), Paragraph(safe(v), styles["ChemBody"])] for k, v in props.items()]
    result_table = Table(rows, colWidths=[9.2*cm, 6.3*cm], repeatRows=1, hAlign="LEFT")
    result_table.setStyle(TableStyle([
        ("BACKGROUND", (0,0), (-1,0), colors.HexColor("#0B6670")),
        ("TEXTCOLOR", (0,0), (-1,0), colors.white),
        ("GRID", (0,0), (-1,-1), 0.4, colors.HexColor("#C5D9DC")),
        ("ROWBACKGROUNDS", (0,1), (-1,-1), [colors.white, colors.HexColor("#F3F8F9")]),
        ("VALIGN", (0,0), (-1,-1), "MIDDLE"),
        ("LEFTPADDING", (0,0), (-1,-1), 7), ("RIGHTPADDING", (0,0), (-1,-1), 7),
        ("TOPPADDING", (0,0), (-1,-1), 3), ("BOTTOMPADDING", (0,0), (-1,-1), 1),
    ]))
    story.append(result_table)
    if note.strip():
        story.extend([Paragraph("6. OBSERVATIONS & INTERPRETATION", styles["ChemHead"]), Paragraph(safe(note).replace("\n", "<br/>"), styles["ChemBody"])])
    story.extend([
        Paragraph("CONCLUSION", styles["ChemHead"]),
        Paragraph("The molecular structure was parsed and selected descriptors were calculated computationally. These values support educational analysis and comparison; they are not experimental measurements and do not independently establish biological activity, safety, or suitability for use.", styles["ChemBody"]),
        Spacer(1, 9),
        Paragraph("ACADEMIC INTEGRITY NOTE", styles["ChemHead"]),
        Paragraph("Review the calculations and interpretation, add any institution-required details, and ensure the report accurately reflects your own work before submission.", styles["ChemSmall"]),
    ])
    def page_decor(canvas, document):
        canvas.saveState()
        width, height = A4
        canvas.setFillColor(colors.HexColor("#0B6670"))
        canvas.rect(0, height - 0.55*cm, width, 0.55*cm, fill=1, stroke=0)
        canvas.setStrokeColor(colors.HexColor("#D6E4E7"))
        canvas.line(doc.leftMargin, 1.25*cm, width-doc.rightMargin, 1.25*cm)
        canvas.setFont("Helvetica", 7.5)
        canvas.setFillColor(colors.HexColor("#526777"))
        canvas.drawString(doc.leftMargin, 0.88*cm, "CHEMLAB STUDIO  •  EDUCATIONAL CHEMINFORMATICS REPORT")
        canvas.drawRightString(width-doc.rightMargin, 0.88*cm, f"Page {document.page}")
        canvas.restoreState()
    doc.build(story, onFirstPage=page_decor, onLaterPages=page_decor)
    buffer.seek(0)
    return buffer.getvalue()

def openai_answer(history, prompt):
    """Call the OpenAI Responses API. The API key is read only from Streamlit secrets/env."""
    try:
        api_key = st.secrets.get("OPENAI_API_KEY", os.getenv("OPENAI_API_KEY", ""))
    except Exception:
        api_key = os.getenv("OPENAI_API_KEY", "")
    if not api_key:
        return None, "OpenAI is not configured. Add OPENAI_API_KEY in Streamlit → Manage app → Settings → Secrets."
    try:
        from openai import OpenAI
        client = OpenAI(api_key=api_key, timeout=35.0, max_retries=1)
        conversation = []
        for item in history[-8:]:
            role = "assistant" if item.get("role") == "assistant" else "user"
            conversation.append({"role": role, "content": str(item.get("content", ""))})
        conversation.append({"role": "user", "content": prompt})
        response = client.responses.create(
            model="gpt-4.1-mini",
            instructions=(
                "You are ChemLab Guide, a patient educational chemistry and cheminformatics tutor. "
                "Explain molecular structures, SMILES, molecular weight, calculated LogP, TPSA, "
                "hydrogen bonding, and chemistry coursework in clear, student-friendly language. "
                "For molecules, discuss structure, properties, common uses, and relevant safety caveats. "
                "Distinguish computational estimates from experimental facts. Do not provide medical diagnoses "
                "or unsafe chemical handling instructions. Use headings and examples when useful."
            ),
            input=conversation,
            max_output_tokens=1200,
        )
        answer = getattr(response, "output_text", "")
        if not answer:
            return None, "OpenAI returned an empty answer. Please try asking the question again."
        return answer, None
    except Exception as e:
        # Show a concise diagnostic but never expose the API key.
        return None, f"OpenAI API error: {str(e)[:350]}"

# -------------------- Styling --------------------
st.markdown("""
<style>
:root { color-scheme: dark; }
.stApp { background: #08121D; color: #EAF4F7; }
.block-container { max-width: 1440px; padding-top: 1.4rem; padding-bottom: 3rem; }
[data-testid="stSidebar"] { background: #101D2A; border-right: 1px solid #243C4D; }
.hero { padding: 1.6rem 1.8rem; border-radius: 20px; border: 1px solid #24505A; background: linear-gradient(115deg,#123543,#102233 68%,#17283A); margin-bottom: 1.2rem; }
.eyebrow { color: #4CE0C1; letter-spacing: .17em; font-weight: 800; font-size: .78rem; }
.hero p { color: #C6DCE3; font-size: 1.05rem; }
.mini { color: #9FB9C7; font-size: .92rem; }
div[data-testid="stMetric"] { background: #102131; border: 1px solid #294354; border-radius: 14px; padding: .9rem; }
div.stButton > button[kind="primary"], div.stDownloadButton > button[kind="primary"] { border-radius: 10px; }
</style>
""", unsafe_allow_html=True)

# -------------------- Sidebar --------------------
st.sidebar.markdown("# 🧪 ChemLab Studio")
st.sidebar.caption("MOLECULAR DISCOVERY WORKSPACE")
page = st.sidebar.radio("Navigate", [
    "Overview",
    "Molecular Explorer",
    "Plain-language Molecule Guide",
    "Property Comparison",
    "Practical Report Builder",
    "AI Chemistry Assistant",
    "Theory & Mini Quiz",
])
st.sidebar.divider()
st.sidebar.markdown("### 👨‍🎓 Student profile")
student_name = st.sidebar.text_input("Student name", key="student_name")
registration = st.sidebar.text_input("Registration number", key="registration", placeholder="e.g. RA2611026050XXX")
course = st.sidebar.text_input("Programme / course", key="course", placeholder="e.g. B.Tech AI & ML")
st.sidebar.caption("Developed by **Sanjey Barathi M I**")
st.sidebar.caption("Artificial Intelligence and Machine Learning (AI & ML)")
st.sidebar.divider()
st.sidebar.caption("Educational tool · calculated values are estimates")

# session state defaults
if "current_molecule" not in st.session_state:
    st.session_state.current_molecule = "Ethanol"
if "current_smiles" not in st.session_state:
    st.session_state.current_smiles = MOLECULES["Ethanol"]["smiles"]
if "ai_history" not in st.session_state:
    st.session_state.ai_history = []

st.markdown('<div class="hero"><div class="eyebrow">CHEMISTRY · VISUALIZATION · DISCOVERY</div><h1>🧪 ChemLab Studio</h1><p>A student-friendly molecular lab with interactive structures, calculated properties, plain-language explanations, and printable practical reports.</p></div>', unsafe_allow_html=True)

if page == "Overview":
    st.subheader("Your molecular workspace")
    c1,c2,c3 = st.columns(3)
    c1.metric("Molecule library", f"{len(MOLECULES)}", "Curated examples")
    c2.metric("Visualization modes", "2D + 3D", "Browser-based viewers")
    c3.metric("Report export", "PDF + CSV", "Printable results")
    st.markdown("### Start exploring")
    left,right = st.columns([1.1, .9])
    with left:
        st.markdown("#### 🧬 Explore a molecule")
        selected = st.selectbox("Choose a molecule", list(MOLECULES.keys()), index=list(MOLECULES.keys()).index(st.session_state.current_molecule) if st.session_state.current_molecule in MOLECULES else 1)
        info = MOLECULES[selected]
        st.write(f"**Category:** {info['category']}")
        st.write(f"**Formula:** {info['formula']}")
        st.write(info["plain"])
        if st.button("Load molecule in Explorer", type="primary"):
            st.session_state.current_molecule = selected
            st.session_state.current_smiles = info["smiles"]
            st.session_state["explorer_selection"] = selected
            st.success("Loaded. Open Molecular Explorer from the sidebar.")
    with right:
        st.markdown("#### 📘 What you can do")
        st.write("• Validate SMILES and inspect canonical notation.")
        st.write("• View a 2D structure and interactive 3D conformer.")
        st.write("• Calculate molecular descriptors and compare compounds.")
        st.write("• Learn properties and everyday uses in simple language.")
        st.write("• Export a student-labelled PDF practical report.")
        st.write("• Ask the AI assistant for help when an OpenAI API key is configured.")
    st.info("Tip: Fill in your student profile in the sidebar before creating your final PDF report.")

elif page == "Molecular Explorer":
    st.header("🧬 Molecular Explorer")
    st.write("Choose a molecule from the library or enter your own SMILES. The same structure powers the visualizations and descriptor analysis.")
    molecule_names = list(MOLECULES.keys())
    default_name = st.session_state.get("explorer_selection", st.session_state.current_molecule)
    if default_name not in molecule_names:
        default_name = "Ethanol"
    selected = st.selectbox("Molecule library", ["Custom SMILES"] + molecule_names, index=(["Custom SMILES"] + molecule_names).index(default_name) if default_name in molecule_names else 0)
    if selected != "Custom SMILES":
        default_smiles = MOLECULES[selected]["smiles"]
        st.session_state.current_molecule = selected
    else:
        default_smiles = st.session_state.current_smiles
    smiles = st.text_input("SMILES notation", value=default_smiles, key="explorer_smiles")
    v1,v2,v3 = st.columns([1,1,1])
    with v1:
        style = st.selectbox("3D display style", ["Ball and stick", "Stick", "Space filling", "Wireframe"])
    with v2:
        surface = st.checkbox("Show van der Waals surface")
    with v3:
        st.write("")
        st.write("")
        analyze = st.button("Generate structure & properties", type="primary", use_container_width=True)
    if analyze:
        try:
            mol = parse_molecule(smiles)
            st.session_state.current_smiles = smiles.strip()
            st.session_state.current_molecule = selected if selected != "Custom SMILES" else "Custom molecule"
            st.session_state.explorer_result = {
                "name": selected if selected != "Custom SMILES" else "Custom molecule",
                "smiles": smiles.strip(),
                "canonical": Chem.MolToSmiles(mol),
                "props": descriptors(mol),
                "molblock": molblock_3d(mol),
            }
        except Exception as e:
            st.error(str(e))
    result = st.session_state.get("explorer_result")
    if result:
        st.success("Molecular structure validated.")
        st.caption(f"Canonical SMILES: {result['canonical']}")
        tab1,tab2,tab3 = st.tabs(["2D structure", "Interactive 3D model", "Calculated descriptors"])
        with tab1:
            components.html(viewer_2d(result["canonical"]), height=375, scrolling=False)
        with tab2:
            st.caption("Drag to rotate and use the scroll wheel or pinch to zoom. This is a generated conformer, not an experimental structure.")
            components.html(viewer_3d(result["molblock"], style, surface), height=445, scrolling=False)
        with tab3:
            props = result["props"]
            metrics = st.columns(3)
            metrics[0].metric("Molecular weight", f"{props['Molecular weight (g/mol)']:.3f} g/mol")
            metrics[1].metric("Calculated LogP", props["Calculated LogP"])
            metrics[2].metric("TPSA", f"{props['Topological polar surface area (Å²)']} Å²")
            st.dataframe(pd.DataFrame([{"Descriptor": k, "Value": v} for k,v in props.items()]), use_container_width=True, hide_index=True)
            numeric = [(k,v) for k,v in props.items() if isinstance(v,(int,float))]
            st.bar_chart(pd.DataFrame(numeric, columns=["Descriptor","Value"]).set_index("Descriptor"))
        st.session_state.report_molecule = result
        st.download_button("Download descriptor CSV", pd.DataFrame([{"Descriptor":k,"Value":v} for k,v in result["props"].items()]).to_csv(index=False).encode(), "molecular_descriptors.csv", "text/csv")

elif page == "Plain-language Molecule Guide":
    st.header("🌱 Understand a molecule in plain language")
    st.write("This page translates chemistry into everyday language. It is educational information, not a safety approval or medical recommendation.")
    name = st.selectbox("Select a molecule", list(MOLECULES.keys()))
    info = MOLECULES[name]
    mol = parse_molecule(info["smiles"])
    props = descriptors(mol)
    a,b = st.columns([1,1])
    with a:
        st.markdown(f"### {name}")
        st.markdown(f"**Formula:** `{props['Molecular formula']}`")
        st.markdown(f"**Category:** {info['category']}")
        st.markdown("#### What is it?")
        st.write(info["plain"])
        st.markdown("#### Where is it used?")
        st.write(info["uses"])
        st.markdown("#### Important caution")
        st.warning(info["caution"])
    with b:
        st.markdown("#### At a glance")
        st.dataframe(pd.DataFrame([{"Property":k,"Value":v} for k,v in props.items()]), use_container_width=True, hide_index=True)
        components.html(viewer_2d(Chem.MolToSmiles(mol)), height=365, scrolling=False)
    with st.expander("What do the properties mean?"):
        st.markdown("""
- **Molecular weight:** approximate mass of one mole of the molecule.
- **Calculated LogP:** a computational estimate of how the neutral molecule partitions between oily and watery environments; it is not a direct solubility measurement.
- **TPSA:** estimated surface area associated with polar atoms; it can help compare molecules but is not a stand-alone biological prediction.
- **Hydrogen-bond donors / acceptors:** counts of structural features that may participate in hydrogen bonding.
- **Rotatable bonds:** a rough indicator of molecular flexibility.
- **Ring count:** number of rings detected by the cheminformatics toolkit.
""")

elif page == "Property Comparison":
    st.header("📊 Compare molecular properties")
    names = list(MOLECULES.keys())
    col1,col2 = st.columns(2)
    with col1: name_a = st.selectbox("Molecule A", names, index=names.index("Ethanol"))
    with col2: name_b = st.selectbox("Molecule B", names, index=names.index("Caffeine"))
    da = descriptors(parse_molecule(MOLECULES[name_a]["smiles"]))
    db = descriptors(parse_molecule(MOLECULES[name_b]["smiles"]))
    rows = []
    for key in da:
        va,vb = da[key],db[key]
        if isinstance(va,(int,float)) and isinstance(vb,(int,float)):
            diff = round(vb-va,3)
        else:
            diff = "—"
        rows.append({"Descriptor":key,name_a:va,name_b:vb,"B − A":diff})
    df = pd.DataFrame(rows)
    st.dataframe(df, use_container_width=True, hide_index=True)
    numeric_keys = [k for k in da if isinstance(da[k],(int,float)) and isinstance(db[k],(int,float))]
    chart = pd.DataFrame({"Descriptor":numeric_keys,name_a:[da[k] for k in numeric_keys],name_b:[db[k] for k in numeric_keys]}).set_index("Descriptor")
    st.markdown("### Descriptor comparison")
    st.bar_chart(chart)
    st.download_button("Download comparison CSV", df.to_csv(index=False).encode(), "chemlab_comparison.csv", "text/csv")
    st.caption("Values are calculated with the same RDKit descriptor functions. They do not by themselves predict safety, effectiveness, or biological activity.")

elif page == "Practical Report Builder":
    st.header("📄 Practical Report Builder")
    st.write("Create a polished, print-ready practical report. Enter the student details below; they will be printed in the report header.")
    student_col, reg_col = st.columns(2)
    with student_col:
        report_student = st.text_input("Student name", value=st.session_state.get("student_name", ""), key="report_student_name", placeholder="Enter full name")
    with reg_col:
        report_registration = st.text_input("Registration number", value=st.session_state.get("registration", ""), key="report_registration", placeholder="Enter registration number")
    report_course = st.text_input("Programme / course", value=st.session_state.get("course", ""), key="report_course")
    st.divider()
    report_name = st.text_input("Molecule name", value=st.session_state.get("current_molecule","Ethanol"), key="report_name")
    report_smiles = st.text_input("SMILES for report", value=st.session_state.get("current_smiles","CCO"), key="report_smiles")
    experiment = st.selectbox("Practical topic", list(EXPERIMENTS.keys()))
    interpretation = st.text_area("Your observations / interpretation (optional)", placeholder="Write what you noticed in the visualization or descriptor results.")
    if st.button("Prepare PDF report", type="primary"):
        try:
            mol = parse_molecule(report_smiles)
            props = descriptors(mol)
            canonical = Chem.MolToSmiles(mol)
            pdf = pdf_report(report_student.strip(), report_registration.strip(), report_course.strip(), experiment, report_name.strip(), report_smiles.strip(), canonical, props, interpretation.strip())
            st.session_state.pdf_bytes = pdf
            st.session_state.pdf_filename = "ChemLab_Practical_Report_" + ("".join(c for c in report_registration if c.isalnum() or c in "-_") or "Student") + ".pdf"
            st.success("Your PDF report is ready.")
            st.dataframe(pd.DataFrame([{"Descriptor":k,"Value":v} for k,v in props.items()]), use_container_width=True, hide_index=True)
        except Exception as e:
            st.error(str(e))
    if "pdf_bytes" in st.session_state:
        st.download_button("⬇️ Download printable PDF", st.session_state.pdf_bytes, st.session_state.pdf_filename, "application/pdf", type="primary", use_container_width=True)
        st.info("Open the downloaded PDF and print it, or choose Print in your PDF viewer. Confirm your college's report requirements before submission.")

elif page == "AI Chemistry Assistant":
    st.header("🤖 ChemLab AI Guide")
    st.write("Ask about molecular structures, SMILES, descriptors, organic chemistry, chemistry coursework, or how to interpret your report.")
    try:
        api_configured = bool(st.secrets.get("OPENAI_API_KEY", os.getenv("OPENAI_API_KEY", "")))
    except Exception:
        api_configured = bool(os.getenv("OPENAI_API_KEY", ""))
    if api_configured:
        st.success("OpenAI API key detected. Responses will be generated by the configured GPT model.")
    else:
        st.warning("OpenAI setup needed: add an API key in Streamlit app settings. A key is not included in the app code.")
    with st.expander("Connect OpenAI (one-time setup)"):
        st.markdown("""
1. Create an API key in the [OpenAI API Platform](https://platform.openai.com/api-keys).
2. Open this app's **Manage app → Settings → Secrets**.
3. Add the following TOML line with your own key:
```toml
OPENAI_API_KEY = "YOUR_ACTUAL_OPENAI_API_KEY"
```
4. Save the secret and allow the app to restart.

**Keep the key private.** ChatGPT subscriptions and API billing are separate; API usage may require billing to be enabled.
""")
    for msg in st.session_state.ai_history:
        with st.chat_message(msg["role"]):
            st.markdown(msg["content"])
    prompt = st.chat_input("Ask ChemLab Guide a question…")
    if prompt:
        st.session_state.ai_history.append({"role": "user", "content": prompt})
        with st.chat_message("user"):
            st.markdown(prompt)
        answer, error = openai_answer(st.session_state.ai_history[:-1], prompt)
        if answer is None:
            fallback = ("I couldn't reach the GPT service for this request. Check the API key, API billing, model access, and then try again. "
                        "You can still explore the Plain-language Molecule Guide and Theory & Mini Quiz pages.")
            if any(w in prompt.lower() for w in ["logp", "hydrophobic"]):
                fallback = "Calculated LogP estimates how a neutral molecule partitions between an oily phase and water. A higher value generally indicates greater hydrophobic character, but it is not a direct measurement of solubility."
            elif any(w in prompt.lower() for w in ["smiles", "canonical"]):
                fallback = "SMILES is a text notation for atoms and bonds. RDKit parses it into a molecular graph and can generate a canonical SMILES representation. Different toolkits can produce different canonical strings for the same molecule."
            elif any(w in prompt.lower() for w in ["tpsa", "polar surface"]):
                fallback = "TPSA is a calculated descriptor that estimates polar surface area associated with polar atoms. It can help compare molecules, but it is not a direct measurement of absorption or solubility."
            st.caption(error or "The AI service did not return an answer.")
            answer = fallback
        st.session_state.ai_history.append({"role": "assistant", "content": answer})
        with st.chat_message("assistant"):
            st.markdown(answer)
    if st.button("Clear chat history"):
        st.session_state.ai_history = []
        st.rerun()

else:
    st.header("📚 Theory & Mini Quiz")
    with st.expander("1. What is cheminformatics?", expanded=True):
        st.write("Cheminformatics combines chemistry, data, and computation to store, represent, search, compare, and analyse chemical structures.")
    with st.expander("2. How does a SMILES string work?"):
        st.write("SMILES describes atoms and bonds in text. For example, `CCO` represents ethanol. A parser turns that text into a molecular graph.")
    with st.expander("3. 2D versus 3D structures"):
        st.write("A 2D depiction shows connectivity clearly. A 3D conformer gives a possible spatial arrangement. Flexible molecules may have many conformations, and a computed conformer is not necessarily the experimentally preferred one.")
    with st.expander("4. Molecular descriptors"):
        st.write("Descriptors summarise a molecular structure. Examples include molecular weight, calculated LogP, TPSA, hydrogen-bond donor/acceptor counts, rotatable bonds, and ring count.")
    st.markdown("### Quick knowledge check")
    quiz = [
        ("What does SMILES mainly represent?", ["Molecular connectivity in text", "A measured melting point", "A microscope image"], 0),
        ("Which descriptor is commonly used as an estimate of hydrophobicity?", ["Ring count", "Calculated LogP", "Hydrogen-bond donors"], 1),
        ("A generated 3D conformer is best described as…", ["Always an experimental structure", "A computational model", "A safety certificate"], 1),
        ("Does a molecular descriptor alone prove a medicine is safe?", ["Yes", "No"], 1),
    ]
    answers = []
    for i,(q,opts,correct) in enumerate(quiz):
        answers.append(opts.index(st.radio(f"{i+1}. {q}", opts, key=f"quiz_{i}")))
    if st.button("Check answers", type="primary"):
        score = sum(answers[i] == q[2] for i,q in enumerate(quiz))
        st.success(f"Score: {score}/{len(quiz)}")
        if score == len(quiz):
            st.balloons()
        else:
            st.info("Review the theory above and try again.")

st.divider()
st.caption("ChemLab Studio · Developed by Sanjey Barathi M I · Educational use only. Computational estimates are not experimental measurements.")
