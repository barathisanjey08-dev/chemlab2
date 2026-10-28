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
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle, PageBreak

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

def _property_interpretations(mol, props):
    """Build readable, cautious explanations for calculated molecular descriptors."""
    from rdkit.Chem import rdMolDescriptors as rdmd
    formula = rdmd.CalcMolFormula(mol)
    exact_mass = Descriptors.ExactMolWt(mol)
    logp = Descriptors.MolLogP(mol)
    tpsa = rdmd.CalcTPSA(mol)
    hbd = rdmd.CalcNumHBD(mol)
    hba = rdmd.CalcNumHBA(mol)
    rot = rdmd.CalcNumRotatableBonds(mol)
    rings = rdmd.CalcNumRings(mol)
    aromatic = rdmd.CalcNumAromaticRings(mol)
    aliphatic = rdmd.CalcNumAliphaticRings(mol)
    heavy = mol.GetNumHeavyAtoms()
    heavy_bonds = sum(1 for b in mol.GetBonds() if b.GetBeginAtom().GetAtomicNum() > 1 and b.GetEndAtom().GetAtomicNum() > 1)
    mol_h = Chem.AddHs(Chem.Mol(mol))
    total_atoms = mol_h.GetNumAtoms()
    total_bonds = mol_h.GetNumBonds()
    pi_bonds = sum(2 if b.GetBondType() == Chem.BondType.TRIPLE else 1 if b.GetBondType() == Chem.BondType.DOUBLE else 0 for b in mol.GetBonds())
    sigma_bonds = total_bonds - pi_bonds
    frac_csp3 = rdmd.CalcFractionCSP3(mol)
    mol_refr = Descriptors.MolMR(mol)
    rows = [
        ("Molecular Formula", formula, "Shows the elemental composition of the molecule and is useful for understanding its basic chemical composition."),
        ("Molecular Weight (g/mol)", f"{Descriptors.MolWt(mol):.3f}", "Represents the average molecular mass. Molecular size and mass can influence diffusion and many physical behaviours."),
        ("Exact Molecular Weight (Da)", f"{exact_mass:.4f}", "Monoisotopic mass calculated from exact isotopic masses; useful in mass-spectrometry comparisons."),
        ("LogP", f"{logp:.3f}", "A calculated estimate of partitioning between oily and aqueous environments. Higher values generally indicate greater lipophilicity; it is not a direct solubility measurement."),
        ("TPSA (Å²)", f"{tpsa:.2f}", "Estimated polar surface area from polar atoms. It helps compare molecular polarity but is not a direct measurement of absorption."),
        ("Hydrogen Bond Donors", str(hbd), f"The structure contains {hbd} hydrogen-bond donor site(s), which can donate hydrogen in suitable intermolecular interactions."),
        ("Hydrogen Bond Acceptors", str(hba), f"The structure contains {hba} hydrogen-bond acceptor site(s), which can accept hydrogen in suitable interactions."),
        ("Rotatable Bonds", str(rot), f"{rot} rotatable bond(s) provide a rough indication of conformational flexibility."),
        ("Ring Count", str(rings), f"The molecule contains {rings} ring(s), which contribute to its overall shape and structural constraints."),
        ("Aromatic Rings", str(aromatic), f"{aromatic} aromatic ring(s) contribute to aromaticity, planarity and delocalized π-electron character where present."),
        ("Aliphatic Rings", str(aliphatic), f"{aliphatic} aliphatic ring(s) contribute cyclic structure without aromatic π-electron delocalization."),
        ("Number of Heavy Atoms", str(heavy), f"Counts non-hydrogen atoms; a useful simple measure of molecular structure size. The hydrogen-expanded model has {total_atoms} total atoms."),
        ("Heavy Atom Bonds", str(heavy_bonds), "Counts connections between non-hydrogen atoms in the molecular graph."),
        ("Sigma (σ) Bonds", str(sigma_bonds), "Approximate count of σ-bond components in the hydrogen-expanded structure; each ordinary bonded pair contributes a σ bond."),
        ("Pi (π) Bonds", str(pi_bonds), "Counts π-bond components associated with double and triple bonds; aromatic bonding is represented by the toolkit's aromatic model."),
        ("Total Bonds (σ + π)", str(total_bonds), "Total bonds in the hydrogen-expanded molecular model. Bond-order conventions can affect simple σ/π summaries."),
        ("Fraction Csp³", f"{frac_csp3:.3f}", "Fraction of carbon atoms assigned sp³ hybridization; a higher value often indicates a more saturated 3D carbon framework."),
        ("Molar Refractivity", f"{mol_refr:.2f}", "A calculated descriptor related to molecular volume and electronic polarizability."),
    ]
    return rows


def _molecule_images(mol):
    """Return in-memory 2D and 3D PNGs for the PDF; 3D coordinates are computational."""
    from rdkit.Chem import Draw
    # 2D structure
    m2 = Chem.Mol(mol)
    Chem.rdDepictor.Compute2DCoords(m2)
    img2 = Draw.MolToImage(m2, size=(950, 560), kekulize=True, wedgeBonds=True)
    b2 = BytesIO()
    img2.save(b2, format="PNG")
    b2.seek(0)

    # A compact 3D ball-and-stick rendering using generated coordinates.
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D  # noqa: F401
    m3 = Chem.AddHs(Chem.Mol(mol))
    params = AllChem.ETKDGv3()
    params.randomSeed = 17
    status = AllChem.EmbedMolecule(m3, params)
    if status != 0:
        raise ValueError("Could not generate a 3D conformer for this molecule.")
    try:
        AllChem.MMFFOptimizeMolecule(m3, maxIters=250)
    except Exception:
        try:
            AllChem.UFFOptimizeMolecule(m3, maxIters=250)
        except Exception:
            pass
    conf = m3.GetConformer()
    coords = [conf.GetAtomPosition(i) for i in range(m3.GetNumAtoms())]
    xyz = [(p.x, p.y, p.z) for p in coords]
    fig = plt.figure(figsize=(7.4, 4.2), dpi=170)
    ax = fig.add_subplot(111, projection="3d")
    element_colors = {1: "#F5F7FA", 6: "#50545A", 7: "#315DCE", 8: "#D83333", 9: "#3B9B63", 15: "#E68A2E", 16: "#E2C13B", 17: "#38A16A"}
    radii = {1: 18, 6: 45, 7: 48, 8: 48, 9: 48, 15: 55, 16: 55, 17: 50}
    for bond in m3.GetBonds():
        a, b = xyz[bond.GetBeginAtomIdx()], xyz[bond.GetEndAtomIdx()]
        ax.plot([a[0], b[0]], [a[1], b[1]], [a[2], b[2]], color="#777D86", linewidth=2.6, zorder=1)
    for i, atom in enumerate(m3.GetAtoms()):
        x, y, z = xyz[i]
        ax.scatter(x, y, z, s=radii.get(atom.GetAtomicNum(), 42), color=element_colors.get(atom.GetAtomicNum(), "#9B59B6"), edgecolors="#33383E", linewidths=0.55, depthshade=True, zorder=3)
    ax.view_init(elev=22, azim=35)
    ax.set_axis_off()
    try:
        ax.set_box_aspect((1, 1, 0.72))
    except Exception:
        pass
    fig.patch.set_facecolor("white")
    plt.tight_layout(pad=0.2)
    b3 = BytesIO()
    fig.savefig(b3, format="PNG", bbox_inches="tight", facecolor="white", transparent=False)
    plt.close(fig)
    b3.seek(0)
    return b2, b3


def pdf_report(student, reg, course, experiment, mol_name, smiles, canonical, props, note=""):
    """Create an A4 molecular analysis report styled like a formal university practical."""
    from reportlab.platypus import Image as RLImage, KeepTogether
    buffer = BytesIO()
    navy = colors.HexColor("#173F70")
    blue = colors.HexColor("#315D9B")
    pale_blue = colors.HexColor("#EAF0F8")
    line_blue = colors.HexColor("#C7D3E2")
    text_col = colors.HexColor("#202B38")
    doc = SimpleDocTemplate(
        buffer, pagesize=A4, leftMargin=1.45*cm, rightMargin=1.45*cm,
        topMargin=1.45*cm, bottomMargin=1.45*cm,
        title=f"Final Student Molecular Analysis Report - {mol_name}",
        author="ChemLab Studio", subject="Educational cheminformatics molecular analysis",
    )
    styles = getSampleStyleSheet()
    styles.add(ParagraphStyle(name="ReportTitleCL", parent=styles["Title"], alignment=TA_CENTER,
        textColor=navy, fontName="Helvetica-Bold", fontSize=18, leading=22, spaceAfter=12))
    styles.add(ParagraphStyle(name="ReportSubtitleCL", parent=styles["Heading2"], alignment=TA_CENTER,
        textColor=navy, fontName="Helvetica-Bold", fontSize=12, leading=15, spaceBefore=0, spaceAfter=13))
    styles.add(ParagraphStyle(name="SectionCL", parent=styles["Heading2"], textColor=navy,
        fontName="Helvetica-Bold", fontSize=12, leading=15, spaceBefore=11, spaceAfter=6, keepWithNext=True))
    styles.add(ParagraphStyle(name="BodyCL", parent=styles["BodyText"], textColor=text_col,
        fontName="Helvetica", fontSize=7.8, leading=9.6, spaceAfter=3))
    styles.add(ParagraphStyle(name="SmallCL", parent=styles["BodyText"], textColor=text_col,
        fontName="Helvetica", fontSize=7.4, leading=9.1))
    styles.add(ParagraphStyle(name="TableHeadCL", parent=styles["BodyText"], textColor=colors.white,
        fontName="Helvetica-Bold", fontSize=7.5, leading=9))
    styles.add(ParagraphStyle(name="PropertyCL", parent=styles["BodyText"], textColor=text_col,
        fontName="Helvetica", fontSize=7.2, leading=8.8))
    styles.add(ParagraphStyle(name="InterpretCL", parent=styles["BodyText"], textColor=text_col,
        fontName="Helvetica", fontSize=6.5, leading=7.7))
    safe = lambda value: escape(str(value or "Not provided"))
    mol = parse_molecule(smiles)
    formula = rdMolDescriptors.CalcMolFormula(mol)
    img2, img3 = _molecule_images(mol)
    exp = EXPERIMENTS.get(experiment, EXPERIMENTS["Molecular descriptor analysis"])
    story = [
        Paragraph("CHEMINFORMATICS VIRTUAL LABORATORY", styles["ReportTitleCL"]),
        Paragraph("FINAL STUDENT MOLECULAR ANALYSIS REPORT", styles["ReportSubtitleCL"]),
    ]
    student_rows = [
        [Paragraph("Student Name", styles["SmallCL"]), Paragraph(safe(student), styles["BodyCL"])],
        [Paragraph("Registration Number", styles["SmallCL"]), Paragraph(safe(reg), styles["BodyCL"])],
        [Paragraph("Input SMILES", styles["SmallCL"]), Paragraph(f"<font name='Courier'>{safe(smiles)}</font>", styles["BodyCL"])],
        [Paragraph("Canonical SMILES", styles["SmallCL"]), Paragraph(f"<font name='Courier'>{safe(canonical)}</font>", styles["BodyCL"])],
        [Paragraph("Molecular Formula", styles["SmallCL"]), Paragraph(safe(formula), styles["BodyCL"])],
    ]
    info = Table(student_rows, colWidths=[4.9*cm, 12.0*cm], hAlign="CENTER")
    info.setStyle(TableStyle([
        ("BACKGROUND", (0,0), (0,-1), pale_blue), ("GRID", (0,0), (-1,-1), 0.45, colors.HexColor("#9AA7B6")),
        ("VALIGN", (0,0), (-1,-1), "MIDDLE"), ("LEFTPADDING", (0,0), (-1,-1), 6),
        ("RIGHTPADDING", (0,0), (-1,-1), 6), ("TOPPADDING", (0,0), (-1,-1), 3.5),
        ("BOTTOMPADDING", (0,0), (-1,-1), 3.5),
    ]))
    story.extend([info, Paragraph("1. Molecular Structures", styles["SectionCL"])])
    struct_tbl = Table([
        [Paragraph("<b>2D Molecular Structure</b>", styles["SmallCL"]), Paragraph("<b>3D Molecular Structure (Ball-and-Stick)</b>", styles["SmallCL"])],
        [RLImage(img2, width=7.4*cm, height=4.35*cm, kind="proportional"), RLImage(img3, width=7.4*cm, height=4.35*cm, kind="proportional")],
    ], colWidths=[8.1*cm, 8.1*cm], hAlign="CENTER")
    struct_tbl.setStyle(TableStyle([
        ("ALIGN", (0,0), (-1,-1), "CENTER"), ("VALIGN", (0,0), (-1,-1), "MIDDLE"),
        ("LEFTPADDING", (0,0), (-1,-1), 2), ("RIGHTPADDING", (0,0), (-1,-1), 2),
        ("TOPPADDING", (0,0), (-1,-1), 2), ("BOTTOMPADDING", (0,0), (-1,-1), 2),
    ]))
    story.extend([struct_tbl, Paragraph("Figure 1. Left: 2D molecular connectivity. Right: generated 3D ball-and-stick geometry showing a possible spatial arrangement of atoms.", styles["SmallCL"]),
        Paragraph("2. Complete Molecular Properties and Interpretation", styles["SectionCL"])])
    property_rows = _property_interpretations(mol, props)
    def make_property_table(rows):
        table_data = [[Paragraph("Property", styles["TableHeadCL"]), Paragraph("Value", styles["TableHeadCL"]), Paragraph("Interpretation", styles["TableHeadCL"])]]
        for name, value, interpretation in rows:
            table_data.append([Paragraph(safe(name), styles["PropertyCL"]), Paragraph(safe(value), styles["PropertyCL"]), Paragraph(safe(interpretation), styles["InterpretCL"])])
        table = Table(table_data, colWidths=[4.0*cm, 2.2*cm, 10.7*cm], repeatRows=1, hAlign="CENTER")
        table.setStyle(TableStyle([
            ("BACKGROUND", (0,0), (-1,0), blue), ("GRID", (0,0), (-1,-1), 0.4, colors.HexColor("#9BA7B5")),
            ("ROWBACKGROUNDS", (0,1), (-1,-1), [colors.white, colors.HexColor("#F3F6FA")]),
            ("VALIGN", (0,0), (-1,-1), "TOP"), ("LEFTPADDING", (0,0), (-1,-1), 4),
            ("RIGHTPADDING", (0,0), (-1,-1), 4), ("TOPPADDING", (0,0), (-1,-1), 2.2),
            ("BOTTOMPADDING", (0,0), (-1,-1), 2.2),
        ]))
        return table
    # Split the descriptor table deliberately, matching the supplied two-page report example.
    split_at = 9
    story.append(make_property_table(property_rows[:split_at]))
    story.extend([PageBreak(), Paragraph("2. Complete Molecular Properties and Interpretation (continued)", styles["SectionCL"]), make_property_table(property_rows[split_at:])])
    story.append(Paragraph("3. Point-wise Molecular Interpretation", styles["SectionCL"]))
    logp = Descriptors.MolLogP(mol)
    tpsa = rdMolDescriptors.CalcTPSA(mol)
    hbd = rdMolDescriptors.CalcNumHBD(mol)
    hba = rdMolDescriptors.CalcNumHBA(mol)
    rot = rdMolDescriptors.CalcNumRotatableBonds(mol)
    rings = rdMolDescriptors.CalcNumRings(mol)
    aromatic = rdMolDescriptors.CalcNumAromaticRings(mol)
    sigma_pi = next((r for r in property_rows if r[0] == "Sigma (σ) Bonds"), None)
    pi_row = next((r for r in property_rows if r[0] == "Pi (π) Bonds"), None)
    fraction = rdMolDescriptors.CalcFractionCSP3(mol)
    bullet_points = [
        ("Molecular identity", f"The submitted structure has molecular formula <b>{safe(formula)}</b> and molecular weight <b>{Descriptors.MolWt(mol):.3f} g/mol</b>."),
        ("Lipophilicity", f"Calculated LogP = <b>{logp:.3f}</b>. Higher values generally indicate greater preference for non-polar environments; this estimate is not a direct solubility measurement."),
        ("Molecular polarity", f"TPSA = <b>{tpsa:.2f} Å²</b>. TPSA reflects contributions from polar atoms and functional groups to the estimated polar surface area."),
        ("Hydrogen bonding", f"The structure contains <b>{hbd}</b> hydrogen-bond donor(s) and <b>{hba}</b> acceptor(s), useful for discussing possible intermolecular interactions."),
        ("Flexibility", f"<b>{rot}</b> rotatable bond(s) indicate the approximate conformational freedom captured by this descriptor."),
        ("Cyclic and aromatic character", f"The molecule contains <b>{rings}</b> ring(s), including <b>{aromatic}</b> aromatic ring(s). Aromatic rings are associated with delocalized π-electron systems."),
        ("Bonding pattern", f"The calculated structure contains <b>{sigma_pi[1] if sigma_pi else '—'} σ bonds</b> and <b>{pi_row[1] if pi_row else '—'} π bonds</b> in the simplified hydrogen-expanded bond count."),
        ("Three-dimensional character", f"Fraction Csp³ = <b>{fraction:.3f}</b>. This descriptor indicates the proportion of carbon atoms assigned sp³ hybridization."),
        ("Overall conclusion", "The calculated descriptors should be interpreted together with the 2D connectivity and generated 3D geometry to understand composition, bonding, polarity, flexibility and structural character."),
    ]
    for title, body in bullet_points:
        story.append(Paragraph(f"• <b>{title}:</b> {body}", styles["BodyCL"]))
    if note and note.strip():
        story.append(Paragraph("4. Student Observations", styles["SectionCL"]))
        story.append(Paragraph(safe(note).replace("\n", "<br/>"), styles["BodyCL"]))
    story.extend([
        Paragraph("5. Learning Note", styles["SectionCL"]),
        Paragraph("The numerical descriptors in this report are computationally calculated from the submitted molecular structure. The generated 3D conformer is one plausible model, not necessarily an experimentally observed geometry. Values should not be interpreted as experimental measurements or as proof of safety, medical suitability, or biological activity.", styles["BodyCL"]),
        Paragraph("Practical details", styles["SectionCL"]),
        Paragraph(f"<b>Practical topic:</b> {safe(experiment)} &nbsp;&nbsp; <b>Course:</b> {safe(course)} &nbsp;&nbsp; <b>Report date:</b> {date.today().strftime('%d %B %Y')}", styles["BodyCL"]),
        Paragraph("Prepared using ChemLab Studio · Educational use only · Review your institution's requirements before submission.", styles["SmallCL"]),
    ])
    def page_decor(canvas, document):
        canvas.saveState()
        width, height = A4
        canvas.setFillColor(navy)
        canvas.rect(0, height-0.32*cm, width, 0.32*cm, fill=1, stroke=0)
        canvas.setStrokeColor(line_blue)
        canvas.line(doc.leftMargin, 1.0*cm, width-doc.rightMargin, 1.0*cm)
        canvas.setFont("Helvetica", 7.2)
        canvas.setFillColor(colors.HexColor("#4B5A6B"))
        canvas.drawString(doc.leftMargin, 0.68*cm, "CHEMINFORMATICS VIRTUAL LABORATORY  •  STUDENT MOLECULAR ANALYSIS")
        canvas.drawRightString(width-doc.rightMargin, 0.68*cm, f"Page {document.page}")
        canvas.restoreState()
    doc.build(story, onFirstPage=page_decor, onLaterPages=page_decor)
    buffer.seek(0)
    return buffer.getvalue()

def openai_answer(history, prompt):
    """Call the OpenAI Responses API and return actionable, key-safe diagnostics."""
    try:
        api_key = st.secrets.get("OPENAI_API_KEY", os.getenv("OPENAI_API_KEY", ""))
    except Exception:
        api_key = os.getenv("OPENAI_API_KEY", "")
    if not api_key or api_key == "YOUR_ACTUAL_OPENAI_API_KEY":
        return None, "OpenAI API key is missing. Add OPENAI_API_KEY in Streamlit → Manage app → Settings → Secrets."
    try:
        from openai import OpenAI, AuthenticationError, PermissionDeniedError, RateLimitError, APIConnectionError, APIStatusError
        client = OpenAI(api_key=api_key, timeout=45.0, max_retries=1)
        conversation = []
        for item in history[-8:]:
            role = "assistant" if item.get("role") == "assistant" else "user"
            conversation.append({"role": role, "content": str(item.get("content", ""))})
        conversation.append({"role": "user", "content": str(prompt)})
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
            return None, "OpenAI returned no text for this request. Please try again with a shorter question."
        return answer, None
    except Exception as e:
        # Exception messages do not contain the secret key; never print the key itself.
        err_name = type(e).__name__
        msg = str(e).replace(api_key, "[API KEY REDACTED]")[:260]
        lower = msg.lower()
        if err_name == "AuthenticationError" or "invalid_api_key" in lower or "incorrect api key" in lower:
            return None, "OpenAI rejected the API key. Create a new key at platform.openai.com/api-keys and update Streamlit Secrets."
        if err_name == "PermissionDeniedError" or "model_not_found" in lower or "do not have access" in lower:
            return None, "Your API project may not have access to gpt-4.1-mini. Check model access in the OpenAI API Platform; the app's model can be changed to one enabled for your project."
        if err_name == "RateLimitError" or "insufficient_quota" in lower or "billing_hard_limit" in lower:
            return None, "OpenAI API quota/billing issue. Check API billing, available credits, and project usage limits in the OpenAI API Platform. ChatGPT subscription billing is separate."
        if err_name == "APIConnectionError":
            return None, "The app could not connect to OpenAI. Try again in a minute and check OpenAI service status."
        if err_name == "APIStatusError":
            status = getattr(e, "status_code", "unknown")
            return None, f"OpenAI returned HTTP {status}: {msg}. Check API billing, model access, and project limits."
        return None, f"OpenAI API error ({err_name}): {msg}"

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
            st.dataframe(pd.DataFrame([{"Descriptor": str(k), "Value": str(v)} for k,v in props.items()]), use_container_width=True, hide_index=True)
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
        st.dataframe(pd.DataFrame([{"Property":str(k),"Value":str(v)} for k,v in props.items()]), use_container_width=True, hide_index=True)
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
    st.dataframe(df.astype(str), use_container_width=True, hide_index=True)
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
            st.dataframe(pd.DataFrame([{"Descriptor":str(k),"Value":str(v)} for k,v in props.items()]), use_container_width=True, hide_index=True)
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
        st.success("OpenAI API key detected. Responses will be generated by GPT-4.1 mini through the OpenAI API.")
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
            st.warning(fallback)
            st.caption(f"Technical detail: {error or 'The AI service did not return an answer.'}")
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
