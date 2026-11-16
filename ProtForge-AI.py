import os
import io
import requests
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Line3DCollection
from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from reportlab.lib.pagesizes import letter
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle, Image
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib import colors

# تلاش برای ایمپورت RDKit (برای تولید SMILES)
try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False
    print("[!] RDKit not found. SMILES files will not be generated.")
    print("    Install via: pip install rdkit")

# ---------------- CONSTANTS & SCALES ---------------- #
KD_SCALE = {
    'A': 1.8, 'R': -4.5, 'N': -3.5, 'D': -3.5, 'C': 2.5, 'Q': -3.5, 'E': -3.5,
    'G': -0.4, 'H': -3.2, 'I': 4.5, 'L': 3.8, 'K': -3.9, 'M': 1.9, 'F': 2.8,
    'P': -1.6, 'S': -0.8, 'T': -0.7, 'W': -0.9, 'Y': -1.3, 'V': 4.2
}
AA_LOG_P = {'A': 0.31, 'R': -1.01, 'N': -0.60, 'D': -0.77, 'C': 1.54, 'Q': -0.22, 'E': -0.64, 'G': 0.0, 'H': -0.40, 'I': 1.80, 'L': 1.70, 'K': -0.99, 'M': 1.23, 'F': 1.79, 'P': 0.72, 'S': -0.04, 'T': 0.26, 'W': 2.25, 'Y': 0.96, 'V': 1.22}
AA_MAP_3TO1 = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E','GLY':'G','HIS':'H','ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S','THR':'T','TRP':'W','TYR':'Y','VAL':'V'}

# مقیاس‌های جدید برای معیارهای Developability
# Boman index (پتانسیل برهم‌کنش پروتئین-پروتئین)
BOMAN_SCALE = {'A': 0.17, 'R': 1.87, 'N': 0.42, 'D': -1.05, 'C': 0.09, 'Q': 0.11, 'E': -0.79, 'G': -0.03, 'H': 0.60, 'I': -1.27, 'L': -1.27, 'K': 1.87, 'M': -1.09, 'F': -1.43, 'P': 0.26, 'S': 0.30, 'T': 0.27, 'W': -1.43, 'Y': -0.79, 'V': -1.27}

# ---------------- HELPER FUNCTIONS ---------------- #

def fetch_esmfold_structure(sequence):
    url = "https://api.esmatlas.com/foldSequence/v1/pdb/"
    headers = {"Content-Type": "application/x-www-form-urlencoded"}
    try:
        print(f"  [*] Predicting 3D structure via ESMFold API...")
        response = requests.post(url, headers=headers, data=sequence, timeout=60)
        if response.status_code == 200 and "ATOM" in response.text:
            return response.text
        else:
            return None
    except Exception as e:
        print(f"  [!] Connection error: {e}")
        return None

def parse_pdb_ss(pdb_data):
    ss_map = {}
    for line in pdb_data.split('\n'):
        if line.startswith("HELIX"):
            start, end = int(line[21:25]), int(line[33:37])
            for i in range(start, end + 1): ss_map[i] = "H"
        elif line.startswith("SHEET"):
            start, end = int(line[22:26]), int(line[33:37])
            for i in range(start, end + 1): ss_map[i] = "E"
    return ss_map

def plot_3d_ribbon(pdb_data, title="3D Structure"):
    atoms = {'CA': []}
    residues = []
    ss_map = parse_pdb_ss(pdb_data)
    use_ss_colors = len(ss_map) > 0
    
    for line in pdb_data.split('\n'):
        if line.startswith("ATOM") and line[12:16].strip() == "CA":
            res_seq = int(line[22:26])
            res_name = line[17:20].strip()
            atoms['CA'].append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
            residues.append((res_seq, AA_MAP_3TO1.get(res_name, 'A')))
    
    if len(atoms['CA']) < 2: return None

    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111, projection='3d')
    
    lines = []
    colors_list = []
    for i in range(len(atoms['CA']) - 1):
        lines.append([atoms['CA'][i], atoms['CA'][i+1]])
        if use_ss_colors:
            res_num = residues[i][0]
            ss = ss_map.get(res_num, "L")
            if ss == "H": colors_list.append('#e74c3c')
            elif ss == "E": colors_list.append('#f1c40f')
            else: colors_list.append('#2ecc71')
        else:
            aa = residues[i][1]
            kd_val = KD_SCALE.get(aa, 0)
            if kd_val > 2: colors_list.append('#e67e22')
            elif kd_val < -1: colors_list.append('#3498db')
            else: colors_list.append('#95a5a6')

    lc = Line3DCollection(lines, colors=colors_list, linewidths=7, alpha=0.9)
    ax.add_collection3d(lc)
    
    xs = [p[0] for p in atoms['CA']]
    ys = [p[1] for p in atoms['CA']]
    zs = [p[2] for p in atoms['CA']]
    ax.scatter(xs, ys, zs, c='white', s=15, alpha=0.6, edgecolors='black')
    
    max_range = np.array([max(xs)-min(xs), max(ys)-min(ys), max(zs)-min(zs)]).max() / 2.0
    mid_x, mid_y, mid_z = (max(xs)+min(xs))*0.5, (max(ys)+min(ys))*0.5, (max(zs)+min(zs))*0.5
    ax.set_xlim(mid_x - max_range, mid_x + max_range)
    ax.set_ylim(mid_y - max_range, mid_y + max_range)
    ax.set_zlim(mid_z - max_range, mid_z + max_range)
    
    ax.set_title(title, fontsize=14, fontweight='bold')
    ax.axis('off')
    
    from matplotlib.lines import Line2D
    if use_ss_colors:
        legend_elements = [Line2D([0], [0], color='#e74c3c', lw=4, label='Alpha Helix'),
                           Line2D([0], [0], color='#f1c40f', lw=4, label='Beta Sheet'),
                           Line2D([0], [0], color='#2ecc71', lw=4, label='Loop')]
    else:
        legend_elements = [Line2D([0], [0], color='#e67e22', lw=4, label='Hydrophobic'),
                           Line2D([0], [0], color='#95a5a6', lw=4, label='Neutral'),
                           Line2D([0], [0], color='#3498db', lw=4, label='Hydrophilic')]
    ax.legend(handles=legend_elements, loc='lower right', fontsize=10)
    
    img_buf = io.BytesIO()
    plt.savefig(img_buf, format='png', bbox_inches='tight', dpi=200)
    img_buf.seek(0)
    plt.close()
    return img_buf

def plot_hydrophobicity(sequence):
    window_size = 9
    scores = [KD_SCALE.get(aa, 0) for aa in sequence]
    if len(scores) < window_size: return None
    smoothed = np.convolve(scores, np.ones(window_size)/window_size, mode='valid')
    
    fig, ax = plt.subplots(figsize=(8, 3))
    ax.plot(range(window_size, len(scores) + 1), smoothed, color='darkblue')
    ax.axhline(0, color='red', linestyle='--', linewidth=1)
    ax.set_title("Hydrophobicity Profile (Kyte-Doolittle)")
    ax.grid(True, alpha=0.3)
    
    img_buf = io.BytesIO()
    plt.savefig(img_buf, format='png', bbox_inches='tight', dpi=150)
    img_buf.seek(0)
    plt.close()
    return img_buf

def find_target_regions(sequence):
    target_regions = []
    window = 10
    for i in range(len(sequence) - window + 1):
        pep = sequence[i:i+window]
        flex_hydro = sum(1 for aa in pep if aa in 'GSQERKDN')
        if flex_hydro >= 7:
            target_regions.append({
                "position": f"{i+1}-{i+window}",
                "sequence": pep,
                "score": flex_hydro,
                "reason": "High flexibility & surface exposure"
            })
    target_regions.sort(key=lambda x: x["score"], reverse=True)
    final_targets = []
    occupied = set()
    for t in target_regions:
        start, end = map(int, t["position"].split('-'))
        if not any(pos in occupied for pos in range(start, end + 1)):
            final_targets.append(t)
            occupied.update(range(start, end + 1))
        if len(final_targets) >= 3: break
    return final_targets

# ---------------- NEW DEVELOPABILITY METRICS ---------------- #

def calculate_aliphatic_index(seq):
    """شاخص آلیفاتیک برای پایداری حرارتی."""
    a = seq.count('A') / len(seq) * 100
    v = seq.count('V') / len(seq) * 100
    i = seq.count('I') / len(seq) * 100
    l = seq.count('L') / len(seq) * 100
    return a + 2.9 * v + 3.9 * (i + l)

def calculate_boman_index(seq):
    """شاخص بومن برای پتانسیل برهم‌کنش پروتئین-پروتئین."""
    return sum(BOMAN_SCALE.get(aa, 0) for aa in seq) / len(seq)

def count_protease_sites(seq, enzyme='trypsin'):
    """شمارش سایت‌های برش پروتئاز."""
    sites = 0
    if enzyme == 'trypsin':
        for i in range(len(seq) - 1):
            if seq[i] in 'KR' and seq[i+1] != 'P':
                sites += 1
    elif enzyme == 'chymotrypsin':
        for i in range(len(seq) - 1):
            if seq[i] in 'FWY' and seq[i+1] != 'P':
                sites += 1
    return sites

def calculate_solubility_score(seq):
    """امتیاز حلالیت بر اساس بار و آمینواسیدهای قطبی (0-1)."""
    polar = sum(seq.count(aa) for aa in 'RKDENQST') / len(seq)
    charge = abs(sum(1 for aa in seq if aa in 'KR') - sum(1 for aa in seq if aa in 'DE'))
    # امتیاز بالاتر = حلالیت بهتر
    score = (polar * 0.7) + (min(charge, 5) / 5 * 0.3)
    return min(score, 1.0)

def calculate_aggregation_score(seq):
    """امتیاز تجمع (0-1) بر اساس پنجره‌های آبگریز."""
    window = 5
    max_hydro = 0
    for i in range(len(seq) - window + 1):
        window_seq = seq[i:i+window]
        hydro = sum(KD_SCALE.get(aa, 0) for aa in window_seq) / window
        max_hydro = max(max_hydro, hydro)
    # نرمال‌سازی (حداکثر آبگریزی ممکن ~4.5)
    return min(max(max_hydro / 4.5, 0), 1.0)

def predict_mhc_ic50_simple(seq):
    """پیش‌بینی ساده IC50 بر اساس موتیف‌های اتصال MHC-I (9-mer)."""
    # این یک پیش‌بینی بسیار ساده است. برای دقت واقعی از API استفاده کنید.
    strong_binders = 0
    for i in range(len(seq) - 8):
        peptide = seq[i:i+9]
        # موتیف ساده: آمینواسید آبگریز در موقعیت 2 و 9
        if peptide[1] in 'AVILMF' and peptide[8] in 'AVILMF':
            strong_binders += 1
    if strong_binders > 0:
        return f"< 500 (Strong: {strong_binders} motifs)"
    else:
        return "> 5000 (Weak)"

def predict_tm_simple(seq):
    """پیش‌بینی ساده Tm بر اساس محتوای مارپیچ و آبگریزی."""
    # تخمین بسیار خام: هرچه آبگریزی و مارپیچ بیشتر، Tm بالاتر
    gravy = sum(KD_SCALE.get(aa, 0) for aa in seq) / len(seq)
    # فرض: Tm پایه 40، با آبگریزی افزایش می‌یابد
    tm = 40 + (gravy * 10) + (len(seq) * 0.05)
    return max(30, min(90, tm))

# ---------------- MAIN PDF GENERATOR ---------------- #

def generate_pdf(seq_record, output_dir="reports"):
    raw_id = seq_record.id
    safe_id = raw_id.split('|')[1] if '|' in raw_id else raw_id
    safe_id = safe_id.replace('/', '_').replace('\\', '_').replace(' ', '_')
    seq_str = str(seq_record.seq).upper()
    seq_len = len(seq_str)
    
    if any(aa not in KD_SCALE for aa in seq_str):
        print(f"  [!] Skipping {safe_id}: Non-standard amino acids.")
        return

    print(f"\n[>] Processing: {safe_id} (Length: {seq_len})")
    analysis = ProteinAnalysis(seq_str)
    
    # --- تولید فایل SMILES ---
    smiles_file = os.path.join(output_dir, f"{safe_id}.smi")
    if RDKIT_AVAILABLE:
        try:
            mol = Chem.MolFromSequence(seq_str)
            if mol:
                smiles_str = Chem.MolToSmiles(mol, isomericSmiles=True)
                with open(smiles_file, 'w') as f:
                    f.write(f"{smiles_str}\t{safe_id}\n")
                print(f"  [+] SMILES saved: {smiles_file}")
            else:
                print(f"  [!] Could not generate SMILES for {safe_id}")
        except Exception as e:
            print(f"  [!] RDKit error for {safe_id}: {e}")
    else:
        print(f"  [!] RDKit not available. Skipping SMILES for {safe_id}.")

    # 1. Electrochemical
    elec_props = {
        "Molecular Weight (Da)": f"{analysis.molecular_weight():.2f}",
        "Isoelectric Point (pI)": f"{analysis.isoelectric_point():.2f}",
        "Instability Index": f"{analysis.instability_index():.2f}",
        "GRAVY (Hydrophobicity)": f"{analysis.gravy():.3f}",
        "Net Charge (pH 7)": f"{sum(1 for aa in seq_str if aa in 'KR') - sum(1 for aa in seq_str if aa in 'DE'):.2f}",
        "Molar Extinction (280nm)": f"{analysis.molar_extinction_coefficient()[0]}",
    }

    # 2. ADMET
    logP_approx = sum(AA_LOG_P.get(aa, 0) for aa in seq_str) / seq_len
    sol_proxy = sum(seq_str.count(aa) for aa in 'RKDENQST') / seq_len
    admet_props = {
        "LogP (Predicted)": f"{logP_approx:.2f}",
        "Aqueous Solubility Proxy": "High" if sol_proxy > 0.4 else ("Moderate" if sol_proxy > 0.2 else "Low"),
        "Blood-Brain Barrier": "Likely" if logP_approx > 1.5 and sol_proxy < 0.2 else "Unlikely",
        "H-Bond Donors": f"{sum(seq_str.count(aa) for aa in 'STNQKRH')}",
        "H-Bond Acceptors": f"{sum(seq_str.count(aa) for aa in 'STNQDE')}",
        "CYP450 Inhibition Risk": "High" if logP_approx > 3.0 else "Low",
    }

    # 3. Target Regions
    target_regions = find_target_regions(seq_str)
    flex_aa = sum(seq_str.count(aa) for aa in ['G', 'S', 'P', 'T', 'Q', 'E', 'K', 'R'])
    binding_props = {
        "Flexibility Index": f"{(flex_aa / seq_len) * 100:.1f}%",
        "Hydrophobic Ratio": f"{(sum(seq_str.count(aa) for aa in ['A', 'V', 'I', 'L', 'M', 'F', 'W', 'Y']) / seq_len) * 100:.1f}%",
        "Charged Residues": f"{(sum(seq_str.count(aa) for aa in ['D', 'E', 'K', 'R', 'H']) / seq_len) * 100:.1f}%",
        "Binding Potential": "High" if flex_aa > seq_len*0.3 else "Moderate",
    }

    # 4. Developability Metrics (جدید)
    aliphatic_index = calculate_aliphatic_index(seq_str)
    boman_index = calculate_boman_index(seq_str)
    trypsin_sites = count_protease_sites(seq_str, 'trypsin')
    chymo_sites = count_protease_sites(seq_str, 'chymotrypsin')
    solubility_score = calculate_solubility_score(seq_str)
    aggregation_score = calculate_aggregation_score(seq_str)
    mhc_ic50 = predict_mhc_ic50_simple(seq_str)
    tm_pred = predict_tm_simple(seq_str)
    
    # Liability sites (جستجوی موتیف‌های ساده)
    liabilities = []
    if "NG" in seq_str: liabilities.append("Deamidation (NG)")
    if "DG" in seq_str: liabilities.append("Isomerization (DG)")
    if "DP" in seq_str: liabilities.append("Acid-labile (DP)")
    liability_str = ", ".join(liabilities) if liabilities else "None detected"

    developability_metrics = {
        "Predicted MHC-I IC50": mhc_ic50,
        "Predicted Tm (°C)": f"{tm_pred:.1f}",
        "Aggregation Propensity Score": f"{aggregation_score:.2f} (0-1)",
        "Liability Sites": liability_str,
        "Aliphatic Index (Stability)": f"{aliphatic_index:.1f}",
        "Boman Index (Interaction)": f"{boman_index:.2f}",
        "Solubility Score (0-1)": f"{solubility_score:.2f}",
        "Protease Sites (Trypsin)": f"{trypsin_sites}",
        "Protease Sites (Chymotrypsin)": f"{chymo_sites}",
    }

    # 5. Fetch 3D
    pdb_data = fetch_esmfold_structure(seq_str) if seq_len <= 400 else None
    
    # SS Percentages safely
    if pdb_data:
        ss_map = parse_pdb_ss(pdb_data)
        if len(ss_map) > 0:
            h_count = sum(1 for v in ss_map.values() if v == 'H')
            e_count = sum(1 for v in ss_map.values() if v == 'E')
            l_count = len(ss_map) - h_count - e_count
            ss_props = {
                "Alpha Helix (H)": f"{h_count} residues ({(h_count/len(ss_map))*100:.1f}%)",
                "Beta Sheet (E)": f"{e_count} residues ({(e_count/len(ss_map))*100:.1f}%)",
                "Loop / Coil (L)": f"{l_count} residues ({(l_count/len(ss_map))*100:.1f}%)",
            }
        else:
            helix, turn, sheet = analysis.secondary_structure_fraction()
            ss_props = {
                "Alpha Helix (Predicted)": f"{helix*100:.1f}%",
                "Beta Sheet (Predicted)": f"{sheet*100:.1f}%",
                "Turn / Coil (Predicted)": f"{turn*100:.1f}%",
                "Note": "3D annotations missing. Using sequence-based prediction."
            }
    else:
        ss_props = {"Note": "3D structure unavailable."}

    # ================= BUILD PDF ================= #
    if not os.path.exists(output_dir): os.makedirs(output_dir)
    pdf_filename = os.path.join(output_dir, f"{safe_id}_report.pdf")
    doc = SimpleDocTemplate(pdf_filename, pagesize=letter)
    styles = getSampleStyleSheet()
    elements = []

    title_style = ParagraphStyle(name='TitleStyle', parent=styles['Heading1'], textColor=colors.HexColor("#0f766e"))
    elements.append(Paragraph(f"ProtForge-AI Report: {safe_id}", title_style))
    elements.append(Spacer(1, 10))

    def add_table(title, data, header_color, col_widths=[180, 300]):
        elements.append(Paragraph(f"<b>{title}</b>", styles['Heading2']))
        t = Table(data, colWidths=col_widths)
        t.setStyle(TableStyle([
            ('BACKGROUND', (0, 0), (-1, 0), header_color),
            ('ALIGN', (0, 0), (-1, -1), 'LEFT'),
            ('GRID', (0, 0), (-1, -1), 0.5, colors.grey),
            ('FONTNAME', (0, 0), (-1, -1), 'Helvetica'),
            ('FONTSIZE', (0, 0), (-1, -1), 9),
            ('VALIGN', (0, 0), (-1, -1), 'TOP'),
        ]))
        elements.append(t)
        elements.append(Spacer(1, 12))

    add_table("1. Electrochemical & Physicochemical", [["Property", "Value"]] + list(elec_props.items()), colors.HexColor("#ccfbf1"))
    add_table("2. ADMET & Drug-Likeness", [["Property", "Value"]] + list(admet_props.items()), colors.HexColor("#fef3c7"))
    
    tr_data = [["Region (Position)", "Amino Acid Sequence", "Design Strategy"]]
    if target_regions:
        for idx, tr in enumerate(target_regions):
            tr_data.append([f"Site {idx+1} ({tr['position']})", tr['sequence'], tr['reason']])
    else:
        tr_data.append(["N/A", "N/A", "No highly flexible surface regions found."])
    
    elements.append(Paragraph("<b>3. Target Regions for Antibody / Peptide Design</b>", styles['Heading2']))
    elements.append(Paragraph("These are highly flexible and hydrophilic surface loops ideal for antibody binding or peptide mimetics.", styles['Normal']))
    elements.append(Spacer(1, 5))
    t_tr = Table(tr_data, colWidths=[100, 150, 230])
    t_tr.setStyle(TableStyle([
        ('BACKGROUND', (0, 0), (-1, 0), colors.HexColor("#e0e7ff")),
        ('ALIGN', (0, 0), (-1, -1), 'LEFT'),
        ('GRID', (0, 0), (-1, -1), 0.5, colors.grey),
        ('FONTNAME', (0, 0), (-1, -1), 'Helvetica'),
        ('FONTSIZE', (0, 0), (-1, -1), 9),
    ]))
    elements.append(t_tr)
    elements.append(Spacer(1, 12))

    add_table("4. Secondary Structure Composition", [["Type", "Details"]] + list(ss_props.items()), colors.HexColor("#fce7f3"))
    add_table("5. Antibody Binding Features", [["Feature", "Assessment"]] + list(binding_props.items()), colors.HexColor("#e0e7ff"))
    
    # جدول جدید Developability
    add_table("6. Drug Developability & Safety Metrics", [["Metric", "Value"]] + list(developability_metrics.items()), colors.HexColor("#fde68a"))

    if pdb_data:
        img_3d = plot_3d_ribbon(pdb_data, title="Predicted 3D Structure (Ribbon-like)")
        if img_3d:
            elements.append(Paragraph("<b>7. Predicted 3D Structure</b>", styles['Heading2']))
            elements.append(Image(img_3d, width=450, height=350))
            elements.append(Spacer(1, 15))

    img_hydro = plot_hydrophobicity(seq_str)
    if img_hydro:
        elements.append(Paragraph("<b>8. Hydrophobicity Profile</b>", styles['Heading2']))
        elements.append(Image(img_hydro, width=500, height=200))
    
    doc.build(elements)
    print(f"  [+] PDF saved: {pdf_filename}")

def main(fasta_file):
    if not os.path.exists(fasta_file):
        print(f"Error: {fasta_file} not found.")
        return
    for record in SeqIO.parse(fasta_file, "fasta"):
        generate_pdf(record)

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", required=True, help="Input FASTA file")
    args = parser.parse_args()
    main(args.input)
