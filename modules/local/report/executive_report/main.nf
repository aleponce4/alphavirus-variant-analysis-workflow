process EXECUTIVE_REPORT {
    label 'process_low'

    container "${params.container_python_reporting}"

    input:
    path summary_tsv
    path cov_plots
    path var_plots
    path sel_table
    path hap_plots
    path software_versions
    path qmd_template

    output:
    path "executive_report.html", emit: html
    path "versions.yml"         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    cat << 'EOF' > render_report.py
import sys, os, glob, base64

dataset = sys.argv[1]

def embed_image(path):
    if os.path.exists(path) and os.path.getsize(path) > 0:
        with open(path, "rb") as f:
            b64 = base64.b64encode(f.read()).decode('utf-8')
        return f'<img src="data:image/png;base64,{b64}" style="max-width:100%; height:auto; border:1px solid #e2e8f0; border-radius:8px; margin:16px 0; box-shadow:0 2px 4px rgba(0,0,0,0.05);"/>'
    return '<p style="color:#718096; font-style:italic;">Plot not available for this run</p>'

def tsv_to_html(filepath):
    if not os.path.exists(filepath) or os.path.getsize(filepath) == 0:
        return "<p>Data not available.</p>"
    try:
        with open(filepath, "r") as f:
            lines = [line.strip().split("\t") for line in f if line.strip()]
        if not lines:
            return "<p>Data empty.</p>"
        headers = lines[0]
        rows = lines[1:]
        th = "".join(f"<th>{h}</th>" for h in headers)
        tr_list = []
        for row in rows:
            td = "".join(f"<td>{cell}</td>" for cell in row)
            tr_list.append(f"<tr>{td}</tr>")
        tbody = "".join(tr_list)
        return f'<table class="styled-table"><thead><tr>{th}</tr></thead><tbody>{tbody}</tbody></table>'
    except Exception as e:
        return f"<p>Error loading table: {e}</p>"

summary_html = tsv_to_html("consolidated_sample_summary.tsv")
cov_img = embed_image(f"{dataset}_coverage_summary.png")
var_img = embed_image(f"{dataset}_variant_density.png")

sel_html = ""
if os.path.exists("selection_gene_key_table.tsv") and os.path.getsize("selection_gene_key_table.tsv") > 0:
    sel_table_html = tsv_to_html("selection_gene_key_table.tsv")
    sel_html = (
        '<section class="card">'
        '<h2>Natural Selection Analysis (SNPGenie)</h2>'
        '<p>Population diversity metrics (&pi;, &pi;N, &pi;S) and gene-level selection statistics.</p>'
        + sel_table_html +
        '</section>'
    )

hap_img = ""
hap_plots = glob.glob("*haplotype*.png")
if hap_plots:
    hap_img = (
        '<section class="card">'
        '<h2>Haplotype Reconstruction</h2>'
        '<p>Quasispecies haplotype reconstruction outputs.</p>'
        + embed_image(hap_plots[0]) +
        '</section>'
    )

ver_html = "<p>No version information available.</p>"
if os.path.exists("software_versions.yml"):
    try:
        with open("software_versions.yml", "r") as f:
            lines = f.readlines()
        rows = ""
        current_proc = ""
        for line in lines:
            line_str = line.strip()
            if not line_str or line_str.startswith("#"):
                continue
            if line.startswith("  ") or line.startswith("\t"):
                if ":" in line_str:
                    parts = line_str.split(":", 1)
                    tool = parts[0].strip()
                    ver = parts[1].strip().replace('"', '').replace("'", "")
                    rows += f"<tr><td><strong>{current_proc}</strong></td><td>{tool}</td><td><code>{ver}</code></td></tr>"
            elif line_str.endswith(":"):
                current_proc = line_str[:-1].strip().replace('"', '').replace("'", "")
        if rows:
            ver_html = f'<table class="styled-table"><thead><tr><th>Process</th><th>Tool</th><th>Version</th></tr></thead><tbody>{rows}</tbody></table>'
    except Exception as e:
        ver_html = f"<p>Error loading versions: {e}</p>"

html_content = (
    '<!DOCTYPE html><html lang="en"><head><meta charset="UTF-8">'
    '<meta name="viewport" content="width=device-width, initial-scale=1.0">'
    f'<title>Executive Summary Report - {dataset}</title>'
    '<style>'
    'body { font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, Helvetica, Arial, sans-serif; background-color: #f8fafc; color: #1e293b; line-height: 1.6; margin: 0; padding: 40px 20px; }'
    '.container { max-width: 1100px; margin: 0 auto; background: #ffffff; padding: 40px; border-radius: 12px; box-shadow: 0 4px 6px -1px rgba(0,0,0,0.1), 0 2px 4px -1px rgba(0,0,0,0.06); }'
    'header { border-bottom: 2px solid #e2e8f0; padding-bottom: 20px; margin-bottom: 30px; }'
    'h1 { color: #0f172a; margin: 0 0 8px 0; font-size: 28px; }'
    '.subtitle { color: #64748b; font-size: 16px; margin: 0; }'
    '.card { background: #ffffff; border: 1px solid #e2e8f0; border-radius: 8px; padding: 24px; margin-bottom: 24px; }'
    'h2 { color: #1e293b; font-size: 20px; margin-top: 0; border-bottom: 1px solid #f1f5f9; padding-bottom: 10px; }'
    '.styled-table { width: 100%; border-collapse: collapse; margin: 16px 0; font-size: 14px; }'
    '.styled-table th { background-color: #f1f5f9; color: #334155; text-align: left; padding: 12px; font-weight: 600; border-bottom: 2px solid #cbd5e1; }'
    '.styled-table td { padding: 10px 12px; border-bottom: 1px solid #e2e8f0; }'
    '.styled-table tr:nth-child(even) { background-color: #f8fafc; }'
    'code { background: #f1f5f9; color: #0f172a; padding: 2px 6px; border-radius: 4px; font-size: 13px; font-family: monospace; }'
    '</style></head><body><div class="container"><header>'
    '<h1>Executive Summary Report</h1>'
    f'<p class="subtitle">Dataset: <strong>{dataset}</strong> | Pipeline: Alphavirus Variant Analysis</p>'
    '</header>'
    '<section class="card"><h2>Workflow Overview & Key Metrics</h2>'
    + summary_html +
    '</section><section class="card"><h2>Coverage & Depth Analysis</h2>'
    + cov_img +
    '</section><section class="card"><h2>Variant Discovery & Density</h2>'
    + var_img +
    '</section>'
    + sel_html + hap_img +
    '<section class="card"><h2>Software Provenance</h2>'
    + ver_html +
    '</section></div></body></html>'
)

with open("executive_report.html", "w") as f:
    f.write(html_content)

EOF

    python3 render_report.py ${params.dataset}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    touch executive_report.html

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """
}
