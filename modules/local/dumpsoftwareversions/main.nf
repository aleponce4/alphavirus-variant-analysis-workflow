process DUMP_SOFTWARE_VERSIONS {
    label 'process_low'

    container "${params.container_multiqc}"

    input:
    path versions_file

    output:
    path "software_versions.yml"    , emit: yml
    path "software_versions_mqc.yml", emit: mqc_yml
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    cat << 'EOF' > parse_versions.py
import sys, yaml

versions = {}
with open(sys.argv[1], 'r') as f:
    text = f.read()

current_proc = None
for line in text.splitlines():
    if not line.strip():
        continue
    if not line.startswith(' ') and ':' in line:
        proc_key = line.split(':')[0].strip().strip('"\\'')
        current_proc = proc_key
        if current_proc not in versions:
            versions[current_proc] = {}
    elif line.startswith(' ') and ':' in line and current_proc:
        parts = line.split(':', 1)
        tool_k = parts[0].strip()
        tool_v = parts[1].strip()
        versions[current_proc][tool_k] = tool_v

with open('software_versions.yml', 'w') as f:
    yaml.dump(versions, f, default_flow_style=False)

table_html = '<table class="table"><thead><tr><th>Process</th><th>Tool</th><th>Version</th></tr></thead><tbody>'
rows = ''
for proc, tools in sorted(versions.items()):
    if isinstance(tools, dict):
        for tool, ver in sorted(tools.items()):
            rows += f'<tr><td>{proc}</td><td>{tool}</td><td>{ver}</td></tr>'

mqc_data = {
    'id': 'software_versions',
    'section_name': 'Software Versions',
    'plot_type': 'html',
    'description': 'Software versions collected during pipeline execution.',
    'data': table_html + rows + '</tbody></table>'
}

with open('software_versions_mqc.yml', 'w') as f:
    yaml.dump(mqc_data, f, default_flow_style=False)
EOF

    python3 parse_versions.py ${versions_file}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
        pyyaml: \$(python3 -c "import yaml; print(yaml.__version__)")
    END_VERSIONS
    """

    stub:
    """
    touch software_versions.yml
    touch software_versions_mqc.yml

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
        pyyaml: \$(python3 -c "import yaml; print(yaml.__version__)")
    END_VERSIONS
    """
}
