# AGENT.md — Executor Specification: Legacy Shell Workflow → Production Nextflow DSL2

> **Document type:** Binding implementation specification.
> **Audience:** Executor Agent (coding model).
> **Deliverable from Executor:** A complete, tested, release-ready Nextflow DSL2 repository built **phase-by-phase, in order**, exactly as specified below.
> **Legacy source:** `run_full_pipeline.sh`, `Variant_discovery_pipeline/Scripts/*.sh`, `Variant_discovery_pipeline/Scripts/Helpers/*`, `envs/*.yml` in this repository.

---

## 0. How the Executor Must Use This Document

1. Execute the 5 phases **strictly sequentially**. Do not start Phase N+1 until every checkpoint in Phase N passes.
2. Each phase has a **Definition of Done (DoD)** table with shell commands. Run them verbatim, capture output, and only proceed when all exit codes are `0`.
3. If a checkpoint fails twice after attempted fixes, **stop and report** the failing command, its output, and your diagnosis. Do not silently skip or weaken a check.
4. **Never** run `git commit`, `git tag`, `git push`, or any other git mutation without explicit user confirmation at the moment it is needed (this is a hard repository rule and applies especially to the Phase 5 release tag).
5. Keep changes minimal and idiomatic. Follow nf-core community conventions where this document does not explicitly override them.
6. If you deviate from this specification for any reason (e.g., a container tag no longer exists), record the deviation in a `DEVIATIONS.md` file at repo root and continue with the closest compliant alternative.

### 0.1 Legacy workflow inventory (source of truth)

| # | Legacy component | Tool(s) | Conda env | Executor disposition |
|---|---|---|---|---|
| 1 | `Scripts/extract_viral_bams.sh` | `samtools view/sort/index` | `lofreq-env` | → module `extract_viral_bam` |
| 2 | `Scripts/run_lofreq.sh` | `lofreq call-parallel`, `lofreq filter`, `bgzip`, `tabix` | `lofreq-env` | → modules `lofreq/call`, `lofreq/filter` |
| 3 | `Scripts/run_ivar.sh` | `samtools mpileup` + `ivar variants`, `ivar consensus` | `ivar_env` | → modules `ivar/variants`, `ivar/consensus` |
| 4 | `Scripts/annotate_all.sh` + `Scripts/ivar_variants_to_vcf.py` | `bcftools csq --local-csq` | `annotation-env` | → module `bcftools/csq` + `bin/ivar_variants_to_vcf.py` |
| 5 | `Scripts/calculate_coverage.sh` | `samtools depth` + custom summary | `lofreq-env` | → modules `coverage/depth`, `coverage/summarize` |
| 6 | `Scripts/run_snpgenie.sh` + Helpers (`summarize_snpgenie_outputs.py`, `analyze_delta_selection.py`, `analyze_delta_limma.R`, `build_compact_selection_tables.py`) | SNPGenie (Perl), Python, R | local clones/envs | → subworkflow `selection` (module `snpgenie/run` + 4 `bin/` scripts) |
| 7 | `Scripts/run_cliquesnv.sh` | CliqueSNV (Java) | `env_cliquesnv` | → module `cliquesnv` |
| 8 | `Scripts/run_viloca.sh` | ShoRAH/VILOCA | `env_viloca` | → module `viloca` |
| 9 | Helpers: `generate_run_summary.py`, `generate_variant_plots.py`, `generate_coverage_plots.py`, `generate_haplotype_plots.py`, `export_excel_variants.py` | Python (pandas, matplotlib, seaborn, openpyxl) | `annotation-env` | → subworkflow `reporting` (`bin/` scripts, argparse CLIs) |
| 10 | Helpers: `parse_gb_to_gff3.py`, `convert_gff3_to_gtf.py` | Python | `annotation-env` | → `bin/` utility scripts + prep module (kept; run manually, not in main DAG) |
| 11 | WSL→Windows result sync block at end of `run_full_pipeline.sh` | `cp`, WSL paths | — | **DROP.** Replaced by `publishDir`. Absolute WSL/Windows paths are forbidden. |
| 12 | `FORCE_RECALL` / skip-if-exists logic in every shell script | bash | — | **DROP.** Replaced by Nextflow `-resume` caching. |

### 0.2 Legacy configuration → Nextflow parameter mapping

Port these defaults verbatim into `nextflow.config` (`params` scope). Names are snake_case Nextflow params.

| Legacy env var (config.sh) | Nextflow param | Default | Notes |
|---|---|---|---|
| `VIRAL_CONTIG` | `params.viral_contig` | `"KP282671.1"` | Per-dataset override via CLI or sample sheet |
| `REFERENCE` | `params.fasta` | `null` (required) | Viral-only reference FASTA |
| `ANNOTATION` | `params.gff` | `null` (required) | GFF3 with CDS features (validated at input check) |
| `IVAR_MIN_VARIANT_DEPTH` | `params.ivar_min_depth` | `1000` | |
| `IVAR_MIN_VARIANT_FREQ` | `params.ivar_min_freq` | `0.01` | |
| `IVAR_MIN_BASE_QUALITY` | `params.ivar_min_bq` | `30` | |
| `IVAR_MIN_CONSENSUS_COVERAGE` | `params.ivar_consensus_min_cov` | `10` | |
| `IVAR consensus -t 0.5` | `params.ivar_consensus_threshold` | `0.5` | Hard-coded in legacy; promote to param |
| `LOFREQ_MIN_VARIANT_DEPTH` | `params.lofreq_min_depth` | `1000` | |
| `LOFREQ_MIN_VARIANT_FREQ` | `params.lofreq_min_freq` | `0.01` | |
| `LOFREQ_MIN_BASE_QUALITY` | `params.lofreq_min_bq` | `30` | |
| `LOFREQ_MIN_MAP_QUALITY` | `params.lofreq_min_mq` | `60` | |
| `LOFREQ_CALL_SIG` | `params.lofreq_sig` | `0.01` | |
| `LOFREQ_ENABLE_INDELQUAL` | `params.lofreq_enable_indelqual` | `false` | Legacy disables for STAR spliced BAMs (crash) — keep `false` default, emit a warning if user enables |
| `LOFREQ_BAQ` | `params.lofreq_enable_baq` | `false` | Same crash rationale as above |
| `THREADS` / `MAX_JOBS` | *(removed)* | — | Replaced by Nextflow `cpus`/`max_cpus`; do not port |
| phase toggles `RUN_PHASE_*` | `params.run_ivar` / `params.run_lofreq` / `params.run_annotation` / `params.run_coverage` / `params.run_snpgenie` / `params.run_haplotype` | `true,true,true,true,false,false` | Mirrors legacy defaults exactly |
| `DATASET` | `params.dataset` | `"mouse_veev"` | Used only for `publishDir` layout; supported: `mouse_veev`, `mouse_eeev`, `rat_veev` |
| — | `params.input` | `null` (required) | Samplesheet CSV (see §4.2) |
| — | `params.outdir` | `"./results"` | |
| — | `params.publish_dir_mode` | `"copy"` | |

---

## 1. System Architecture & Directory Map

### 1.1 Design overview

The pipeline processes **per-sample RNA-seq BAMs** (from nf-core/rnaseq STAR output) against a small (~11.5 kb) alphavirus reference:

```
samplesheet.csv ──► INPUT_CHECK ──► EXTRACT_VIRAL_BAM ──┬─► IVAR_VARIANTS ─┐
                    (validate)      (samtools contig)    ├─► IVAR_CONSENSUS │
                                                         ├─► LOFREQ_CALL ──► LOFREQ_FILTER ─┐
                                                         └─► COVERAGE_DEPTH ─► COVERAGE_SUMMARIZE
                                                                                            │
              ivar_variants_to_vcf.py ◄── IVAR_VARIANTS ──► (TSV)                           │
                                                                                            ▼
                                              BCFTOOLS_CSQ (annotation, both callers) ◄── ANNOTATE
                                                                                            │
              LOFREQ_FILTER VCFs ──► SNPGenie ──► summarize/analyze/tables/limma   [optional]
              viral BAMs ──► CLIQUESNV │ VILOCA                                  [optional]
              everything ──► REPORTING (run summary, variant/coverage/haplotype plots, Excel export)
```

Design rules: one tool invocation per process; per-sample parallelism comes free from channels; optional branches (`snpgenie`, `haplotype`) are gated by params **inside the subworkflow wiring**, never inside a process `script:` block.

### 1.2 Mandatory directory map

Create **exactly** this tree. Files marked `(stub)` are created in Phase 1 as compilable skeletons and filled in later phases.

```
alphavirus-variant-analysis-workflow/
├── AGENT.md                          # this specification (do not delete)
├── README.md                         # Phase 5
├── CHANGELOG.md                      # Phase 1 (init), Phase 5 (finalize)
├── main.nf                           # entry workflow (stub → Phase 4)
├── nextflow.config                   # params, manifest, profile includes (Phase 1)
├── nextflow_schema.json              # nf-core-style parameter schema (Phase 1, minimal; Phase 5, complete)
├── nf-test.config                    # nf-test harness config (Phase 1)
├── .gitignore                        # add: work/, .nextflow*, results*/, .nf-test/
├── .editorconfig
│
├── assets/
│   ├── samplesheet.test.csv          # 2-sample test manifest
│   └── snpgenie/                     # vendored SNPGenie Perl scripts, pinned commit (Phase 3)
│       └── PINNED_COMMIT             # single line: <full commit SHA>  <tag/branch>
│
├── bin/                              # ALL custom scripts; chmod +x; shebang; strict CLI parsing
│   ├── ivar_variants_to_vcf.py       # Phase 3
│   ├── summarize_coverage.py         # Phase 3
│   ├── summarize_snpgenie.py         # Phase 3
│   ├── analyze_delta_selection.py    # Phase 3
│   ├── build_selection_tables.py     # Phase 3
│   ├── analyze_delta_limma.R         # Phase 3 (optparse)
│   ├── generate_run_summary.py       # Phase 4/5 reporting
│   ├── plot_variants.py              # sliding-window variant plots (Phase 4)
│   ├── plot_coverage.py              # Phase 4
│   ├── plot_haplotypes.py            # Phase 4
│   ├── export_excel_variants.py      # Phase 4
│   ├── parse_gb_to_gff3.py           # utility (Phase 3)
│   └── convert_gff3_to_gtf.py        # utility (Phase 3)
│
├── conf/
│   ├── base.config                   # shared process config: labels, errorStrategy, publishDir defaults
│   ├── modules.config                # withName: process-specific args + publishDir layout
│   ├── containers.config             # SINGLE SOURCE OF TRUTH for all pinned container images
│   ├── test.config                   # profile: test (Docker, tiny fixtures, resource caps)
│   ├── slurm.config                  # profile: slurm (SLURM executor + Singularity)
│   └── awsbatch.config               # profile: awsbatch (AWS Batch executor + Docker)
│
├── modules/local/
│   ├── samtools/faidx/main.nf
│   ├── extract_viral_bam/main.nf
│   ├── ivar/variants/main.nf
│   ├── ivar/consensus/main.nf
│   ├── lofreq/call/main.nf
│   ├── lofreq/filter/main.nf
│   ├── bcftools/csq/main.nf
│   ├── coverage/depth/main.nf
│   ├── coverage/summarize/main.nf
│   ├── snpgenie/run/main.nf
│   ├── cliquesnv/main.nf
│   ├── viloca/main.nf
│   ├── report/run_summary/main.nf
│   ├── report/variant_plots/main.nf
│   ├── report/coverage_plots/main.nf
│   ├── report/haplotype_plots/main.nf
│   └── report/excel_export/main.nf
│
├── subworkflows/
│   ├── input_check/main.nf
│   ├── variant_calling/main.nf       # extract → ivar + lofreq branches
│   ├── annotation/main.nf
│   ├── coverage_qc/main.nf
│   ├── selection/main.nf             # snpgenie + downstream tables
│   ├── haplotype/main.nf             # cliquesnv + viloca
│   └── reporting/main.nf
│
├── tests/
│   ├── data/
│   │   ├── README.md                 # how fixtures were generated + how to regenerate
│   │   ├── generate_fixtures.sh      # deterministic fixture generator (containerized)
│   │   ├── viral_ref.test.fasta      # small VEEV-derived reference
│   │   ├── viral_ref.test.gff3       # CDS annotation for the fixture reference
│   │   ├── sampleA.test.bam{,.bai}   # synthetic BAM with known seeded variants
│   │   └── sampleB.test.bam{,.bai}
│   ├── modules/                      # nf-test specs, mirroring modules/local tree
│   │   └── local/...  (one <name>.nf.test per module + .nf.test.snap snapshots)
│   ├── subworkflows/                 # nf-test specs per subworkflow
│   └── e2e/
│       └── main.nf.test              # full-pipeline test on fixture data
│
└── .github/workflows/
    ├── ci.yml                        # lint + container-pin checks + nf-test matrix
    └── containers.yml                # (only if custom images needed) build & push to ghcr.io
```

### 1.3 Legacy directories — final disposition

`envs/`, `Haplotype/`, `SNPGenie/`, `Variant_discovery_pipeline/`, `run_full_pipeline.sh`, `requirements.txt` are **reference-only** during migration. Do **not** delete them until Phase 5, when (after user confirmation) they are moved into `legacy/` or removed in the release commit. The new pipeline must never read from them at runtime.

---

## 2. Strict Technical Guardrails (binding for ALL phases)

### G1 — Nextflow DSL2 only
- `nextflow.enable.dsl = 2` is set (or defaulted) in `nextflow.config`; no DSL1 constructs anywhere.
- Every process lives in a file under `modules/local/.../main.nf` and declares `input`/`output`/`script` blocks. No inline processes in `main.nf` or subworkflows.
- All sample-scoped channels carry a **`meta` map** as the first element: `tuple val(meta), path(bam), path(bai)`. `meta` must contain at least `id` (sample name). Use `meta.id` for output file prefixes.
- No absolute paths, no `$HOME`, no environment-dependent paths anywhere in `modules/`, `subworkflows/`, `bin/`. References arrive only via `params` and channel inputs.

### G2 — Zero local dependencies (containers everywhere)
- **Every** process must declare a `container` directive resolving to an explicit image from `quay.io/biocontainers/<tool>:<version>--<build_hash>` with a **fully pinned tag**. The strings `latest` and untagged image names are forbidden.
- All image names are defined **once**, in `conf/containers.config` (params scope, e.g. `params.container_samtools = 'quay.io/biocontainers/samtools:1.21--h50ea8bc_0'`), and modules reference them via `${params.container_samtools}`. This gives one auditable pin file.
- **Pin resolution procedure (mandatory, Phase 2/3):** for each tool, look up valid tags at `https://quay.io/repository/biocontainers/<tool>?tab=tags` (or `skopeo list-tags docker://quay.io/biocontainers/<tool>`), pick the newest tag matching the **required version** in the table below, and record it in `conf/containers.config`. Verify each pin exists before use:
  ```bash
  docker manifest inspect "quay.io/biocontainers/<tool>:<pinned-tag>" >/dev/null && echo "PIN OK"
  ```
- **Minimum container set** (versions required by feature parity; resolve full tags per procedure above):

  | Purpose | Tool | Required version | Image |
  |---|---|---|---|
  | BAM ops, mpileup, depth, faidx | samtools | 1.21 | `quay.io/biocontainers/samtools:<pin>` |
  | Annotation (`csq`), bgzip/tabix (via htslib in bcftools img) | bcftools | 1.21 | `quay.io/biocontainers/bcftools:<pin>` |
  | LoFreq calling | lofreq | 2.1.5 | `quay.io/biocontainers/lofreq:<pin>` |
  | iVar variants/consensus | ivar | 1.4.4 | `quay.io/biocontainers/ivar:<pin>` |
  | CliqueSNV haplotypes | cliquesnv | 2.0.3 | `quay.io/biocontainers/cliquesnv:<pin>` |
  | VILOCA haplotypes | shorah (provides VILOCA) | 1.99.2 | `quay.io/biocontainers/shorah:<pin>` |
  | Python reporting (pandas, matplotlib, seaborn, openpyxl) | mulled | n/a | `quay.io/biocontainers/mulled-v2-<hash>:<pin>` (resolve via bioconda mulled index) |
  | R analysis (limma + base R) | bioconductor-limma | ≥ 3.58 | `quay.io/biocontainers/bioconductor-limma:<pin>` |
  | SNPGenie (Perl; **not** on bioconda) | perl-bioperl + vendored scripts | 1.7.8 (perl-bioperl) | `quay.io/biocontainers/perl-bioperl:<pin>`; scripts vendored under `assets/snpgenie/` at pinned commit, mounted into the process |

- If (and only if) no biocontainer/mulled image satisfies a dependency set, the fallback is a custom image: a `containers/<name>/Dockerfile` with a **pinned base image and pinned package versions**, built and pushed to `ghcr.io/<org>/<repo>/<name>:<semver>` by `.github/workflows/containers.yml`, then pinned in `conf/containers.config`. Record any such fallback in `DEVIATIONS.md`.
- **CI enforcement (Phase 1):** `tests/lint_containers.sh` must fail the build if any `container` directive outside `conf/containers.config` exists, if any image is not `quay.io/biocontainers/*` or the approved `ghcr.io` fallback, or if any tag is `latest`/missing.

### G3 — `bin/` executables
- All custom Python/R helpers live in `bin/`, are executable (`chmod +x`, shebang `#!/usr/bin/env python3` / `#!/usr/bin/env Rscript`), and are on the process `PATH` automatically (Nextflow stages `bin/`).
- **Strict CLI parsing:** Python → `argparse` (stdlib) or `click`; R → `optparse`. Every script supports `--help` and fails with exit code `2` on missing/invalid args. No positional-argument guessing; no reading of environment variables for logic.
- Every script must have a `--version` flag printing a semantic version, and must run on the pinned containers from G2 with no other dependencies.
- Windows dev note: after creating each script, run `git update-index --chmod=+x bin/<script>` (no commit) so the executable bit survives on a win32 checkout.

### G4 — Profiles (three, mandatory)

| Profile | File | Executor | Container runtime | Purpose |
|---|---|---|---|---|
| `test` | `conf/test.config` | `local` | **Docker** (`docker.enabled = true`) | CI + laptop runs on `tests/data` fixtures; `process.resourceLimits = [cpus: 2, memory: 4.GB, time: 1.h]`; `params.input = "$projectDir/assets/samplesheet.test.csv"`; fixture `params.fasta`/`params.gff`; all optional branches ON |
| `slurm` | `conf/slurm.config` | `slurm` | **Singularity** (`singularity.enabled = true`, `singularity.autoMounts = true`, `singularity.cacheDir` param) | HPC production; `process.queue`, `executor.queueSize = 50`; params for real reference paths |
| `awsbatch` | `conf/awsbatch.config` | `awsbatch` | **Docker** | Cloud production; `process.queue = '<batch-queue>'` (param `params.aws_queue`), `aws.region` (param `params.aws_region`), `aws.batch.cliPath = '/home/ec2-user/miniconda/bin/aws'` overridable, work dir on `s3://` via `-w` CLI |

Profile composition: `nextflow.config` includes `conf/base.config` and `conf/modules.config` unconditionally, and each profile file adds its layer. `-profile test,slurm` style combination is not required; each profile is self-sufficient.

### G5 — Process configuration standards (`conf/base.config` + `conf/modules.config`)
- Resource labels only: `label 'process_low'` (1 cpu / 2.GB), `'process_medium'` (4 cpu / 8.GB), `'process_high'` (8 cpu / 16.GB). No per-process hardcoded `cpus`/`memory` literals.
- `process.errorStrategy = { task.exitStatus in ((130..145) + 104) ? 'retry' : 'finish' }`, `process.maxRetries = 2`, `process.maxErrors = '-1'`.
- `publishDir` rules live **only** in `conf/modules.config` via `withName:`, using `params.publish_dir_mode`, and reproduce the legacy output layout:
  ```
  ${params.outdir}/${params.dataset}/LoFreq/<sample>/variants.filtered.vcf.gz{,.tbi}, qc_stats.txt
  ${params.outdir}/${params.dataset}/Ivar/<sample>/variants.tsv, consensus.fa
  ${params.outdir}/${params.dataset}/Annotated_variants/{LoFreq,Ivar}/...
  ${params.outdir}/${params.dataset}/Coverage/...
  ${params.outdir}/${params.dataset}/SNPGenie/...
  ${params.outdir}/${params.dataset}/Haplotypes/{CliqueSNV,VILOCA}/<sample>/...
  ${params.outdir}/${params.dataset}/Reports/{tables,Plots}/...
  ```
- Every process emits a `versions.yml` topic channel output (`topic: versions`) capturing tool name + version for the run summary.

### G6 — Testing
- All tests use **nf-test ≥ 0.9** with `nf-test.config` at repo root and `tests/nextflow.config` enabling Docker inside the harness.
- Every module gets: (a) a `-stub` test, and (b) a real-data test on fixtures asserting success + snapshot (`assert snapshot(process.out).match()`). Binary outputs (BAM) are asserted by existence + non-zero size, never by md5; text outputs (VCF/TSV) by snapshot of content with stochastic lines (headers containing paths/dates) filtered out.
- Seed test fixtures deterministically (see Phase 1 task) so snapshots are stable.

### G7 — Naming & style
- Processes: `SCREAMING_SNAKE_CASE`. Modules files/dirs: lowercase. Params/vars: `snake_case`. Channels: `ch_<noun>`.
- One logical step per process; composition only in subworkflows/`main.nf`.
- `nextflow.config` `manifest` block: `name = 'alphavirus-variant-analysis'`, `author`, `description`, `mainScript = 'main.nf'`, `nextflowVersion = '!>=24.04.0'`, `version = '1.0.0dev'` (bumped in Phase 5).

---

## 3. Phase-by-Phase Execution Roadmap

---

### PHASE 1 — Repo Scaffolding, `nextflow.config` profiles, GitHub Actions CI

**Goal:** A compilable, CI-guarded skeleton. Nothing scientific runs yet, but everything parses and lints green.

**Tasks:**
1. Create the full directory tree from §1.2, including empty-but-valid stubs:
   - `main.nf`: DSL2 header, `include` lines commented as TODO, `workflow { }` that only logs params.
   - Every module/subworkflow file: valid DSL2 skeleton (process/workflow with TODO script body that just runs `echo` and `touch`es declared outputs so `-stub-run` semantics work later).
   - All `bin/` scripts: argparse/optparse skeletons with `--help` and `--version`, no logic yet.
2. Write `nextflow.config`: `params` block with **all** defaults from §0.2; `manifest` block (G7); `includeConfig 'conf/base.config'`, `includeConfig 'conf/modules.config'`, `includeConfig 'conf/containers.config'`; `profiles { test { includeConfig 'conf/test.config' } slurm { ... } awsbatch { ... } }`.
3. Write `conf/base.config` (labels, errorStrategy, retries, `report`/`timeline`/`trace`/`dag` enabled with unique-name overwrite), `conf/modules.config` (publishDir skeleton), `conf/containers.config` (placeholder pins to be finalized in Phase 2 — but syntactically final), and the three profile files per G4.
4. `nf-test.config` + `tests/nextflow.config` (Docker enabled, work dir under `.nf-test/`).
5. **Test fixtures (deterministic):** write `tests/data/generate_fixtures.sh` that, inside pinned containers, (a) writes the small VEEV-derived `viral_ref.test.fasta` (~11.5 kb) + matching `viral_ref.test.gff3` with CDS features, and (b) synthesizes two BAMs with known seeded SNVs at known allele fractions (e.g., `dwgsim`/`wgsim` pinned container) + indexes them. Run it once, commit the outputs (`tests/data/*.bam{,.bai}`, fasta, gff3) plus `assets/samplesheet.test.csv`. Document provenance in `tests/data/README.md`. Fixture size budget: < 5 MB total.
6. Write `tests/lint_containers.sh` implementing the G2 CI enforcement (grep-based; exits non-zero on violation), and `tests/lint_structure.sh` asserting the mandatory tree exists.
7. Write `.github/workflows/ci.yml` (runs on push + PR to `main`):
   - Job `lint`: `actions/checkout@v4` → `nf-core/setup-nextflow@v2` (Nextflow `24.04.4`, Java 17) → run `tests/lint_structure.sh`, `tests/lint_containers.sh`, `nextflow config -profile test > /dev/null` (parse check), and `nextflow run . -profile test -stub-run` (must succeed on stubs).
   - Job `nftest` (Phase 1: allowed-empty matrix `["scaffold"]`): install `nf-test` (`wget .../nf-test releases/0.9.2`, chmod, add to PATH), run `nf-test test --profile docker` — passes trivially until Phase 2 adds specs; wire the matrix now (`["modules", "subworkflows", "e2e"]`) with `if: hashFiles(format('tests/{0}/**', matrix.suite)) != ''` guards.
   - Job `profile-parse`: `nextflow config -profile slurm >/dev/null` and `nextflow config -profile awsbatch >/dev/null` (validation only, no execution).
8. Write `.gitignore`, `.editorconfig`, minimal `CHANGELOG.md` (`## [Unreleased]`).

**Phase 1 DoD — run every command, all must exit 0:**

| # | Checkpoint command | Pass criterion |
|---|---|---|
| 1.1 | `bash tests/lint_structure.sh` | `STRUCTURE OK` |
| 1.2 | `bash tests/lint_containers.sh` | `CONTAINERS OK` (no `latest`, no foreign registries, no scattered `container` directives) |
| 1.3 | `nextflow config -profile test >/dev/null && nextflow config -profile slurm >/dev/null && nextflow config -profile awsbatch >/dev/null` | All three profiles parse |
| 1.4 | `nextflow run . -profile test -stub-run` | `Pipeline completed` with 0 errors on stubs |
| 1.5 | `for f in bin/*; do "$f" --help >/dev/null || exit 1; done; echo CLI-OK` | Every script has a working CLI |
| 1.6 | `nf-test --version` | `0.9.x` available |
| 1.7 | Push branch; GitHub Actions `ci.yml` | All jobs green |

---

### PHASE 2 — Core Variant Calling Modules + nf-test snapshots

**Goal:** Real, tested per-sample modules: reference indexing, viral-read extraction, iVar, LoFreq. Freeze all container pins.

**Tasks:**
1. **Finalize container pins** in `conf/containers.config` for samtools, bcftools, lofreq, ivar using the G2 resolution procedure. Run the `docker manifest inspect` verification for each and paste the resolved tags (with digests, e.g. `@sha256:...` recorded in a comment) into the config.
2. `modules/local/samtools/faidx/main.nf` — `samtools faidx`; emits `tuple val(meta_ref), path(fasta), path(fai)`.
3. `modules/local/extract_viral_bam/main.nf` — faithful port of legacy step 1: `samtools view -b <bam> <params.viral_contig> | samtools sort` → `viral_only.bam` + `samtools index`; emits `tuple val(meta), path("*.viral_only.bam"), path("*.viral_only.bam.bai")`. Handle the zero-viral-read edge: if extracted BAM has 0 reads (`samtools view -c`), still emit the (empty) BAM but set `meta.viral_reads = 0`; downstream callers must not crash on it (guard with `task.ext` or channel `filter`).
4. `modules/local/ivar/variants/main.nf` — `samtools mpileup -aa -A -d 0 -Q 0 -r <contig> --reference <fasta> <bam> | ivar variants -p <prefix> -q ${params.ivar_min_bq} -t ${params.ivar_min_freq} -m ${params.ivar_min_depth} -r <fasta> -g <gff>`; emits `variants.tsv` + `versions.yml`. `label 'process_medium'`.
5. `modules/local/ivar/consensus/main.nf` — `samtools mpileup -aa -A -d 0 -Q 0 --reference <fasta> <bam> | ivar consensus -p <prefix> -m ${params.ivar_consensus_min_cov} -t ${params.ivar_consensus_threshold} -q ${params.ivar_min_bq}`; emits `consensus.fa`.
6. `modules/local/lofreq/call/main.nf` — `lofreq call-parallel --pp-threads ${task.cpus} -f <fasta> --min-cov ${params.lofreq_min_depth} --min-bq/--min-alt-bq ${params.lofreq_min_bq} --min-mq ${params.lofreq_min_mq} --sig ${params.lofreq_sig} -o variants.vcf <bam>`. Wrap `lofreq viterbi` / `lofreq indelqual --dindel` as **optional pre-steps controlled by `task.ext.when`-style config from `params.lofreq_enable_baq` / `params.lofreq_enable_indelqual`** (implement as separate `task.ext.args`-driven command blocks inside the same module; default OFF, log a warning when ON because legacy documents STAR-BAM crashes).
7. `modules/local/lofreq/filter/main.nf` — `lofreq filter -i variants.vcf -o variants.filtered.vcf --snvqual-thresh 20 --indelqual-thresh 20` → `bgzip -f` + `tabix -f -p vcf` → emits `tuple val(meta), path("*.variants.filtered.vcf.gz"), path("*.variants.filtered.vcf.gz.tbi")`. Also compute legacy QC: viral read count + filtered variant count → `qc_stats.txt` (same 5-line format as legacy).
8. **nf-test specs** under `tests/modules/local/...` for every module above:
   - `stub` block test for each.
   - Real test using `tests/data/sampleA.test.bam`: for `lofreq/call`+`filter` assert `process.success` and snapshot the filtered VCF body (strip `##fileDate`/path-bearing header lines before `md5` or use `snapshot(...).match()` on curated lines); assert at least one seeded fixture variant is recovered (`assert vcf.contains('KP282671.1')` style content assertion via `path(...).linesGzip`).
   - iVar variants: assert TSV header + ≥1 data line; consensus: assert FASTA exists, length ≈ reference length.
9. Update `ci.yml` `nftest` matrix guards so `tests/modules/**` now runs.

**Phase 2 DoD:**

| # | Checkpoint command | Pass criterion |
|---|---|---|
| 2.1 | `grep -E "quay.io/biocontainers/(samtools|bcftools|lofreq|ivar):" conf/containers.config \| wc -l` | `≥ 4` distinct pinned images, zero `latest` |
| 2.2 | `nf-test test tests/modules/local/extract_viral_bam/main.nf.test --profile docker` | PASS |
| 2.3 | `nf-test test tests/modules/local/ivar --profile docker` | All iVar specs PASS |
| 2.4 | `nf-test test tests/modules/local/lofreq --profile docker` | All LoFreq specs PASS |
| 2.5 | `nf-test test tests/modules --profile docker` | Whole module suite green, snapshots committed (`*.snap` files exist) |
| 2.6 | `bash tests/lint_containers.sh` | Still `CONTAINERS OK` |
| 2.7 | CI run on pushed branch | `lint` + `nftest[modules]` green |

---

### PHASE 3 — Downstream Analytics Modules (annotation, coverage, SNPGenie, haplotypes, `bin/` scripts)

**Goal:** Port every remaining analytical step; all `bin/` scripts become strict-CLI, container-compatible tools.

**Tasks:**
1. **Port `bin/` scripts with strict CLIs** (translate legacy logic 1:1; improve only argument handling and path-independence):
   - `bin/ivar_variants_to_vcf.py` — `--input-tsv --output-vcf --reference-fasta` (argparse).
   - `bin/summarize_coverage.py` — consumes `samtools depth` output; `--depth-file --sample-id --output-tsv [--windowsize 500]`.
   - `bin/summarize_snpgenie.py`, `bin/analyze_delta_selection.py`, `bin/build_selection_tables.py` — ports of the three legacy Python helpers; each takes explicit `--input-dir/--output-dir` (+ `--manifest` where legacy used one).
   - `bin/analyze_delta_limma.R` — port of legacy R script using `optparse`; pin behavior to the `bioconductor-limma` container.
   - `bin/parse_gb_to_gff3.py`, `bin/convert_gff3_to_gtf.py` — utility ports (not wired into the DAG; documented in README).
2. `modules/local/bcftools/csq/main.nf` — `bcftools csq -f <fasta> -g <gff> --local-csq <vcf> -o <prefix>.csq.vcf`; used for LoFreq VCFs directly and for iVar via a preceding `IVAR_TSV_TO_VCF` module that calls `ivar_variants_to_vcf.py` (same G2 Python container as reporting).
3. `modules/local/coverage/depth/main.nf` (`samtools depth -a`) + `modules/local/coverage/summarize/main.nf` (`summarize_coverage.py`).
4. **SNPGenie (vendored, pinned):** SNPGenie is not on bioconda. Vendor `snpgenie.pl` and its `bin/` Perl dependencies from `github.com/chasewnelson/SNPGenie` into `assets/snpgenie/`, record the exact commit SHA in `assets/snpgenie/PINNED_COMMIT`, and add a SHA-256 checksum file (`checksums.sha256`) verified by `tests/lint_structure.sh`. Module `modules/local/snpgenie/run/main.nf` runs it inside the pinned `perl-bioperl` container: `perl ${projectDir}/assets/snpgenie/snpgenie.pl --vcfformat=2 --snpreport=<in.vcf> --fastafile=<fasta> --gtffile=<gtf>` (port legacy flags exactly; GTF produced by `bin/convert_gff3_to_gtf.py` in a tiny prep module). Emits per-sample `*_results.tsv`/`product_results` files mirroring legacy outputs.
5. `modules/local/cliquesnv/main.nf` — port of `run_cliquesnv.sh` invocation (Java jar call exactly as legacy, using `task.cpus` for threads and `task.memory` for `-Xmx`); emits haplotype FASTA/TSV outputs.
6. `modules/local/viloca/main.nf` — port of `run_viloca.sh` (`shorah shotgun` VILOCA mode flags exactly as legacy); emits `cooccurring_mutations.csv`, `coverage.txt`, haplotype FASTAs.
7. **nf-test specs** for every Phase 3 module + a CLI smoke test for every `bin/` script (`--help` exit 0, `--version` exit 0, invalid args exit 2). SNPGenie/CliqueSNV/VILOCA tests run on the tiny fixture inputs and assert output existence + key content lines; mark them `tag 'slow'` if runtime > 3 min so CI can select with `--tag`/exclude as needed — but they must pass in CI on the fixtures.

**Phase 3 DoD:**

| # | Checkpoint command | Pass criterion |
|---|---|---|
| 3.1 | `for f in bin/*; do "$f" --help >/dev/null 2>&1 \|\| { echo "FAIL $f"; exit 1; }; "$f" --badflag >/dev/null 2>&1 && { echo "NO-STRICT $f"; exit 1; }; done; echo BIN-OK` | All CLIs strict |
| 3.2 | `cat assets/snpgenie/PINNED_COMMIT && (cd assets/snpgenie && sha256sum -c checksums.sha256)` | Vendored SNPGenie integrity verified |
| 3.3 | `nf-test test tests/modules/local/bcftools tests/modules/local/coverage --profile docker` | PASS |
| 3.4 | `nf-test test tests/modules/local/snpgenie tests/modules/local/cliquesnv tests/modules/local/viloca --profile docker` | PASS |
| 3.5 | `nf-test test tests/modules --profile docker` | Entire module suite green |
| 3.6 | `bash tests/lint_containers.sh && bash tests/lint_structure.sh` | Both OK |
| 3.7 | CI run | `lint` + `nftest[modules]` green |

---

### PHASE 4 — Subworkflow Assembly & Channel Plumbing in `main.nf`

**Goal:** Wire everything into a coherent, resumable DAG driven by a samplesheet.

**Tasks:**
1. `subworkflows/input_check/main.nf` — parse `params.input` CSV via `Channel.fromSamplesheet` (nf-core style, header `sample,bam,bai`); validate: files exist, BAM index present, `params.fasta`/`params.gff` provided, GFF3 contains `CDS` (fail fast with `error()` and a remediation message — port of legacy `annotate_all.sh` guard); emit `ch_samples = [meta, bam, bai]`, `ch_fasta`, `ch_gff`.
2. `subworkflows/variant_calling/main.nf` — `EXTRACT_VIRAL_BAM` → fan-out to `IVAR_VARIANTS`, `IVAR_CONSENSUS` (if `params.run_ivar`), `LOFREQ_CALL → LOFREQ_FILTER` (if `params.run_lofreq`). Zero-read samples bypass callers via channel `filter { meta, bam, bai -> meta.viral_reads > 0 }` with a logged warning. Emit per-caller output channels + `versions`.
3. `subworkflows/annotation/main.nf` — iVar TSV → VCF conversion, then `BCFTOOLS_CSQ` on both callers' VCFs (if `params.run_annotation`); `mix` the two annotated streams with `meta.caller = 'ivar'|'lofreq'`.
4. `subworkflows/coverage_qc/main.nf` — depth + summarize on viral BAMs (if `params.run_coverage`).
5. `subworkflows/selection/main.nf` — (if `params.run_snpgenie`) GTF prep → `SNPGENIE_RUN` per sample → `summarize_snpgenie.py` (collect) → `analyze_delta_selection.py` → `analyze_delta_limma.R` → `build_selection_tables.py`; manifest channel carries `meta.dpi`/`meta.condition` columns from the samplesheet (extend samplesheet schema: `sample,bam,bai,condition,dpi`).
6. `subworkflows/haplotype/main.nf` — (if `params.run_haplotype`) `CLIQUESNV` and `VILOCA` in parallel branches on viral BAMs.
7. `subworkflows/reporting/main.nf` — run summary (`generate_run_summary.py`), sliding-window variant plots (`plot_variants.py`), coverage plots (`plot_coverage.py`), haplotype plots (`plot_haplotypes.py` — input channel may be empty when haplotypes off; script must exit 0 with "no data" notice), consolidated Excel export (`export_excel_variants.py`). Use `collect()`/`groupTuple()` so reporting runs once per dataset, not per sample.
8. `main.nf` — the only wiring file:
   ```groovy
   workflow {
       INPUT_CHECK() → VARIANT_CALLING → ANNOTATION → COVERAGE_QC
                    ↘ (optional) SELECTION   ↘ (optional) HAPLOTYPE
                    ↘ REPORTING (gathers everything above)
       // software_versions.yml aggregation via Channel.topic('versions')
   }
   ```
   Print the legacy-style startup banner (dataset, reference, contig, enabled phases) via `log.info`.
9. `conf/modules.config` — complete all `withName:` `publishDir` rules to reproduce the §G5 layout; add `ext.args`/`ext.prefix` where modules consume them.
10. **nf-test**: one spec per subworkflow under `tests/subworkflows/`, plus `tests/e2e/main.nf.test` running the full pipeline on fixtures with **all branches on** (`run_snpgenie=true, run_haplotype=true`) asserting: pipeline success, expected published files exist (spot-check ~10 paths across the layout), filtered LoFreq VCF for `sampleA` contains a seeded fixture variant, `software_versions.yml` lists every tool.

**Phase 4 DoD:**

| # | Checkpoint command | Pass criterion |
|---|---|---|
| 4.1 | `nextflow run . -profile test -stub-run` | Full DAG stub-runs clean |
| 4.2 | `nextflow run . -profile test --outdir results_test` | Real end-to-end run on fixtures: exit 0 |
| 4.3 | `test -f results_test/mouse_veev/LoFreq/sampleA/variants.filtered.vcf.gz && test -f results_test/mouse_veev/Annotated_variants/LoFreq/sampleA_filtered.vcf && test -f results_test/mouse_veev/Reports/tables/consolidated_sample_summary.tsv && echo LAYOUT-OK` | Published layout matches spec |
| 4.4 | `nextflow run . -profile test --outdir results_test2 -resume` | Resume run completes with ≥1 `cached` process and 0 failures |
| 4.5 | `nf-test test tests/subworkflows tests/e2e --profile docker` | PASS |
| 4.6 | `nf-test test --profile docker` | **Entire** test suite green |
| 4.7 | `nextflow run . -profile test --run_snpgenie false --run_haplotype false -stub-run` | Optional branches correctly excluded |
| 4.8 | CI run | All jobs green including `nftest[subworkflows, e2e]` |

---

### PHASE 5 — Documentation, Dry-Runs, and v1.0.0 Release Preparation

**Goal:** Hand-off quality. A new user can run the pipeline from README alone; release artifacts prepared.

**Tasks:**
1. **README.md** (repo root) with: badges (CI status), one-paragraph summary, DAG diagram (embed the Phase 4 `dag` PNG or a Mermaid diagram), Quick start (`nextflow run . -profile test`), full usage section (samplesheet format table incl. `condition`/`dpi`, all `params` from §0.2 with defaults in a table), the three profiles with copy-paste commands (`-profile test`, `-profile slurm`, `-profile awsbatch -w s3://bucket/work`), output directory layout, container policy note, how to regenerate test fixtures, migration note mapping legacy scripts → modules (from §0.1), and citation/license stubs.
2. Complete `nextflow_schema.json` for all params (helps `--help` and tower/launchpad validation).
3. Finalize `CHANGELOG.md`: move `Unreleased` → `## [1.0.0] - <date>` with a full feature list.
4. Bump `manifest.version` in `nextflow.config`: `1.0.0dev` → `1.0.0`.
5. **Dry-run matrix (no real cluster needed):** validate the production profiles parse and stub-resolve:
   - `nextflow run . -profile slurm -stub-run` with `docker` disabled and Singularity absent → must at minimum pass `nextflow config` and `-stub-run` on a machine with Singularity, or document the limitation; run `nextflow config -profile slurm` and confirm `process.executor = 'slurm'` and `singularity.enabled = true` appear in resolved output.
   - Same for `awsbatch`: confirm `executor = 'awsbatch'`, `docker.enabled`, `queue`, `region` resolve.
6. Fresh-clone reproducibility check: `git clone` (or `git archive`) the repo into a clean temp dir and run the full test profile there — catches accidental dependence on untracked local files.
7. **Release tag preparation (requires explicit user confirmation before any git mutation):** prepare, but do not execute without confirmation, the exact commands:
   ```bash
   git add -A && git commit -m "chore(release): 1.0.0"
   git tag -a v1.0.0 -m "v1.0.0: Nextflow DSL2 port of alphavirus variant analysis workflow"
   git push origin main --tags
   ```
   Present them to the user for approval. After approval and execution, verify `git describe --tags` prints `v1.0.0`.
8. Update `AGENT.md` §1.3 once legacy directories are archived per user decision.

**Phase 5 DoD:**

| # | Checkpoint command | Pass criterion |
|---|---|---|
| 5.1 | `grep -q "version = '1.0.0'" nextflow.config && echo VERSION-OK` | Release version set |
| 5.2 | `nextflow config -profile slurm \| grep -E "executor = 'slurm'" && nextflow config -profile awsbatch \| grep -E "executor = 'awsbatch'"` | Both production profiles resolve correctly |
| 5.3 | `tmp=$(mktemp -d); git archive HEAD \| tar -x -C "$tmp"; (cd "$tmp" && nextflow run . -profile test --outdir results_clean) && echo CLEAN-CLONE-OK` | Pipeline runs from pristine export |
| 5.4 | `nf-test test --profile docker` | Final full suite green |
| 5.5 | `bash tests/lint_containers.sh && bash tests/lint_structure.sh` | Final lints green |
| 5.6 | README renders: all Quick-start commands copy-paste runnable (spot-check `-profile test` from README verbatim) | OK |
| 5.7 | `git describe --tags` (after user-approved tagging) | `v1.0.0` |
| 5.8 | CI on `main` post-merge | All jobs green |

---

## 4. Cross-Cutting Reference

### 4.1 Global Definition of Done (project level)

- [ ] All 5 phase DoD tables passed, in order.
- [ ] `nf-test test --profile docker` green from a clean clone.
- [ ] Zero occurrences of `latest`, zero non-biocontainer/approved-ghcr images (`tests/lint_containers.sh` enforces in CI).
- [ ] Zero absolute paths / conda references in `main.nf`, `modules/`, `subworkflows/`, `conf/`, `bin/`.
- [ ] All three profiles (`test`, `slurm`, `awsbatch`) parse and are documented in README.
- [ ] Legacy parameter defaults from §0.2 preserved exactly.
- [ ] Output layout matches §G5 for the `mouse_veev` test dataset.
- [ ] `v1.0.0` tag created **only after explicit user confirmation**.

### 4.2 Samplesheet schema (`params.input`)

```csv
sample,bam,bai,condition,dpi
sampleA,tests/data/sampleA.test.bam,tests/data/sampleA.test.bam.bai,infected,3
sampleB,tests/data/sampleB.test.bam,tests/data/sampleB.test.bam.bai,infected,7
```
`condition`/`dpi` are required only when `params.run_snpgenie = true`; `INPUT_CHECK` must validate this conditionally.

### 4.3 Commands the Executor runs most often

```bash
nextflow run . -profile test -stub-run                      # fast DAG sanity
nextflow run . -profile test --outdir results_test          # local e2e on fixtures
nf-test test --profile docker                               # full test suite
nf-test test tests/modules/local/lofreq --profile docker    # focused suite
bash tests/lint_containers.sh && bash tests/lint_structure.sh
nextflow config -profile slurm | less                       # inspect resolved profile
```

### 4.4 Standing prohibitions (repeat of hard rules)

- No DSL1 syntax. No `channel.from` deprecated forms.
- No `conda`, no local tool assumptions, no `latest` tags, no unpinned images.
- No per-process resource literals (labels only). No absolute paths.
- No git mutations without explicit, at-the-moment user confirmation.
- No skipping a DoD checkpoint to "save time". If blocked, stop and report.

---

*End of specification. Executor: begin at Phase 1, Task 1.*
