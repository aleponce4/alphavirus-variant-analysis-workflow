# AGENT_UPDATE.md — Executor Specification: Phases 4–7 (Production Hardening → Subworkflow Toggles → MultiQC Release → Generalization & UX)

> **Document type:** Binding implementation specification.
> **Audience:** Executor Agent (coding model).
> **Relationship to `AGENT.md`:** This document **extends and partially supersedes** `AGENT.md`. The original Phases 1–3 (scaffolding, core callers, downstream analytics) and the first-pass subworkflow assembly are **complete and released as `v1.0.0`** (tag exists). This document defines the **next four phases (4–7)** on top of that codebase. Where this document conflicts with `AGENT.md`, **this document wins**. Where it is silent, `AGENT.md` guardrails (G1–G7) still apply in full.
> **Hard rule carried over:** **Never** run `git commit`, `git tag`, `git push`, or any other git mutation without explicit user confirmation at the moment it is needed.

---

## 0. Standing Rules (carried over from AGENT.md §0, still binding)

1. Execute Phases 4 → 5 → 6 → 7 **strictly sequentially**. Do not start Phase N+1 until every checkpoint in Phase N passes with exit code `0`.
2. Each phase ends with a **Definition of Done (DoD)** table of shell commands. Run them verbatim, capture output, proceed only when all pass.
3. If a checkpoint fails twice after attempted fixes, **stop and report** the failing command, output, and diagnosis. Do not weaken or skip checks.
4. **Minimal diffs.** Retrofit and extend — never rewrite working modules/subworkflows wholesale. Existing committed nf-test snapshots (`*.snap`) must keep passing unless a checkpoint in this document explicitly authorizes regenerating them.
5. Record every forced deviation from this spec in `DEVIATIONS.md` at repo root.
6. Windows dev note: after creating or moving any `bin/` script, run `git update-index --chmod=+x bin/<script>` (no commit) so the executable bit survives on win32 checkouts.

### 0.1 Current-State Assessment (verified against repo @ `v1.0.0`)

The Executor must treat this table as ground truth. Do not re-implement anything marked ✅; extend it.

| Requirement (this spec) | Current repo state | Action |
|---|---|---|
| `check_max()` in `nextflow.config` | ⚠️ Defined in `conf/base.config` instead; works, but wrong home (nf-core canonical = `nextflow.config`) | **Move** function to `nextflow.config`; keep call sites in `conf/base.config` |
| Retry on exit 130–145 + 104, dynamic `task.attempt` scaling | ✅ Already in `conf/base.config` (labels `process_low/medium/high`) | **Verify only**, do not touch logic |
| `nf-schema@2.1.1` plugin validation | ❌ No `plugins` block; legacy draft-07 `nextflow_schema.json`; hand-rolled `splitCsv` in `INPUT_CHECK` | **New** (Phase 4) |
| Samplesheet `sample,fastq_1,fastq_2,treatment` | ❌ Current sheet is `sample,bam,bai,condition,dpi` (BAM ingress) | **Migrate to FASTQ ingress** (Phase 4): new `FASTQC`/`BWA_INDEX`/`BWA_MEM`/`SAMTOOLS_STATS`/`IVAR_TRIM` modules + `read_preprocessing` subworkflow; downstream `VARIANT_CALLING` interface **unchanged** (`[meta, bam, bai]`) |
| Conditional `--protocol` / `--primer_bed` | ❌ Not present | **New** params + schema conditional + `IVAR_TRIM` branch gated by `params.protocol == 'amplicon'` |
| Modules emit `versions.yml` | ✅ All 21 module files emit `path "versions.yml", emit: versions` (named emit, `END_VERSIONS` heredoc) | **Audit + standardize** (Phase 4): every process must emit it in `script:` **and** `stub:`; every subworkflow must re-emit a mixed `versions` channel (`INPUT_CHECK` currently does **not** — fix) |
| `DUMP_SOFTWARE_VERSIONS` aggregation | ❌ `main.nf` currently **ignores all `versions` outputs** | **New module + plumbing** (Phase 6) |
| `MULTIQC` terminal step | ❌ Not present | **New** (Phase 6), fed by FastQC/BWA/Samtools/caller logs created in Phase 4 |
| `subworkflows/local/selection.nf`, `subworkflows/local/haplotype.nf` | ⚠️ Logic exists at `subworkflows/selection/main.nf` and `subworkflows/haplotype/main.nf`, flag-gated internally by `params.run_snpgenie` / `params.run_haplotype` | **Relocate** (Phase 5) — move files, do not rewrite logic |
| `.collect()` / `.groupTuple()` wiring in `main.nf` | ⚠️ `.collect()` exists *inside* `SELECTION`; `main.nf` has no grouping logic | **Add** treatment `groupTuple()` + collect wiring (Phase 5) |
| Executive HTML report (Quarto) | ❌ `REPORTING` emits raw PNG/TSV/Excel artifacts only; no consolidated human-readable summary | **New** (Phase 7): `assets/executive_report.qmd` + terminal `EXECUTIVE_REPORT` step wrapping existing report outputs — complements MultiQC (QC) with a findings report (biology) |
| Hosted "tiny" test dataset | ⚠️ Fixtures are committed locally under `tests/data/` (< 5 MB) — works after clone, but couples demo data to the repo | **New** (Phase 7): ~50k-read public dataset in a dedicated repo, SHA-pinned raw URLs hardcoded into `conf/test.config`; local fixtures retained for hermetic nf-test |
| Pipeline identity / scope | ⚠️ Branded `alphavirus-variant-analysis`, but the design is virus-agnostic (reference FASTA/GFF are user params; nothing alphavirus-specific in the DAG) | **Rename** to generic intrahost-virus identity + add compatibility note: validated on alphaviruses, **not tested for highly recombinant/hypervariable viruses** (Phase 7) |
| Release tag `v1.0.0` | ⚠️ **Already exists** (pre-hardening cut) | **Deviation:** cut **`v1.1.0`** in Phase 6; never move/retag `v1.0.0`. Log in `DEVIATIONS.md` |

---

## 1. Target Architecture (post-update)

```
samplesheet.csv ──► INPUT_CHECK ──► READ_PREPROCESSING ──► VARIANT_CALLING ──► ANNOTATION
 (fastq_1/fastq_2,   (nf-schema      ├ FASTQC (QC)         ├ EXTRACT_VIRAL_BAM   (bcftools csq)
  treatment)          validation)    ├ BWA_INDEX→BWA_MEM   ├ IVAR_VARIANTS/CONSENSUS
                                     ├ IVAR_TRIM*          ├ LOFREQ_CALL→FILTER ──► COVERAGE_QC
                                     └ SAMTOOLS_STATS      │        │
                                      *amplicon only       │        ├─► SELECTION   (--run_snpgenie)
                                                           │        │   SNPGenie→delta→limma→tables
                                                           │        └─► HAPLOTYPE   (--run_haplotype)
                                                           │            CliqueSNV ∥ VILOCA
                                                           ▼
                                              DUMP_SOFTWARE_VERSIONS ──► MULTIQC (terminal)
                                                           ▼
                                              EXECUTIVE_REPORT (Quarto, self-contained HTML; Phase 7)
```

New/changed files only (everything else stays as-is):

```
assets/schema_input.json                 # Phase 4 — nf-schema samplesheet schema
assets/multiqc_config.yml                # Phase 6
assets/samplesheet.test.csv              # Phase 4 — REWRITE to fastq/treatment columns
modules/local/fastqc/main.nf             # Phase 4
modules/local/bwa/index/main.nf          # Phase 4
modules/local/bwa/mem/main.nf            # Phase 4
modules/local/samtools/stats/main.nf     # Phase 4
modules/local/ivar/trim/main.nf          # Phase 4
modules/local/dumpsoftwareversions/main.nf  # Phase 6
modules/local/multiqc/main.nf            # Phase 6
subworkflows/local/read_preprocessing.nf # Phase 4
subworkflows/local/selection.nf          # Phase 5 (moved from subworkflows/selection/main.nf)
subworkflows/local/haplotype.nf          # Phase 5 (moved from subworkflows/haplotype/main.nf)
tests/data/sample{A,B}.test_{1,2}.fastq.gz  # Phase 4 fixtures
tests/data/primers.test.bed              # Phase 4 fixture (amplicon mode)
tests/modules/local/{fastqc,bwa/*,samtools/stats,ivar/trim,dumpsoftwareversions,multiqc}/main.nf.test
```

---

## 2. PHASE 4 — Production Hardening (Retrofitting)

**Goal:** Cloud-native resiliency in the configs, strict nf-schema input validation, FASTQ ingress with a read-preprocessing branch, and a fully standardized software-provenance surface — without disturbing the released variant-calling logic.

### Task 4.1 — Centralize `check_max()` in `nextflow.config`

1. **Cut** the entire `def check_max(obj, type) { ... }` function from `conf/base.config` and **paste it at the bottom of `nextflow.config`** (after the `dag { }` block), unchanged. Rationale: nf-core canonical layout; one definition, visible to all included configs at task-evaluation time.
2. Do **not** reorder `includeConfig 'conf/base.config'` — it must remain where it is. The closures referencing `check_max` in `conf/base.config` evaluate lazily at task runtime, so a later definition in the including file is correct.
3. `conf/base.config` must end with the `process { }` block and nothing else.

### Task 4.2 — Verify (do not redesign) retry & dynamic scaling in `conf/base.config`

Confirm each of the following is present and unchanged; fix only if missing:

- `process.errorStrategy = { task.exitStatus in ((130..145) + 104) ? 'retry' : 'finish' }` — retries SIGTERM/OOM-kill range (130–145) and Spot/preemption (104); `finish` otherwise.
- `process.maxRetries = 2`, `process.maxErrors = '-1'`.
- Dynamic escalation on all three labels, exactly of the form (example): `memory = { check_max( 8.GB * task.attempt, 'memory' ) }`, `cpus = { check_max( 4 * task.attempt, 'cpus' ) }`, `time = { check_max( 6.h * task.attempt, 'time' ) }`.

### Task 4.3 — Strict input validation with `nf-schema@2.1.1`

1. **`nextflow.config`** — add at top level, immediately after `nextflow.enable.dsl = 2`:
   ```groovy
   plugins {
       id 'nf-schema@2.1.1'
   }
   ```
2. **New params** in `nextflow.config` `params { }` block (with the other input/output options):
   ```groovy
   protocol      = 'metagenomic'   // 'metagenomic' | 'amplicon'
   primer_bed    = null            // required when protocol == 'amplicon'
   ```
3. **`assets/schema_input.json`** (new file, JSON Schema draft 2020-12):
   ```json
   {
       "$schema": "https://json-schema.org/draft/2020-12/schema",
       "$id": "https://raw.githubusercontent.com/aleponce4/alphavirus-variant-analysis-workflow/main/assets/schema_input.json",
       "title": "alphavirus-variant-analysis samplesheet schema",
       "description": "Schema for the file provided with --input",
       "type": "array",
       "items": {
           "type": "object",
           "properties": {
               "sample":    { "type": "string", "pattern": "^\\S+$", "meta": ["id"], "errorMessage": "Sample name must be provided and cannot contain spaces" },
               "fastq_1":   { "type": "string", "format": "file-path", "pattern": "^\\S+\\.f(ast)?q(\\.gz)?$", "errorMessage": "Reads 1 FastQ must exist and end in .fq/.fastq (optionally .gz)" },
               "fastq_2":   { "type": "string", "format": "file-path", "pattern": "^\\S+\\.f(ast)?q(\\.gz)?$", "errorMessage": "Reads 2 FastQ must exist and end in .fq/.fastq (optionally .gz)" },
               "treatment": { "type": "string", "pattern": "^\\S+$", "errorMessage": "Treatment group must be provided and cannot contain spaces" }
           },
           "required": ["sample", "fastq_1", "fastq_2", "treatment"]
       },
       "uniqueEntries": ["sample"]
   }
   ```
4. **`nextflow_schema.json`** — update in place (keep all existing definitions):
   - Bump top-level `"$schema"` to `"https://json-schema.org/draft/2020-12/schema"`.
   - On the `input` property: add `"schema": "assets/schema_input.json"` and change the description to *"Path to comma-separated samplesheet with columns `sample,fastq_1,fastq_2,treatment`."*
   - Add a new definition `protocol_options` with:
     - `protocol`: `{"type": "string", "default": "metagenomic", "enum": ["metagenomic", "amplicon"], "description": "Sequencing protocol. 'amplicon' enables iVar primer trimming."}`
     - `primer_bed`: `{"type": "string", "format": "file-path", "description": "Primer BED file. Required when --protocol amplicon."}`
   - Enforce the conditional with draft 2020-12 `if/then` at the top level of the schema:
     ```json
     "allOf": [
         {
             "if":   { "properties": { "protocol": { "const": "amplicon" } } },
             "then": { "required": ["primer_bed"] }
         }
     ]
     ```
   - Add `protocol` and `primer_bed` to the top-level `properties` map via `$ref` like the other params.
5. **`main.nf`** — activate the plugin (top of file, after `nextflow.enable.dsl`):
   ```groovy
   include { validateParameters; paramsSummaryLog; samplesheetToList } from 'plugin/nf-schema'
   ```
   and as the first statements inside `workflow { }` (before the banner):
   ```groovy
   validateParameters()
   log.info paramsSummaryLog(workflow)
   ```
   Add `Protocol` and `Primer BED` lines to the startup banner.
6. **`subworkflows/input_check/main.nf`** — replace the hand-rolled `splitCsv` block with validated parsing (keep the existing `params.fasta`/`params.gff`/CDS guards and `SAMTOOLS_FAIDX` wiring untouched):
   ```groovy
   ch_samplesheet = Channel
       .fromList(samplesheetToList(params.input, "${projectDir}/assets/schema_input.json"))
       .map { row ->
           def meta = [:]
           meta.id        = row.sample
           meta.treatment = row.treatment
           return [ meta, file(row.fastq_1, checkIfExists: true), file(row.fastq_2, checkIfExists: true) ]
       }
   ```
   Add belt-and-braces runtime guards (schema is primary; these give actionable CLI errors):
   ```groovy
   if (params.protocol == 'amplicon' && !params.primer_bed) {
       error "Parameter --primer_bed is required when --protocol amplicon."
   }
   if (params.protocol == 'metagenomic' && params.primer_bed) {
       log.warn "--primer_bed supplied but --protocol is 'metagenomic'; primer trimming is OFF."
   }
   ```
   Change the `samples` emit to `fastqs` (`[meta, fastq_1, fastq_2]`) and **add a `versions` emit** mixing `SAMTOOLS_FAIDX.out.versions` (currently missing — provenance gap).
7. **Samplesheet + scripts column migration:** the selection `bin/` scripts (`analyze_delta_selection.py`, `analyze_delta_limma.R`, `build_selection_tables.py`, `summarize_snpgenie.py`) read the manifest's `condition`/`dpi` columns. Update them minimally: prefer `treatment`, fall back to `condition` (alias); treat `dpi` as optional — when absent, skip dpi-stratified logic and emit the non-stratified outputs. Each script's `--help`/tests must still pass. **Do not** change their analytical logic.

### Task 4.4 — FASTQ ingress: read-preprocessing modules + subworkflow

All new modules follow existing conventions exactly: `tag "$meta.id"`, resource label (no literals), `container "${params.container_<tool>}"`, `when: task.ext.when == null || task.ext.when`, `def prefix = task.ext.prefix ?: "${meta.id}"`, and the `cat <<-END_VERSIONS > versions.yml` heredoc in **both** `script:` and `stub:`.

1. **Container pins** — append to `conf/containers.config`, resolving tags per the G2 procedure (`docker manifest inspect ... && echo "PIN OK"`, record digest in a comment). Suggested (verify before use): fastqc `0.12.1`, bwa `0.7.18`. Params: `container_fastqc`, `container_bwa`. (`SAMTOOLS_STATS` reuses `container_samtools`; `IVAR_TRIM` reuses `container_ivar`.)
2. `modules/local/fastqc/main.nf` — `fastqc --quiet --threads ${task.cpus} ${fastq_1} ${fastq_2}`; emits `tuple val(meta), path("*.html"), emit: html`, `tuple val(meta), path("*.zip"), emit: zip`, `versions`. Label `process_low`.
3. `modules/local/bwa/index/main.nf` — `bwa index ${fasta}`; emits `tuple val(meta_ref), path(fasta), path("bwa")` (stage index files into a `bwa/` subdir) + versions. Label `process_low`.
4. `modules/local/bwa/mem/main.nf` — `bwa mem -t ${task.cpus} ${index_dir}/${fasta} ${fastq_1} ${fastq_2} 2> ${prefix}.bwa.log | samtools sort -@ ${task.cpus} -o ${prefix}.bam -` then `samtools index ${prefix}.bam`; emits `tuple val(meta), path("*.bam"), path("*.bam.bai"), emit: bam`, `tuple val(meta), path("*.bwa.log"), emit: log`, `versions`. Label `process_medium`.
5. `modules/local/ivar/trim/main.nf` — canonical amplicon step: `ivar trim -i ${bam} -b ${primer_bed} -p ${prefix}.trimmed -q ${params.ivar_min_bq} -m 30 -e` → `samtools sort -o ${prefix}.trimmed.sorted.bam ${prefix}.trimmed.bam` → `samtools index`; emits `tuple val(meta), path("*.trimmed.sorted.bam"), path("*.trimmed.sorted.bam.bai"), emit: bam`, versions. Label `process_medium`.
6. `modules/local/samtools/stats/main.nf` — on the **final per-sample aligned BAM** (post-trim when amplicon, else BWA output): `samtools stats`, `samtools flagstat`, `samtools idxstats`; emits one channel `tuple val(meta), path("*.stats"), path("*.flagstat"), path("*.idxstats"), emit: stats` + versions. Label `process_low`.
7. **`subworkflows/local/read_preprocessing.nf`** (new; flat `subworkflows/local/` convention):
   ```groovy
   workflow READ_PREPROCESSING {
       take:
       ch_fastqs   // [ val(meta), path(fastq_1), path(fastq_2) ]
       ch_fasta    // [ val(meta_ref), path(fasta), path(fai) ]
       main:
       ch_versions = Channel.empty()
       FASTQC(ch_fastqs)
       BWA_INDEX(ch_fasta.map { meta, fasta, fai -> [ meta, fasta ] })
       BWA_MEM(ch_fastqs, BWA_INDEX.out.index)
       ch_bams = BWA_MEM.out.bam
       if (params.protocol == 'amplicon') {
           IVAR_TRIM(BWA_MEM.out.bam, file(params.primer_bed, checkIfExists: true))
           ch_bams = IVAR_TRIM.out.bam
           ch_versions = ch_versions.mix(IVAR_TRIM.out.versions)
       }
       SAMTOOLS_STATS(ch_bams)
       ch_versions = ch_versions
           .mix(FASTQC.out.versions)
           .mix(BWA_INDEX.out.versions)
           .mix(BWA_MEM.out.versions)
           .mix(SAMTOOLS_STATS.out.versions)
       emit:
       bams        = ch_bams                  // [ meta, bam, bai ]
       fastqc_zip  = FASTQC.out.zip
       bwa_log     = BWA_MEM.out.log
       stats       = SAMTOOLS_STATS.out.stats
       versions    = ch_versions
   }
   ```
   Gate `IVAR_TRIM` in the wiring (as shown), never inside the process `script:` (AGENT.md §1.1 design rule).
8. **`main.nf` rewiring** (minimal):
   ```groovy
   include { READ_PREPROCESSING } from './subworkflows/local/read_preprocessing'
   ...
   INPUT_CHECK(file(params.input, checkIfExists: true))
   READ_PREPROCESSING(INPUT_CHECK.out.fastqs, INPUT_CHECK.out.fasta)
   VARIANT_CALLING(READ_PREPROCESSING.out.bams, INPUT_CHECK.out.fasta, INPUT_CHECK.out.gff)
   ```
   Every other `main.nf` line stays as-is. `VARIANT_CALLING`, `EXTRACT_VIRAL_BAM`, and all downstream modules/tests keep their `[meta, bam, bai]` interface — **zero churn below this seam**.
9. **Fixtures (deterministic, < 5 MB total new):**
   - Extend `tests/data/generate_fixtures.sh` to additionally emit `sampleA.test_1.fastq.gz`, `sampleA.test_2.fastq.gz`, `sampleB.test_1.fastq.gz`, `sampleB.test_2.fastq.gz` **derived from the committed synthetic BAMs** (`samtools fastq -1/-2` inside the pinned samtools container) so e2e results stay consistent with the BAM-fixture module tests. Run once, commit outputs, update `tests/data/README.md`.
   - Write `tests/data/primers.test.bed`: 2–3 synthetic primer rows spanning the fixture reference's seeded-variant regions (6-column BED).
   - **Rewrite** `assets/samplesheet.test.csv`:
     ```csv
     sample,fastq_1,fastq_2,treatment
     sampleA,tests/data/sampleA.test_1.fastq.gz,tests/data/sampleA.test_2.fastq.gz,infected
     sampleB,tests/data/sampleB.test_1.fastq.gz,tests/data/sampleB.test_2.fastq.gz,infected
     ```
   - Keep all existing `*.bam` fixtures — module-level nf-tests for iVar/LoFreq/coverage/etc. continue to use them unchanged.
10. **nf-test specs (new):** one per new module (`stub` + real-fixture test; `IVAR_TRIM` uses `sampleA.test.bam` + `primers.test.bed`; `BWA_MEM` uses the new FASTQ fixtures + fixture reference). **Rewrite** `tests/subworkflows/input_check/main.nf.test` for the FASTQ samplesheet. Update `tests/e2e/main.nf.test` to the FASTQ ingress (assert `results_test/mouse_veev/Bwa/sampleA/sampleA.bam` or equivalent published path exists) — snapshot regeneration of the **e2e spec only** is authorized at this checkpoint (`nf-test test tests/e2e --update-snapshot`).
11. **`conf/modules.config`** — add `withName:` publishDir rules extending the G5 layout:
    ```
    ${params.outdir}/${params.dataset}/FastQC/<sample>/...
    ${params.outdir}/${params.dataset}/Bwa/<sample>/*.bam{,.bai}, *.bwa.log
    ${params.outdir}/${params.dataset}/Samtools_stats/<sample>/...
    ${params.outdir}/${params.dataset}/Ivar/<sample>/*.trimmed.sorted.bam{,.bai}   (IVAR_TRIM only)
    ```
12. **`tests/lint_structure.sh`** — extend the mandatory tree with every new file listed in §1 (Phase 4 rows).

### Task 4.5 — Software-provenance audit (retrofit, no rewrites)

1. Verify **every** file under `modules/local/**/main.nf` (including new Phase 4 modules) emits `path "versions.yml", emit: versions` from **both** `script:` and `stub:` blocks. Add only where missing.
2. Verify **every** subworkflow mixes its modules' versions into an `emit: versions` channel. Known gap to fix: `INPUT_CHECK` (add it, see 4.3.6). Audit `ANNOTATION`, `COVERAGE_QC`, `REPORTING` too.
3. Do **not** convert to `topic:` channels — CI pins Nextflow `24.04.4`; keep named emits + explicit `mix` (aggregation arrives in Phase 6).

### Phase 4 DoD — all must exit 0

| # | Checkpoint command | Pass criterion |
|---|---|---|
| 4.1 | `grep -q "def check_max" nextflow.config && ! grep -q "def check_max" conf/base.config && echo CHECKMAX-OK` | Single definition, in `nextflow.config` only |
| 4.2 | `grep -q "130..145" conf/base.config && grep -q "+ 104" conf/base.config && grep -q "task.attempt" conf/base.config && echo RETRY-OK` | Retry + dynamic scaling intact |
| 4.3 | `nextflow config -profile test >/dev/null && nextflow config -profile slurm >/dev/null && nextflow config -profile awsbatch >/dev/null && echo PARSE-OK` | All profiles parse with the `plugins` block present |
| 4.4 | `grep -q "nf-schema@2.1.1" nextflow.config && echo PLUGIN-OK` | Plugin pinned |
| 4.5 | `nextflow run . -profile test -stub-run` | Full DAG stub-runs clean (plugin downloads on first run) |
| 4.6 | `nextflow run . -profile test --protocol amplicon -stub-run 2>&1 \| grep -q "primer_bed is required"; test $? -eq 0 && echo CONDITIONAL-OK` | Missing `--primer_bed` under amplicon protocol **fails loudly** |
| 4.7 | `printf 'sample,fastq_1,fastq_2\nsampleX,tests/data/sampleA.test_1.fastq.gz,tests/data/sampleA.test_2.fastq.gz\n' > /tmp/bad_sheet.csv; nextflow run . -profile test --input /tmp/bad_sheet.csv -stub-run >/dev/null 2>&1; test $? -ne 0 && echo SCHEMA-FAIL-OK` | Samplesheet missing `treatment` is rejected by nf-schema |
| 4.8 | `nextflow run . -profile test --outdir results_test` | Real e2e on FASTQ fixtures, exit 0 |
| 4.9 | `test -f results_test/mouse_veev/LoFreq/sampleA/variants.filtered.vcf.gz && echo INGRESS-OK` | Downstream seam intact after FASTQ migration |
| 4.10 | `nf-test test tests/modules/local/fastqc tests/modules/local/bwa tests/modules/local/samtools/stats tests/modules/local/ivar --profile docker` | New + touched module specs PASS |
| 4.11 | `nf-test test --profile docker` | **Entire** suite green (pre-existing snapshots unchanged except authorized e2e update) |
| 4.12 | `bash tests/lint_containers.sh && bash tests/lint_structure.sh` | Both OK with new pins/files |
| 4.13 | Push branch; GitHub Actions | All jobs green |

---

## 3. PHASE 5 — Subworkflow Architecture (Feature Toggles)

**Goal:** Relocate the advanced-science branches into the flat `subworkflows/local/` convention as cleanly isolated, flag-gated units, and make the fan-in/fan-out channel logic explicit in `main.nf`. **Move, don't rewrite.**

### Task 5.1 — Relocate SELECTION and HAPLOTYPE

1. `git mv subworkflows/selection/main.nf subworkflows/local/selection.nf` and `git mv subworkflows/haplotype/main.nf subworkflows/local/haplotype.nf` (file moves preserve history; stage with `git add`, **no commit**).
2. The relative include depth is identical (`subworkflows/local/` → `../../modules/local/...`), so module `include` lines inside the moved files need **no changes**. Verify by parse (DoD 5.1).
3. Update `main.nf` includes:
   ```groovy
   include { SELECTION } from './subworkflows/local/selection'
   include { HAPLOTYPE } from './subworkflows/local/haplotype'
   ```
4. Keep the existing internal guards `if (params.run_snpgenie) { ... }` / `if (params.run_haplotype) { ... }` — they are the isolation mechanism. Optional branches must remain inert (no processes, empty emits) when flags are off.
5. Update the `script` path in `tests/subworkflows/selection/main.nf.test` and `tests/subworkflows/haplotype/main.nf.test` to point at the new locations; remove the now-empty old directories; update `tests/lint_structure.sh`.

### Task 5.2 — Explicit channel logic in `main.nf`

Replace the current inline `SELECTION(...)` call arguments with named, documented channels (keep everything else in `main.nf` untouched):

```groovy
// 5. Selection Analysis (optional, --run_snpgenie)
// Per-sample LoFreq VCFs (drop .tbi); aggregate fan-in (.collect()) lives inside SELECTION
ch_snpgenie_vcfs = VARIANT_CALLING.out.lofreq_vcf
    .map { meta, vcf, tbi -> [ meta, vcf ] }

// Treatment-grouped sample manifest for limma contrast validation.
// groupTuple() fans samples into one tuple per treatment group.
ch_treatment_groups = INPUT_CHECK.out.fastqs
    .map { meta, fastq_1, fastq_2 -> [ meta.treatment, meta.id ] }
    .groupTuple()
    .map { treatment, ids -> [ treatment: treatment, n_replicates: ids.size() ] }

// Hard guard: limma needs >= 2 treatments, each with >= 2 replicates
if (params.run_snpgenie) {
    ch_treatment_groups
        .filter { it.n_replicates < 2 }
        .subscribe { error "SNPGenie/limma requires >=2 replicates per treatment; group '${it.treatment}' has ${it.n_replicates}." }
}

SELECTION(
    ch_snpgenie_vcfs,
    INPUT_CHECK.out.fasta,
    INPUT_CHECK.out.gff,
    file(params.input),
    ch_treatment_groups          // new 5th arg: grouped manifest
)
```

1. Extend `subworkflows/local/selection.nf` `take:` with `ch_treatment_groups` and log the groups at startup (`ch_treatment_groups.subscribe { log.info "SELECTION treatment group: ${it.treatment} (n=${it.n_replicates})" }`). The legacy per-sample `.collect()` fan-in to `SNPGENIE_SUMMARIZE` stays as-is.
2. `HAPLOTYPE(...)` wiring is unchanged (`VARIANT_CALLING.out.viral_bams`, `INPUT_CHECK.out.fasta`); only its include path changes.
3. Confirm empty-channel safety: with flags off, `SELECTION.out.selection_tables` / `HAPLOTYPE.out.haplotypes` are empty and `REPORTING` still completes (existing behavior — verify, don't fix unless broken).

### Task 5.3 — Test & fixture alignment for the toggle matrix

1. Update `tests/e2e/main.nf.test` (or add a second e2e spec) to cover the **toggle matrix** at stub level: `--run_snpgenie false --run_haplotype false` excludes `SNPGENIE`/`CLIQUESNV`/`VILOCA` process names from the trace preview.
2. Ensure `tests/subworkflows/selection/main.nf.test` passes a `treatment`-column samplesheet consistent with the Phase 4 schema.

### Phase 5 DoD — all must exit 0

| # | Checkpoint command | Pass criterion |
|---|---|---|
| 5.1 | `test -f subworkflows/local/selection.nf && test -f subworkflows/local/haplotype.nf && test ! -d subworkflows/selection && test ! -d subworkflows/haplotype && echo LAYOUT-OK` | Flat layout, old dirs gone |
| 5.2 | `nextflow run . -profile test -stub-run` | DAG parses/stub-runs after relocation |
| 5.3 | `nextflow run . -profile test --run_snpgenie true --outdir results_snpgenie` | Selection branch runs e2e, exit 0 |
| 5.4 | `find results_snpgenie/mouse_veev/SNPGenie -name "*.tsv" \| head -1 \| xargs test -f && echo SNPGENIE-OUT-OK` | Legacy SNPGenie output layout intact |
| 5.5 | `nextflow run . -profile test --run_haplotype true --outdir results_haplo` | Haplotype branch runs e2e, exit 0 |
| 5.6 | `find results_haplo/mouse_veev/Haplotypes -type f \| head -1 \| xargs test -f && echo HAPLO-OUT-OK` | CliqueSNV/VILOCA outputs published |
| 5.7 | `nextflow run . -profile test --run_snpgenie false --run_haplotype false -stub-run 2>&1 \| grep -cE "SNPGENIE\|CLIQUESNV\|VILOCA"; test $? -eq 1 && echo TOGGLES-OFF-OK` | Zero optional-branch processes when flags off (`grep -c` exits 1 on no match) |
| 5.8 | `grep -q "groupTuple" main.nf && grep -q "collect" main.nf && echo CHANNEL-LOGIC-OK` | Explicit grouping/collection wiring present |
| 5.9 | `nf-test test tests/subworkflows tests/e2e --profile docker` | Subworkflow + e2e suites PASS |
| 5.10 | `nf-test test --profile docker` | Full suite green |
| 5.11 | Push branch; GitHub Actions | All jobs green |

---

## 4. PHASE 6 — MultiQC & Final DAG Assembly

**Goal:** Software-provenance aggregation, a terminal MultiQC dashboard, final DAG plumbing, docs, and the **`v1.1.0`** release cut.

### Task 6.1 — `modules/local/dumpsoftwareversions/main.nf`

1. New module, label `process_low`, `container "${params.container_multiqc}"` (MultiQC image ships Python + PyYAML — reuse it rather than adding a pin).
2. Input: `path versions` — the collated set of all per-process `versions.yml` files (staged into `versions/` via `stageAs: 'versions/*'`).
3. Script: Python heredoc that glob-reads `versions/*.yml`, merges into one nested map keyed by process name, and writes:
   - `software_versions.yml` — canonical provenance artifact;
   - `software_versions_mqc.yml` — the same data wrapped as MultiQC custom content (`# plot_type: 'html'` section header format used by nf-core `custom/dumpsoftwareversions`).
4. Emits `path "software_versions.yml", emit: yml`, `path "software_versions_mqc.yml", emit: mqc_yml`, plus its own `versions.yml`. Publish to `${params.outdir}/pipeline_info/` via `conf/modules.config`.
5. **`main.nf` plumbing** (this closes the Phase 4 provenance loop):
   ```groovy
   ch_versions = Channel.empty()
       .mix(INPUT_CHECK.out.versions)
       .mix(READ_PREPROCESSING.out.versions)
       .mix(VARIANT_CALLING.out.versions)
       .mix(ANNOTATION.out.versions)
       .mix(COVERAGE_QC.out.versions)
       .mix(SELECTION.out.versions)
       .mix(HAPLOTYPE.out.versions)
       .mix(REPORTING.out.versions)

   DUMP_SOFTWARE_VERSIONS(ch_versions.unique().collect())
   ```

### Task 6.2 — `modules/local/multiqc/main.nf` + config

1. Resolve and pin `container_multiqc` in `conf/containers.config` (MultiQC ≥ 1.21, e.g. `1.25.x`; G2 `docker manifest inspect` verification + digest comment mandatory).
2. `assets/multiqc_config.yml` (new): `title`, `subtitle`, `report_comment`, and `extra_fn_clean_exts` entries that strip `.test_1`, `.trimmed.sorted`, `.viral_only` etc. so report sample names match `meta.id`.
3. New module: label `process_low`; input `path multiqc_files` + `path multiqc_config`; script `multiqc --force --config ${multiqc_config} .`; emits `path "*multiqc_report.html", emit: report`, `path "*_data", emit: data`, `versions.yml`. Publish to `${params.outdir}/${params.dataset}/MultiQC/`.

### Task 6.3 — Terminal wiring in `main.nf`

```groovy
// MultiQC file aggregation — collect() fans every QC artifact into one invocation
ch_multiqc_files = Channel.empty()
    .mix(READ_PREPROCESSING.out.fastqc_zip.map { meta, zip  -> zip })
    .mix(READ_PREPROCESSING.out.bwa_log.map    { meta, log  -> log })
    .mix(READ_PREPROCESSING.out.stats.map      { meta, stats, flagstat, idxstats -> [ stats, flagstat, idxstats ] })
    .mix(VARIANT_CALLING.out.lofreq_qc.map     { meta, qc   -> qc })
    .mix(DUMP_SOFTWARE_VERSIONS.out.mqc_yml)
    .collect()

MULTIQC(ch_multiqc_files, file("${projectDir}/assets/multiqc_config.yml"))
```

- Coverage: **FastQC** (zips), **BWA** (alignment logs staged in the search dir; machine metrics come from), **Samtools** (`stats`/`flagstat`/`idxstats` — parsed natively), **callers** (LoFreq `qc_stats.txt` staged for the report directory), **software versions** (custom content). MultiQC parses what it has modules for; all staged files ship inside the report data dir regardless.
- `MULTIQC` must be the **last** invocation in `workflow { }`.

### Task 6.4 — Tests, CI, docs

1. nf-test specs for both new modules (stub + real; `MULTIQC` real test asserts `multiqc_report.html` exists and non-empty; `DUMP_SOFTWARE_VERSIONS` asserts the yml contains `lofreq:` and `ivar:` keys). Add to the CI `nftest` matrix automatically (they live under `tests/modules/**`).
2. Extend `tests/e2e/main.nf.test`: assert `MultiQC/multiqc_report.html` exists and `pipeline_info/software_versions.yml` lists every tool (spot-grep `samtools`, `lofreq`, `ivar`, `multiqc`).
3. Add a CI job `e2e-full` to `.github/workflows/ci.yml` (ubuntu-latest, docker available): checkout → setup-nextflow → `chmod +x bin/*` → `nextflow run . -profile test --outdir results_ci` → assert `results_ci/mouse_veev/MultiQC/multiqc_report.html` exists. Also add the amplicon smoke: `nextflow run . -profile test --protocol amplicon --primer_bed tests/data/primers.test.bed -stub-run`.
4. **README.md**: update samplesheet table to `sample,fastq_1,fastq_2,treatment`; document `--protocol`/`--primer_bed` with an amplicon example; add MultiQC to the output-layout section and DAG diagram; note nf-schema validation behavior.
5. **`AGENT.md`**: add a one-line banner at the top — *"Phases 4–6 of the production roadmap are specified in `AGENT_UPDATE.md`, which supersedes this file where they conflict."*
6. **`CHANGELOG.md`**: new `## [1.1.0] - <date>` section listing: nf-schema validation, FASTQ ingress (FastQC/BWA/iVar trim), `check_max` relocation, provenance aggregation, MultiQC, subworkflow relocation.

### Task 6.5 — Release `v1.1.0` (deviation from the original "`v1.0.0`" instruction)

1. `v1.0.0` **already exists** as a tag (pre-hardening cut). Do **not** move, delete, or re-tag it. Record in `DEVIATIONS.md`: *"Phase 6 release cut as `v1.1.0` because `v1.0.0` was tagged before the hardening phases."*
2. Bump `manifest.version` in `nextflow.config` to `'1.1.0'`.
3. Prepare — but **do not execute without explicit user confirmation** — and present verbatim:
   ```bash
   git add -A && git commit -m "chore(release): 1.1.0"
   git tag -a v1.1.0 -m "v1.1.0: production hardening (nf-schema, resiliency), FASTQ ingress, MultiQC, provenance aggregation"
   git push origin main --tags
   ```
4. After approved execution, verify `git describe --tags` prints `v1.1.0`.

### Phase 6 DoD — all must exit 0

| # | Checkpoint command | Pass criterion |
|---|---|---|
| 6.1 | `nf-test test tests/modules/local/dumpsoftwareversions tests/modules/local/multiqc --profile docker` | New module specs PASS |
| 6.2 | `nextflow run . -profile test --outdir results_v11` | Full DAG incl. terminal MULTIQC, exit 0 |
| 6.3 | `test -f results_v11/mouse_veev/MultiQC/multiqc_report.html && echo MQC-OK` | Dashboard published |
| 6.4 | `for t in samtools lofreq ivar multiqc fastqc bwa; do grep -qi "$t" results_v11/pipeline_info/software_versions.yml \|\| { echo "MISSING $t"; exit 1; }; done; echo VERSIONS-OK` | Provenance artifact complete |
| 6.5 | `nextflow run . -profile test --protocol amplicon --primer_bed tests/data/primers.test.bed --outdir results_amp` | Amplicon branch e2e (incl. IVAR_TRIM), exit 0 |
| 6.6 | `nextflow run . -profile test --outdir results_v11 -resume 2>&1 \| grep -q "cached" && echo RESUME-OK` | Resume still works post-plumbing |
| 6.7 | `nf-test test --profile docker` | Entire suite green |
| 6.8 | `tmp=$(mktemp -d); git archive HEAD \| tar -x -C "$tmp"; (cd "$tmp" && nextflow run . -profile test --outdir results_clean) && echo CLEAN-CLONE-OK` | Pristine-export reproducibility |
| 6.9 | `bash tests/lint_containers.sh && bash tests/lint_structure.sh` | Final lints green |
| 6.10 | `grep -q "version = '1.1.0'" nextflow.config && echo VERSION-BUMP-OK` | Release version set |
| 6.11 | `git describe --tags` (after user-approved tagging) | `v1.1.0` |
| 6.12 | CI on `main` post-merge | All jobs green incl. new `e2e-full` |

---

## 5. PHASE 7 — Executive HTML Reporting & Public "Tiny" Test Dataset (UX Polish)

> **Guiding principle (binding):** The biological stack is **frozen**. LoFreq + iVar + SNPGenie + CliqueSNV/VILOCA is already state-of-the-art for this use case — do **not** add another variant caller, assembler, or haplotype tool to compete with it. This phase exclusively ports **engineering and user-experience** patterns from nf-core/viralrecon and V-pipe: self-contained HTML reporting and zero-friction test data.

**Goal:** (1) Wrap the existing report artifacts into one self-contained HTML executive summary per run. (2) Generalize the pipeline identity from alphavirus-specific to virus-agnostic, with an honest compatibility note. (3) Let anyone run the pipeline immediately after cloning — or straight from GitHub — against a public minimal dataset.

### Task 7.1 — Automated executive HTML report (Quarto)

The existing `modules/local/report/*` scripts keep emitting their PNG/TSV/Excel artifacts **untouched**. A new terminal reporting step wraps them into a single `.html` executive summary. Positioning vs. MultiQC: **MultiQC = QC dashboard** (was the data good?); **executive report = findings** (what did we find?).

1. **`assets/executive_report.qmd`** (new):
   - YAML header: `format: html` with `embed-resources: true` (one self-contained file, figures base64-embedded), plus a `params:` block (`dataset`, `has_selection`, `has_haplotype`).
   - Engine: **Jupyter / python3** (pandas) — deliberately avoids an R/knitr runtime dependency; all heavy computation already happened upstream in the existing modules.
   - Sections: run-summary KPI table (from `consolidated_sample_summary.tsv`), sliding-window variant plots, coverage plots, selection key tables (**conditional** on `params.has_selection`), haplotype plots (**conditional** on `params.has_haplotype`), software-versions table (from `software_versions.yml`).
   - Conditional sections are driven by `quarto render -P key:value` parameter overrides — never by editing the template at runtime.
2. **Container resolution (G2 decision tree, in order):** (a) check `https://quay.io/repository/biocontainers/quarto?tab=tags` for a pinned quarto image with python; (b) try a `mulled-v2` combination via the bioconda mulled index; (c) **fallback** — custom image `containers/executive_report/Dockerfile`: pinned base image, pinned Quarto release tarball (URL + sha256 verification inside the Dockerfile), pinned `pandas`/`jupyter`/`ipykernel`; built and pushed to `ghcr.io/aleponce4/alphavirus-variant-analysis-workflow/executive_report:1.0.0` by a new `.github/workflows/containers.yml`, pinned as `params.container_executive_report` in `conf/containers.config`, and recorded in `DEVIATIONS.md` per G2. Whichever path resolves, verify with `docker manifest inspect`. Use the **new pipeline slug from Task 7.2** in the ghcr path (`ghcr.io/aleponce4/viral-intrahost-variant-analysis/executive_report:1.0.0`) — if Task 7.2 has not landed yet, leave a `TODO(7.2)` and finalize the pin there.
3. **`modules/local/report/executive_report/main.nf`** — process `EXECUTIVE_REPORT`, label `process_low`:
   - Inputs (all staged into one work dir): run-summary TSV, variant-plot PNGs, coverage-plot PNGs, haplotype-plot PNGs (may be empty), selection key table (may be empty), `software_versions.yml`, and the qmd template.
   - Script computes `-P has_selection:true/false` and `-P has_haplotype:true/false` from staged-file presence, then `quarto render executive_report.qmd --output executive_report.html -P dataset:${params.dataset} ...`.
   - Emits `path("executive_report.html")` + `versions.yml` (quarto version). Publish to `${params.outdir}/${params.dataset}/Reports/` via `conf/modules.config`.
4. **Wiring** — minimal diff in `subworkflows/reporting/main.nf`: `mix` + `.collect()` the upstream report artifacts, with `.ifEmpty([])` guards on the optional-branch channels so the report still renders under `--run_snpgenie false` / `--run_haplotype false`; pass `DUMP_SOFTWARE_VERSIONS.out.yml` down from `main.nf`. `EXECUTIVE_REPORT` is the **last step of `REPORTING`**.
5. **nf-test:** stub test + real test asserting `executive_report.html` exists, size > 100 KB (proves resources are embedded), and contains the marker strings `Executive Summary` and the dataset id; extend `tests/e2e/main.nf.test` with the same published-path assertion; one additional run with both optional branches off asserting exit 0 (graceful degradation).

### Task 7.2 — Scope generalization: rename + compatibility note

The workflow is **alphavirus-validated but virus-agnostic by design**: the reference FASTA/GFF3 and target contig are user-supplied params, and nothing in the DAG is alphavirus-specific. Rebrand accordingly so users don't assume an alphavirus-only tool.

1. **New identity:** working slug **`viral-intrahost-variant-analysis`** — confirm the final slug with the user before executing (it appears in public URLs; the GitHub repo rename itself is a user action in the GitHub UI, after which GitHub auto-redirects the old URL — all in-repo references must still be updated). Update:
   - `nextflow.config` → `manifest.name`, `manifest.description`, `manifest.homePage`.
   - `main.nf` → startup banner (replace the `A L P H A V I R U S   V A R I A N T   A N A L Y S I S` block with the new name).
   - `nextflow_schema.json` and `assets/schema_input.json` → `title` and `$id`.
   - `README.md` → title, badges, quick-start slug; new **"Scope & Limitations"** section (item 2 below).
   - `CHANGELOG.md` → under `[1.2.0]`: *"Renamed from `alphavirus-variant-analysis` to `viral-intrahost-variant-analysis`; first release under the generalized identity."*
   - All Phase 7 references going forward (ghcr path in Task 7.1, test-datasets repo in Task 7.3, DoD URLs).
   - **Do NOT rename** `params.dataset = 'mouse_veev'` or fixture/file names — the test data *is* legitimately alphavirus-derived; that stays accurate.
2. **Compatibility note (binding text)** — add verbatim to the new README "Scope & Limitations" section, and as a comment block above `manifest` in `nextflow.config`:
   > Validated on alphaviruses (VEEV/EEEV). The workflow is reference-driven and is expected to work for other viruses given an appropriate reference FASTA + GFF3. **Not tested for highly recombinant or hypervariable viruses** (e.g., HIV, HCV): single-reference alignment, haplotype reconstruction (CliqueSNV/VILOCA), and SNPGenie selection statistics all assume a representative reference and low within-host recombination. For such datasets, interpret `--run_haplotype` / `--run_snpgenie` outputs with caution or leave those branches disabled.

### Task 7.3 — Public "tiny" test dataset (clone-and-run demo)

Mirror `nextflow run nf-core/viralrecon -profile test`: zero local data required.

1. **Generate deterministically (documented):** extend `tests/data/generate_fixtures.sh` with an optional `--export-tiny <dir>` mode that down-samples the existing synthetic FASTQ fixtures to **~25,000 read pairs (~50k reads total)** via `seqtk sample -s 42` (pinned `quay.io/biocontainers/seqtk` container; resolve per G2), and copies the fixture reference FASTA/GFF3 plus a `samplesheet.tiny.csv` beside them. **Validation gate before publishing:** run the pipeline once locally on the tiny dataset and confirm ≥ 1 known seeded fixture variant appears in the LoFreq VCF.
2. **Host publicly:** create the public repo `aleponce4/viral-intrahost-test-datasets` (requires user confirmation — it is a new public artifact; named after the Task 7.2 slug for discoverability — its README must state the synthetic data is **alphavirus-derived**, which is fine and honest: the demo data is alphavirus, the tool is general). Commit the tiny dataset + `checksums.sha256` + a README with the exact regeneration commands, and record the **commit SHA**. All pipeline references must use **SHA-pinned raw URLs** — `https://raw.githubusercontent.com/aleponce4/viral-intrahost-test-datasets/<SHA>/...` — never branch names (immutability).
3. **`conf/test.config`:** add `params.test_data_base = 'https://raw.githubusercontent.com/aleponce4/viral-intrahost-test-datasets/<SHA>'` and point `params.input = "${params.test_data_base}/samplesheet.tiny.csv"`, `params.fasta`, `params.gff` at the hosted files. The hosted samplesheet lists the remote FASTQ URLs (the Phase 4 nf-schema pattern `^\S+\.f(ast)?q(\.gz)?$` already accepts them).
4. **Keep hermetic offline tests:** the committed `assets/samplesheet.test.csv` and `tests/data/*` fixtures **stay**. `tests/nextflow.config` (nf-test harness) overrides `params.input`/`params.fasta`/`params.gff` back to the local files so `nf-test test --profile docker` never depends on the network. Rule: `-profile test` = public demo; nf-test = hermetic CI. Document this split in the README.
5. **README quick-start** becomes:
   ```bash
   nextflow run aleponce4/viral-intrahost-variant-analysis -profile test --outdir results_demo
   ```
   (plus the clone-then-run form), with a Docker requirement note and expected runtime on the tiny dataset.

### Task 7.4 — Release `v1.2.0`

1. `CHANGELOG.md`: `## [1.2.0] - <date>` — generalized identity/rename (Task 7.2), executive HTML report, public tiny dataset, `containers.yml` (if the custom-image fallback was used).
2. `manifest.version` → `'1.2.0'`; prepare — but **do not execute without explicit user confirmation** — and present verbatim:
   ```bash
   git add -A && git commit -m "chore(release): 1.2.0"
   git tag -a v1.2.0 -m "v1.2.0: generalized virus-intrahost identity, executive HTML reporting, public tiny test dataset"
   git push origin main --tags
   ```
3. After approved execution, verify `git describe --tags` prints `v1.2.0`.

### Phase 7 DoD — all must exit 0

| # | Checkpoint command | Pass criterion |
|---|---|---|
| 7.1 | `nextflow run . -profile test --outdir results_v12` | Full run incl. `EXECUTIVE_REPORT`, exit 0 |
| 7.2 | `test -f results_v12/mouse_veev/Reports/executive_report.html && test $(stat -c%s results_v12/mouse_veev/Reports/executive_report.html) -gt 102400 && grep -q "Executive Summary" results_v12/mouse_veev/Reports/executive_report.html && echo REPORT-OK` | Self-contained HTML rendered |
| 7.3 | `nextflow run . -profile test --run_snpgenie false --run_haplotype false --outdir results_v12_min` | Report still renders with optional branches off, exit 0 |
| 7.4 | `nf-test test tests/modules/local/report tests/e2e --profile docker` | Report module + e2e specs PASS |
| 7.5 | `grep -q "name = 'viral-intrahost-variant-analysis'" nextflow.config && grep -qi "viral.intrahost" main.nf && echo RENAME-OK` | Manifest + banner carry the new identity |
| 7.6 | `grep -qi "hypervariable" README.md && grep -qi "alphavirus" README.md && echo SCOPE-NOTE-OK` | Compatibility note present (validated on alphaviruses; not tested for highly recombinant/hypervariable viruses) |
| 7.7 | `curl -sI "https://raw.githubusercontent.com/aleponce4/viral-intrahost-test-datasets/<SHA>/samplesheet.tiny.csv" \| head -1 \| grep -q "200" && echo REMOTE-OK` | Tiny dataset publicly reachable at SHA-pinned URL |
| 7.8 | `tmp=$(mktemp -d); cd "$tmp" && nextflow run aleponce4/viral-intrahost-variant-analysis -profile test --outdir results_demo` | Zero-local-files public demo run (requires the user-confirmed repo rename), exit 0 |
| 7.9 | `zgrep -q "KP282671.1" results_demo/mouse_veev/LoFreq/sampleA/variants.filtered.vcf.gz && echo TINY-VARIANT-OK` | Seeded variant recovered from the tiny dataset |
| 7.10 | `nf-test test --profile docker` | Entire suite green (offline/hermetic — no network dependency) |
| 7.11 | `bash tests/lint_containers.sh && bash tests/lint_structure.sh` | Lints green (including approved `ghcr.io` fallback if used) |
| 7.12 | `grep -q "version = '1.2.0'" nextflow.config && echo VERSION-OK`; after approved tagging: `git describe --tags` | Version set; tag prints `v1.2.0` |
| 7.13 | CI on `main` post-merge | All jobs green |

---

## 6. Global Definition of Done (project level, after Phase 7)

- [ ] All Phase 4/5/6 DoD tables passed, in order.
- [ ] `nf-test test --profile docker` green from a pristine `git archive` export.
- [ ] `check_max()` defined exactly once (in `nextflow.config`); retry covers exit 130–145 + 104; all resources scale with `task.attempt`.
- [ ] `nf-schema@2.1.1` active: invalid params and malformed samplesheets fail fast with actionable messages; `--protocol amplicon` without `--primer_bed` fails.
- [ ] Every process emits `versions.yml`; `pipeline_info/software_versions.yml` aggregates all of them; `MultiQC/multiqc_report.html` is produced on every full run.
- [ ] `--run_snpgenie` / `--run_haplotype` branches are fully inert when off and fully functional when on.
- [ ] Zero `latest` tags / foreign registries (`tests/lint_containers.sh`); zero absolute paths outside configs.
- [ ] `DEVIATIONS.md` records the `v1.1.0`-instead-of-`v1.0.0` release decision.
- [ ] `v1.1.0` tag created **only after explicit user confirmation**.
- [ ] All Phase 7 DoD checkpoints passed.
- [ ] `executive_report.html` renders self-contained, with and without the optional branches.
- [ ] `-profile test` runs end-to-end using only the SHA-pinned public tiny dataset; nf-test remains fully offline/hermetic.
- [ ] Biological stack unchanged (no new callers/assemblers/haplotype tools) — Phase 7 changed reporting, data hosting, identity, and config only.
- [ ] Pipeline renamed to the generalized intrahost-virus identity everywhere (manifest, banner, schema, README, public URLs); README "Scope & Limitations" carries the recombinant/hypervariable caveat.
- [ ] `v1.2.0` tag created **only after explicit user confirmation**.

## 7. Standing prohibitions (unchanged, repeated)

- No DSL1 syntax. No `conda`. No local tool assumptions. No per-process resource literals (labels only). No absolute paths in `main.nf`, `modules/`, `subworkflows/`, `bin/`.
- No git mutations without explicit, at-the-moment user confirmation (applies especially to Tasks 5.1 staging and 6.5 tagging).
- No rewriting working modules/subworkflows; retrofit per the minimal-diff rule.
- No skipping a DoD checkpoint. If blocked, stop and report.

---

*End of update specification. Executor: begin at Phase 4, Task 4.1.*
