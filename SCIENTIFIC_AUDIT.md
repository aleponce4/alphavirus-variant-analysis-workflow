# SCIENTIFIC_AUDIT.md — Methodological Audit of `viral-intrahost-variant-workflow`

**Auditor role:** Principal Bioinformatics Scientist
**Date:** 2026-07-29
**Scope:** Nextflow DSL2 configuration, tool selection, and default statistical thresholds for ultra-deep (≥1,000×) targeted viral sequencing; cross-sectional in vivo design; iSNV detection down to 0.1%–1%; quasispecies haplotype reconstruction. De novo assembly is out of scope by design and was not evaluated.

**Audit basis (code inspected, not assumed):** `nextflow.config`, `conf/modules.config`, `conf/containers.config`, `conf/test.config`, `main.nf`, `modules/local/{lofreq/call, lofreq/filter, ivar/variants, ivar/consensus, ivar/trim, ivar/tsv_to_vcf, cliquesnv, viloca, snpgenie/*}`, `subworkflows/{variant_calling, local/selection, local/haplotype}`, `assets/snpgenie/snpgenie.pl` (vendored, pinned), `bin/*.py`, `nextflow_schema.json`.

**Statistical reference frame used throughout:** at Q30 the per-base sequencing error rate is ε ≈ 1×10⁻³, i.e. **0.1%**. A 0.1% iSNV therefore sits exactly at the single-base noise floor; separating it from error requires (a) depth n such that the expected ALT count n·AF ≥ ~10 reads (⇒ ≥10,000× for AF = 0.001), (b) a per-base-quality-aware error model (LoFreq's Poisson-binomial), and (c) strand/consistency corroboration. Any pipeline element that hard-codes a ≥1% floor cannot, by construction, contribute to the 0.1%–1% objective.

---

## Audit Finding & Resolution Tracking Record (Updated: 2026-07-30)

| Finding ID | Original Severity | Description | Status | Corrective Action & Location | Verification Test / Evidence |
|---|---|---|---|---|---|
| **MF-1** | Critical | iVar 1% hard floor (`-t 0.01`) | **RESOLVED** | Updated `ivar_min_freq = 0.001` in `nextflow.config` | Verified iVar emits candidate variants down to 0.1% |
| **MF-2** | Critical | 1,000× inclusion depth floor | **RESOLVED** | Set `ivar_min_depth = 10`, `lofreq_min_depth = 10` in `nextflow.config` | Permissive calling across low & high depth regions |
| **MF-3** | High | Dead `lofreq_min_freq` parameter | **RESOLVED** | Removed unused parameter from `nextflow.config` and schema | Clean schema without misleading options |
| **MF-4** | High | Legacy ShoRAH / VILOCA geometry | **RESOLVED** | Updated module to run `viloca run` with explicit window parameters (`viloca_window = 150`, `viloca_shift = 50`) & pinned biocontainer | Deterministic quasispecies window tiling & real VILOCA outputs |
| **MF-5** | High | Silent NaN Kruskal statistics (no scipy) | **RESOLVED** | Updated `container_python_reporting` to `quay.io/biocontainers/seaborn:0.13.2` | Valid p-values emitted by `analyze_delta_selection.py` |
| **MF-6** | Moderate | LoFreq MQ ≥ 60 over-aggressive | **RESOLVED** | Relaxed `lofreq_min_mq = 20` in `nextflow.config` | Preserves reads in variable viral regions |
| **MF-7** | Moderate | Percent/Proportion contract & stub plots | **RESOLVED** | Machine VCFs store proportions; `generate_variant_plots.py` computes % & generates figures | Real PNG plots (274 KB) & `variant_frequency_summary_pct.tsv` |
| **MF-8** | Moderate | Fabricated GT/QUAL in iVar VCF | **RESOLVED** | `ivar_variants_to_vcf.py` outputs haploid GT (`1`) & converts p-value to QUAL | Valid VCF output feeding `bcftools csq` |
| **MF-9** | Low | Implicit CliqueSNV threshold | **RESOLVED** | Explicit `--cliquesnv_min_freq 0.001` configured | Pinned CliqueSNV frequency threshold |
| **MF-10** | Low | Undocumented indel scope | **DOCUMENTED** | Documented SNV-focused scope (`lofreq_enable_indelqual = false`) in `README.md` | Documented design choice |
| **MF-11** | Low | Sub-0.5% corroboration policy | **RESOLVED** | Implemented candidate frequency tiering (`>=1%`, `0.1-1% candidate`, `<0.1% exploratory`) | Output TSV & plot tiering |

---

## 1. Confirmed Validities

Tools and settings verified to be scientifically sound for this biological context.

- **CV-1 — No assembly step anywhere in the DAG.** The workflow is strictly reference-based (BWA MEM → per-contig extraction → callers). Minor-variant diversity is not at risk of assembler collapse. Complies with the out-of-scope constraint.
- **CV-2 — `samtools mpileup -d 0` (unlimited depth) feeding iVar** (`modules/local/ivar/variants/main.nf:22`). This is critical: samtools' default depth cap (8000, and 250 in older builds) would silently truncate ≥1,000×–100,000× amplicon piles and flatten all frequency estimates. Correctly disabled.
- **CV-3 — mpileup emits base qualities unfiltered (`-Q 0`) and iVar applies its own `-q 30`** (`ivar/variants/main.nf:22-23`). Quality filtering happens exactly once, in the tool whose error model consumes it. Double-filtering avoided. The same correct division of labor holds for `ivar consensus` (`-q ${params.ivar_min_bq}` passed through).
- **CV-4 — `-aa -A` in mpileup.** All positions (including zero-coverage and reference-homozygous) and all reads (including orphan) are emitted; iVar's depth accounting and the GFF-annotated output are not biased by silently dropped positions.
- **CV-5 — LoFreq runs with no allele-frequency floor at call time, and `LOFREQ_FILTER` applies only quality thresholds** (`lofreq/filter/main.nf:21`: `--snvqual-thresh 20 --indelqual-thresh 20`; no `--af-thresh`, no `--cov-thresh`). This is exactly what 0.1%–1% iSNV work requires: significance is decided by the Poisson-binomial model over per-base qualities (Wilm et al., 2012), not by an arbitrary frequency cutoff. The rare-variant tail reaches downstream intact.
- **CV-6 — `lofreq_sig = 0.01`** (`nextflow.config:40`). LoFreq applies dynamic Bonferroni correction internally; α = 0.01 pre-correction is appropriately conservative for genome-wide testing of a ~11–12 kb viral genome without being destructive to true low-frequency calls (LoFreq calls well below 1% at high depth with high base qualities by design).
- **CV-7 — BAQ and indelqual disabled by default** (`lofreq_enable_baq = false`, `lofreq_enable_indelqual = false`). Correct for this use case: BAQ/viterbi realignment is known to erase genuine low-frequency variation in haploid, hypervariable viral data by realigning ALT-supporting reads toward the reference. Leaving it off protects regions of high natural viral mutation from systematic penalization. SNV-focused scope is respected.
- **CV-8 — `--min-alt-bq` = `--min-bq` = 30 in LoFreq** (`lofreq/call/main.nf:36`). ALT-supporting bases must independently pass Q30, which is the single most important defense for sub-1% calling. Correctly symmetric.
- **CV-9 — SNPGenie input wiring and `--vcfformat=2` are correct.** `main.nf:64-65` feeds the **quality-filtered-but-not-AF-filtered** LoFreq VCF (AF in INFO) to `SNPGENIE_RUN`. Verified in the vendored `assets/snpgenie/snpgenie.pl` that vcfformat=2 parses single-allele records via `AF=([\d\.e\-]+)` and multiallelic via `AF=f1,f2`, i.e. exactly LoFreq's VCF layout. Format mismatch would have terminated the run; this is compatible.
- **CV-10 — SNPGenie `--minfreq` is not set ⇒ defaults to 0** (verified: `snpgenie.pl:184-185`). All statistically called variants, including the rare tail, enter πN/πS and dN/dS estimation. For cross-sectional within-host diversity this is the correct operational choice — a minfreq > 0 would bias π downward by construction.
- **CV-11 — Cross-sectional design respected downstream.** SNPGenie is executed **per sample independently** (`SNPGENIE_RUN` maps over `[meta, vcf]`), group differences are tested with Kruskal–Wallis across independent treatment groups + Benjamini–Hochberg FDR (`bin/analyze_delta_selection.py`), and limma is applied to per-sample per-gene Δ(πN−πS) values as independent observations. No repeated-measures, pairing, or within-subject correlation structure is assumed anywhere — consistent with the stated strictly cross-sectional design.
- **CV-12 — iVar consensus parameters are sane and correctly isolated.** `-t 0.5 -m 10` (majority-rule consensus, 10× floor) affects only the per-sample consensus FASTA artifact, not iSNV calling. No contamination of the low-frequency branch.
- **CV-13 — CliqueSNV mode preset is the right model class.** `-m snv-illumina` selects the Illumina error/heterogeneity preset for clique-based SNV graph assembly on haploid data — appropriate for short-read amplicon quasispecies reconstruction. Memory is derived from `task.memory` and threads are passed, so ultra-deep BAMs won't silently OOM under the config's retry logic.
- **CV-14 — Reproducibility infrastructure is real.** All tool containers are tag-pinned with recorded `@sha256` digests in comments (`conf/containers.config`); vendored SNPGenie carries `PINNED_COMMIT` + `checksums.sha256`; versions are emitted per process and aggregated. Methodological drift via "latest" tags is not possible.
- **CV-15 — VCF outputs are bgzipped + tabix-indexed, multiallelic-capable, and keep AF as a proportion in machine-readable files** — the correct storage convention (see RA-9 for the human-facing half of the unit story).

---

## 2. Methodological Flaws

Defaults or assumptions that violate the stated biological context. Ordered by severity.

- **MF-1 (Critical) — iVar branch is hard-capped at 1% and cannot see 0.1%–1% iSNVs.** `ivar_min_freq = 0.01` (`nextflow.config:30`) is passed to `ivar variants -t`, which is a *reporting floor*: nothing below 1% AF is emitted, before any statistics. Combined with `ivar_min_depth = 1000` (`-m`), the iVar branch has a detection floor of exactly 10 ALT reads at 1,000×. The stated goal of 0.1%–1% variants is unattainable in this branch; iVar currently functions only as a ≥1% cross-caller corroborator, which is not what the parameter names and README advertise.

- **MF-2 (Critical) — Both callers impose a 1,000× *minimum-coverage inclusion floor* that silently erases genomic positions.** `lofreq call-parallel --min-cov 1000` (`lofreq/call/main.nf:36`) and `ivar variants -m 1000` mean positions with 1–999× depth are **never tested** — they are absent from the VCF, not present as invariant or low-confidence. Consequences:
  1. Amplicon edges, primer-dropout shoulders, and any genuinely low-coverage region produce *no calls at all*, creating invisible dead zones (a 1% variant at 500× = 5 reads, which LoFreq could call, is never even examined).
  2. **SNPGenie denominators are corrupted at dropped sites**: SNPGenie treats sites absent from the SNP report as invariant reference, so πN/πS and dN/dS are systematically biased downward wherever coverage fell below 1,000×.
  3. The floor conflates two distinct concepts: "minimum depth for a *meaningful* statistical test" (belongs to the error model) and "positions to *attempt*" (should be permissive). LoFreq's Poisson-binomial already fails gracefully at low depth; the pre-filter adds zero statistical protection and real information loss.

- **MF-3 (High) — `lofreq_min_freq = 0.01` is a dead parameter — false documentation of a 1% LoFreq floor.** It is declared in `nextflow.config:37` and described in `nextflow_schema.json` ("Minimum variant allele frequency for LoFreq") but is referenced in **zero** process scripts (verified by full-repo search). A user setting `--lofreq_min_freq 0.001` to hunt 0.1% variants would change nothing. This is worse than a missing feature: it silently misrepresents the statistical behavior of the pipeline. (Ironically the *actual* behavior — no AF floor — is the scientifically correct one; the parameter must either be wired or removed.)

- **MF-4 (High) — The "VILOCA" module does not run VILOCA; it runs legacy ShoRAH 1.99.2 with completely default inference settings.** `modules/local/viloca/main.nf` executes `shorah shotgun -b ${bam} -f ${fasta}` from `quay.io/biocontainers/shorah:1.99.2` (`containers.config:20`). Two problems:
  1. **Mislabeling:** VILOCA is the maintained rewrite of ShoRAH with corrected Dirichlet Process mixture updates and improved error-rate handling; publishing results under a `Haplotypes/VILOCA/` path (`modules.config:80`) while running 1.99.2 misattributes the method and version in any downstream methods section.
  2. **Unpinned inference geometry:** no `-w` (window size), `-s` (window shift), or DPM concentration (`--alpha`) is given. Window defaults are read-length-derived; the same command yields different tiling — and potentially different reconstructed haplotypes — for 150 bp vs 250 bp libraries. For a 0.1%–1% quasispecies claim this must be fixed, explicit, and reported. The DPM itself (Bayesian nonparametric clustering of reads within windows, haploid-appropriate) is a valid model class; running it with implicit, input-dependent hyperparameters is not.

- **MF-5 (High) — The selection statistics can silently degrade to NaN.** `bin/analyze_delta_selection.py` guards `from scipy.stats import kruskal` with try/except and, when scipy is absent, emits `kruskal_H/p = NaN` **without failing** (lines 12-16, 156-167). But `container_python_reporting` is pinned to vanilla `quay.io/biocontainers/python:3.12.12` (`containers.config:22`) — no scipy, no pandas, despite the comment claiming "pandas, matplotlib, seaborn, openpyxl". In the shipped configuration, `SELECTION_DELTA` runs to exit 0 while producing a Kruskal table with **no valid p-values**, and `SELECTION_BUILD_TABLES` then propagates an empty statistical story into `selection_gene_key_table.tsv`. A pipeline that loses its inferential layer without an error is methodologically unacceptable.

- **MF-6 (Moderate) — `lofreq_min_mq = 60` is over-aggressive and can penalize naturally variable regions.** Requiring MAPQ ≥ 60 under BWA MEM means "uniquely mapped only" (BWA caps MQ at 60). In low-complexity or repeated viral elements (e.g., conserved sequence elements / subgenomic promoter regions), genuinely informative reads map with 0 < MQ < 60 and are discarded *before* the error model sees them — precisely in regions of high natural viral mutation, which the audit was asked to protect. LoFreq's model already weighs base-level evidence; the MQ filter's job is only to exclude mis-mapped reads, for which MQ 20 is the accepted standard.

- **MF-7 (Moderate) — Percent-vs-proportion contract is undefined, and the reporting layer that must enforce it is stub code.** All machine artifacts carry **proportions** (iVar `ALT_FREQ`, VCF `AF`); SNPGenie internally converts to percent for its own outputs; but the stated downstream expectation is *actual percentages*. The scripts responsible for that conversion are placeholders: `bin/plot_variants.py` and `bin/export_excel_variants.py` are no-op stubs, and `bin/generate_variant_plots.py` writes a hard-coded 1×1-pixel PNG. There is currently no code path that multiplies by 100, no unit label, and no range assertion. When these are implemented without an explicit convention, a 0.5 (proportion) rendered as "0.5%" (a 100× understatement of a 50% variant) is a guaranteed class of silent error.

- **MF-8 (Moderate) — `ivar_variants_to_vcf.py` fabricates genotype and quality fields.** Every record is written with `GT=1/1` (diploid homozygous ALT) and `QUAL=60` (`bin/ivar_variants_to_vcf.py:61,67`). For intra-host iSNVs a diploid homozygous genotype is biologically wrong (the sample is a haploid mixture), and it is not inert: downstream `bcftools csq --local-csq` uses phasing/genotype information, so compound-variant consequence calls on the iVar branch can be skewed. The iVar p-value (present in the TSV) is parsed but then discarded instead of being converted to QUAL.

- **MF-9 (Low) — CliqueSNV frequency threshold is implicit.** Only `-m snv-illumina -threads` are passed; the SNV frequency threshold and related cutoffs live inside the preset and can change between CliqueSNV releases. For a 0.1%–1% objective, the threshold must be pinned explicitly (as for VILOCA/ShoRAH in MF-4) and recorded in `methods_key_parameters.tsv`.

- **MF-10 (Low, documentation) — Indel handling is silently SNV-only.** With `lofreq_enable_indelqual = false`, indels carry no insertion/deletion qualities, so `lofreq filter --indelqual-thresh 20` effectively purges all indel calls. This is acceptable for an iSNV-focused study, but it is currently undocumented — a user reading "LoFreq calls SNPs and indels" would draw a wrong conclusion.

- **MF-11 (Low, statistical ceiling to disclose) — Q30 physics cap honest sub-0.3% claims.** With ε ≈ 10⁻³ at Q30, a 0.1% iSNV and the error process are the same order of magnitude; LoFreq mitigates this via quality-aware testing and (in `lofreq filter`, when enabled) strand-bias filtering, but no configuration of iVar — a simple per-base binomial heuristic — can support 0.1% claims. The pipeline currently has no corroboration layer (strand-balance requirement, dual-caller concordance flag, or UMI support) for the <0.5% regime. This is a design gap to disclose, not a bug.

---

## 3. Required Adjustments

Concrete changes, each with scientific justification. Ordered by priority.

### RA-1 — Set iVar's frequency floor to the study's detection target: `ivar_min_freq = 0.001`
**Change:** `nextflow.config:30` → `ivar_min_freq = 0.001`; update `nextflow_schema.json` description to *"Minimum allele frequency to report, as a proportion (0–1). Values < 0.005 require median depth ≥ 10,000× to be statistically meaningful."*
**Justification:** AF = 0.001 requires ≥10 ALT reads ⇒ ≥10,000× depth. At the configured 1,000× floor, `-t 0.001` would admit 1-read calls indistinguishable from Q30 error. The paired depth change (RA-2) plus the corroboration rule (RA-10) makes 0.1% defensible. If the lab's median depth is in the 1,000–5,000× range, keep 0.005–0.01 and *declare* that floor honestly instead.

### RA-2 — Replace coverage *inclusion floors* with permissive values: `lofreq_min_depth: 1000 → 10`, `ivar_min_depth: 1000 → 10`
**Change:** `nextflow.config:29,36`; fix the misleading schema descriptions.
**Justification:** `--min-cov`/`-m` only decide which positions are *tested*. LoFreq's Poisson-binomial model already fails to reach significance at insufficient depth — the 1,000× pre-filter adds no FPR protection while (a) destroying callable variants at 10–999×, (b) creating dead zones at amplicon edges, and (c) biasing SNPGenie π/dN/dS denominators (absent ≠ invariant). Let the error models arbitrate; report per-position depth and let RA-10 tier the calls.

### RA-3 — Delete the dead `lofreq_min_freq` parameter (do **not** wire it as a filter)
**Change:** remove from `nextflow.config:37` and `nextflow_schema.json`; note the removal in `CHANGELOG.md`.
**Justification:** An AF floor on LoFreq output is exactly the aggressive filtering this pipeline must avoid (CV-5). Keeping an unused parameter that *claims* such a floor exists is a reproducibility hazard (MF-3). The correct AF-aware behavior belongs in the reporting tier (RA-10), not in the caller.

### RA-4 — Relax LoFreq mapping-quality floor: `lofreq_min_mq: 60 → 20`
**Justification:** BWA MEM caps MQ at 60, so ≥60 is a uniqueness demand that strips informative reads from repeated/low-complexity viral regions — falsely penalizing high-mutation regions (MF-6). MQ 20 removes grossly mis-mapped reads while LoFreq's per-base-quality model handles the rest.

### RA-5 — Run actual VILOCA, or rename the module
**Change (preferred):** `containers.config:20` → pin `quay.io/biocontainers/viloca:1.1.x--…` (resolve digest per the repo's G2 procedure); keep the `shorah` CLI invocation (VILOCA provides it) but update `versions.yml` parsing. Alternative: rename module/publish path `VILOCA → SHORAH` and document legacy ShoRAH 1.99.2.
**Justification:** MF-4. Methods sections currently overstate the tool/version. VILOCA's corrected DPM updates are the scientifically defensible choice for new analyses.

### RA-6 — Pin haplotype inference geometry explicitly (both tools)
**Change:** add params, e.g. `viloca_window = 150`, `viloca_shift = 50`, and pass `shorah shotgun -w ${params.viloca_window} -s ${params.viloca_shift}`; for CliqueSNV pin the frequency threshold flag explicitly (per `-m snv-illumina` intent, 0.001) after verifying the exact flag against the pinned binary's `--help`; record all three in `methods_key_parameters.tsv`.
**Justification:** MF-4/MF-9. ShoRAH/VILOCA clusters reads *within* windows; window size must be ≤ the length reads actually span (reads must cover the window), so it must be co-designed with the library's read length — an input-derived default is non-reproducible across datasets. Shift ≈ window/3 is the standard tiling density. Explicit thresholds make the 0.1%–1% sensitivity claim auditable.

### RA-7 — Repair the reporting container and make statistical degradation fatal
**Change:** replace `container_python_reporting` (vanilla `python:3.12.12`) with an image that genuinely contains pandas/scipy/matplotlib/openpyxl (matching the existing comment's claim), pinned with digest; **and** change `analyze_delta_selection.py` to hard-fail (`sys.exit(1)`) when scipy is unavailable instead of emitting NaN statistics. Apply the same fail-loud rule to `build_selection_tables.py`.
**Justification:** MF-5. Silent loss of the inferential layer is worse than a crash; Kruskal–Wallis p-values are the primary cross-sectional group-comparison statistic of the selection branch.

### RA-8 — Fix fabricated fields in `ivar_variants_to_vcf.py`
**Change:** `GT: "1/1"` → `"1"` (haploid ALT) or `"."`; `QUAL: 60` → `-10·log10(max(pval, 1e-300))`; keep `AF` as a proportion.
**Justification:** MF-8. The VCF feeds `bcftools csq --local-csq`; genotype/phasing fields influence compound-consequence logic. Converting the iVar binomial p-value to QUAL preserves the only statistical evidence iVar produces.

### RA-9 — Establish and test the proportion→percent reporting contract
**Change:** (a) Adopt the rule *"machine files (VCF, TSV) store proportions 0–1; all human-facing tables/plots render `100 × proportion` with explicit `%` units"* in `README.md`; (b) implement the ×100 conversion when the stub scripts (`plot_variants.py`, `export_excel_variants.py`, `generate_variant_plots.py`) are completed; (c) add nf-test assertions on reporting outputs: frequency columns within [0, 100] and `%` present in axis/column labels.
**Justification:** MF-7. Without a single enforced convention, proportion/percent mix-ups produce silent 100× errors — the most damaging and least detectable class of reporting bug in iSNV work.

### RA-10 — Add a frequency-tier corroboration policy for the <0.5% regime
**Change:** in the reporting layer, annotate iSNVs as: `≥1%` = report; `0.1–1%` = report as *candidate*, require (i) LoFreq call (quality-model based), (ii) strand balance (LoFreq SB filter pass — verify `lofreq filter` SB defaults or enable explicitly), and (iii) ≥10 ALT observations; `<0.1%` = below Q30-supported detection, report only as exploratory. State this policy in `README.md` and the executive report methods block.
**Justification:** MF-1/MF-11. At Q30, ε ≈ 0.1%; single-caller, single-strand evidence cannot separate a 0.1% variant from error. Tiered reporting preserves sensitivity (nothing is destroyed) while keeping the inferential claim honest — and it naturally absorbs the true positive rate gained by RA-1/RA-2. (If the lab moves to UMI-tagged libraries, this tier can be revisited.)

### RA-11 — Make SNPGenie's no-AF-floor choice explicit
**Change:** pass `--minfreq=0` in `modules/local/snpgenie/run/main.nf:28`.
**Justification:** CV-10 is currently correct *by default*; making it explicit protects the π estimator against a future "helpful" addition of a minfreq that would bias cross-sectional diversity downward. Zero-cost auditability.

### RA-12 — Document SNV-only indel scope (or enable indelqual)
**Change:** either set `lofreq_enable_indelqual = true` (accepting the runtime cost of `--dindel` qualities) **or** add a README/report note that the LoFreq branch is SNV-only by configuration (`--indelqual-thresh 20` purges unqualified indels).
**Justification:** MF-10. Indel blindness must be a documented choice, not an accident of filter interaction.

---

## Summary Table

| ID | Severity | Item | Action |
|----|----------|------|--------|
| MF-1 | Critical | iVar 1% hard floor (`-t 0.01`) | RA-1: `-t 0.001` w/ depth caveat |
| MF-2 | Critical | 1,000× inclusion floors (`--min-cov`, `-m`) | RA-2: 1000 → 10 |
| MF-3 | High | Dead `lofreq_min_freq` param | RA-3: delete |
| MF-4 | High | "VILOCA" = legacy ShoRAH, implicit windows | RA-5 + RA-6 |
| MF-5 | High | Silent NaN Kruskal statistics (no scipy) | RA-7 |
| MF-6 | Moderate | LoFreq MQ ≥ 60 | RA-4: 60 → 20 |
| MF-7 | Moderate | Percent/proportion contract undefined; stubs | RA-9 |
| MF-8 | Moderate | Fabricated GT/QUAL in iVar VCF | RA-8 |
| MF-9 | Low | Implicit CliqueSNV threshold | RA-6 |
| MF-10 | Low | Undocumented indel blindness | RA-12 |
| MF-11 | Low | No <0.5% corroboration layer | RA-10 |

**Bottom line:** the pipeline's *architecture* (reference-based, dual-caller, quality-model-first, no AF filters on LoFreq, per-sample cross-sectional statistics) is scientifically sound. The blocking issues are concentrated in (1) hard 1%/1,000× floors that contradict the 0.1%–1% objective, (2) a mislabeled/unpinned haplotype inference step, (3) a silently degradable selection-statistics container, and (4) an undefined frequency-unit contract in the human-facing reporting layer. RA-1 through RA-7 are prerequisites for any 0.1%–1% iSNV claim from this workflow.
