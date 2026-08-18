<div align="center">

# GRASP — Command-Line Interface

`version 14-Aug-2026`

</div>

---

## Table of Contents

- [What can GRASP do?](#what-can-grasp-do)
- [Quick Setup — Download & Run](#quick-setup--download--run)
- [Building from Source — Maven](#building-from-source--maven)
- [Building from Source — IntelliJ](#building-from-source--intellij)
- [Running GRASP](#running-grasp)
- [MIP Indel Inference](#mip-indel-inference)
- [Important Command-Line Options](#important-command-line-options)
- [Full CLI Reference](#full-cli-reference)
- [Developing in IntelliJ](#developing-in-intellij-no-jar-build)

---

## What can GRASP do?

GRASP accepts an alignment (**FASTA** or **Clustal**) and a phylogenetic tree (**Newick**) with matching labels, and infers ancestor sequences by joint or marginal maximum-likelihood reconstruction. Along the way it also infers insertion/deletion (indel) events, represented internally as **partial-order graphs (POGs)**, and identifies the single most-supported path through each ancestor's graph.

It can save:
- All ancestor sequences (joint reconstruction), or one sequence with character-state distributions (marginal reconstruction)
- Partial-order graphs as JSON or DOT (viewable in GraphViz)
- The tree, re-saved with assigned ancestor labels

> **Note:** GRASP was built primarily for protein sequences. The CLI also supports DNA models, but DNA-specific functionality (codon-centric analysis, custom background stats, etc.) hasn't been extensively tested or developed yet.

### Main indel inference methods in GRASP

| Approach | Behavior |
|---|---|
| **Bi-directional edge (BE)** — *default* | Efficient; can represent ambiguous indel histories in a single POG |
| **Mixed Integer Programming (MIP)** | More conservative and computationally heavier, but guarantees a globally optimal solution across all sites simultaneously |

On top of BE and MIP, indels can also be coded as **Position Specific (PS)** or **Simple Indel Coding (SIC)**, and inferred by **Parsimony (P)** or **Maximum Likelihood (ML)** — six combined methods in total, plus the two MIP solvers, **SCIP** and **Gurobi** (same optimal answer; Gurobi is commercial and substantially faster).

> Under 1,000 sequences usually finishes in under 10 minutes. Up to ~10,000 sequences is doable on a server (a "gappy" alignment slows things down). A clean alignment can run on a modern laptop with 16GB RAM in under a day.

---

## Quick Setup — Download & Run

1. **Install Java 11+.** Any OS works — macOS, Windows, Linux.
2. **Download the jar** from the [bnkit releases page](https://github.com/bodenlab/bnkit/releases).
3. **Sanity-check it:**

   ```console
   java -jar bnkit.jar -h
   ```

   You should see the full help text print to the console.

---

## Building from Source — Maven

1. Clone [bnkit](https://github.com/bodenlab/bnkit) in full. (JUnit 5 is only needed if you're a developer running the test suite.)
2. Confirm Maven is installed:

   ```console
   mvn -h
   ```

   If that fails, install it from [maven.apache.org](https://maven.apache.org/download.cgi), or via a package manager (`apt` on Ubuntu/Debian, `brew` on macOS).

3. From the repo root, build the jar:

   ```console
   mvn package -f poms/pom.xml -DskipTests
   ```

   `-DskipTests` skips the test suite which is fine unless you're developing. The jar is placed in `target/GRASP-<version>.jar`.

---

## Building from Source — IntelliJ

1. **File → New → Project From Version Control**
2. URL: `http://github.com/bodenlab/bnkit`
3. Choose a save directory
4. **File → Project Structure → Artifacts** → click `+` → *Jar → From modules with dependencies*
5. Select module `asr`, main class `asr.GRASP` → **OK** → **Apply**
6. **Build → Build Artifacts → `bnkit:jar` → Build**
7. Output lands in `out/artifacts/bnkit_jar`

---

## Running GRASP

1. **Get the jar** — build it or download it (above). You can also run it directly:

   ```console
   java -jar ~/Downloads/grasp.jar
   ```

2. **Wrap it in a launcher script** called `grasp`:

   ```bash
   #!/bin/sh
   java -jar -Xmx16g </path/to/grasp.jar> $@
   ```

   > 💡 **Memory tip:** `-Xmx` sets the max memory heap. For large and/or gappy alignments budget more — **60GB** (`-Xmx60000m`) is a reasonable recommendation.

3. **Make it executable:**

   ```console
   chmod 755 grasp
   ```

4. **Put it on your PATH:**

   ```console
   mv grasp /usr/local/bin
   ```

5. **Confirm it works:**

   ```console
   grasp -h
   ```

### Example invocation

```console
grasp --aln 500_2112_dhad_18032019.aln \
      --nwk r_500_2112_dhad_18032019.nwk \
      --output-folder recon_0500 \
      --verbose --threads 5
```

---

## MIP Indel Inference

By default GRASP uses **bi-directional edge parsimony (BEP)** — fast and accurate on large datasets, but the paths through ancestral sequences are not guaranteed to be optimal. **MIP** considers every site and branchpoint under one objective function to guarantee a globally optimal solution.

| Solver | Type | Notes |
|---|---|---|
| **SCIP** | Open-source, bundled | Use `--indel-method SCIP` — no extra install needed |
| **Gurobi** | Commercial | Substantially faster, handles very large datasets. Free academic licenses at [gurobi.com](https://www.gurobi.com/). Install per [Gurobi's setup guide](https://support.gurobi.com/hc/en-us/articles/14799677517585-Getting-Started-with-Gurobi-Optimizer), then use `--indel-method Gurobi` — GRASP auto-detects it |

---

## Important Command-Line Options

**Site-specific evolutionary rates** — `--rates-file <filename>`
Improves accuracy on diverse families using per-site rates (e.g. from IQ-TREE2). Tab-separated:

```
# Columns are tab-separated with following meaning:
#   Site:   Alignment site ID
#   Rate:   Site rate estimated by maximum likelihood
Site    Rate
1       2.51550
2       12.89129
3       34.31350
4       2.44313
```

**Indel inference method** — `--indel-method <methodname>`
Default `BEP`. See the [method table](#indel-methods) below.

**Substitution model** — `--substitution-model <modelname>`
Default `JTT`. `Dayhoff`, `LG`, `WAG` also available for protein; `JC` and `Yang` for DNA.

**Empirical base frequencies** — `--empirical-freqs <filename>`
TSV of stationary character frequencies:

```
Character Proportion
A 0.25
C 0.25
T 0.25
G 0.25
```

**Threads** — `--threads <number>`
Default `4`. More threads only help up to a point.

---

## Full CLI Reference

### Required

| Flag | Description |
|---|---|
| `-a, --aln <file>` | Multiple-sequence alignment, FASTA or CLUSTAL format |
| `-n, --nwk <file>` | Phylogenetic tree, Newick format, with labels matching the alignment |

### Common optional flags

| Flag | Description |
|---|---|
| `-o, --output-folder` | Where output files are written (default: cwd, or input folder if set) |
| `-i, --input-folder` | Skip indel inference; load a prior reconstruction |
| `-pre, --prefix` | Stub prepended to result filenames (default: alignment file's prefix) |
| `-sa, --save-as` | Which output formats to generate (see [table](#output-formats)) |
| `--save-all` | Generate every output format |
| `-rf, --rates-file` | Position-specific substitution rates (recommended for diverse/distant trees) |
| `-ef, --empirical-freqs` | Stationary character frequencies for the substitution model |
| `-j` / `-m <branchpoint-id>` | Joint reconstruction (default) or marginal at a specific branchpoint |
| `--onlyindel` | Skip character-state inference; infer indels only |
| `--include-extants` | Include extants in output files, where the format allows |
| `--nogap` | Exclude the gap character from output, where allowed |
| `--nonibble` | Keep POG indices that can't form a start-to-end path (don't trim them) |
| `--orphans` | Keep orphaned indel trees (don't remove them) |
| `--exclude-noedge` | Remove "non-existing edge" as a parsimony option in BEP |
| `--solver-time-limit` | Minutes before the MIP solver gives up and falls back to BEP |
| `--supported-path` | `DIJKSTRA` (default) or `ASTAR` |
| `--save-tree` | Skip inference; re-save the tree with GRASP's depth-first ancestor labels (starting `N0`) |
| `--save-poag { branchpoint-id }` | Skip inference; save the input alignment as a POAG under a given ancestor (default `N0`) |
| `--seed` | Random seed |
| `--verbose` / `--time` | Print progress / print total runtime |
| `-h, --help` | Print help |

### Indel methods

| Method | What it is |
|---|---|
| `BEP` *(default)* | Bi-directional edge, maximum parsimony |
| `BEML` | Bi-directional edge, maximum likelihood (uniform / JC-like model) |
| `SICP` | Simple indel-coding, maximum parsimony (Simmons & Ochoterena) |
| `SICML` | Simple indel-coding, maximum likelihood (uniform model) |
| `PSP` | Position-specific, maximum parsimony |
| `PSML` | Position-specific, maximum likelihood (uniform model) |
| `SCIP` | Globally optimal, distance-sensitive parsimony via open-source [SCIP](https://www.scipopt.org/). No multi-threading support |
| `Gurobi` | Same guarantee as SCIP, via commercial Gurobi solver — requires local install |

> Append `*` to a method name for a less conservative setting (where available), or a simple gap-opening penalty of 2 for MIP solvers.

### Substitution models

| Model | Type | Notes |
|---|---|---|
| `JTT` | Protein | Jones–Taylor–Thornton — default |
| `Dayhoff` | Protein | Dayhoff–Schwartz–Orcutt |
| `LG` | Protein | Le–Gasquel |
| `WAG` | Protein | Whelan–Goldman |
| `JC` | DNA | Jukes–Cantor |
| `Yang` | DNA | Yang's general reversible process model |

### Output formats

| Format | Contains |
|---|---|
| `FASTA` | Most-preferred-path sequences per ancestor, gapped or not |
| `CLUSTAL` | Most-preferred-path sequences per ancestor, gapped |
| `TREE` | Phylogenetic tree with ancestor nodes labelled |
| `DISTRIB` | Per-position character distributions (marginal reconstruction only), indexed by POG |
| `ASR` | Full reconstruction as JSON — POGs of ancestors + extants, and tree (`ASR.json`) |
| `DOT` | Ancestor partial-order graphs in DOT format (view with GraphViz) |
| `TREES` | Position-specific trees with ancestor states labelled |

### Full usage synopsis

```text
Usage: asr.GRASP 
        [-a | --aln <filename>]
        [-n | --nwk <filename>]
        {-o | --output-folder <foldername>}
        {-i | --input-folder <foldername>}
        {-pre | --prefix <stub>}
        {-rf | --rates-file <filename>}
        {-ef | --empirical-freqs <filename>}
        {-s | --substitution-model <JTT(default)|Dayhoff|LG|WAG|JC|Yang>}
        {-t | --threads <number>}
        {-j | --joint (default)}
        {-m | --marginal <branchpoint-id>}
        {--indel-method <methodname>}   (BEP default, BEML, SICP, SICML, PSP, PSML, SCIP, Gurobi)
        {--supported-path <methodname>} (DIJKSTRA default, ASTAR)
        {--nogap} {--seed <seed>} {--nonibble} {--exclude-noedge}
        {--save-as <list-of-formats>}   (FASTA CLUSTAL TREE DISTRIB ASR DOT TREES)
        {--save-all}
        {--save-tree}
        {--save-poag { <branchpoint-id> }}
        {--time} {--verbose} {--help}
```

---

## Developing in IntelliJ

It can be useful to run GRASP from within IntelliJ IDEA without having to build a
jar file. This is especially useful for debugging and development. The following steps will help you set up GRASP in IntelliJ:

1. **File → New → Project From Version Control** — URL: `http://github.com/bodenlab/bnkit`
2. **Run → Edit Configurations → `+` → Application** — name it (e.g. *GRASP*), set main class to `asr.GRASP`
3. **Run it** — green arrow, `Ctrl-R`, or Run → Run GRASP
4. **Fix the SDK if needed** — requires Java 11+ in project settings
5. **Add CLI arguments** — Run → Edit Configurations → "Program arguments" textbox

---
