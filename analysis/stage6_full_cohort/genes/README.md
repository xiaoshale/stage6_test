# Full-cohort gene identity preparation

Prepared 2026-11-16. Scope: native-feature identity and coverage evidence for human review.
Original count matrices are not opened or changed by these modules. Source counting references remain distinct from the common identity target.

## Frozen common references

`references_20261116/reference_manifest.json` records four official assets, exact URLs, HTTP metadata, acquisition time, byte counts and SHA-256:

| Asset | Bytes | Verified content |
|---|---:|---|
| Ensembl BioMart registry | 1,897 | `ensembl_mart_116` from the June 2026 official archive |
| Ensembl 116 human gene axis | 4,772,369 | 86,411 unique gene IDs; version, symbol, chromosome and biotype |
| HGNC approved complete set | 17,007,798 | 45,395 approved HGNC records |
| HGNC withdrawn table | 259,101 | 5,293 withdrawn symbols and replacement reports where supplied |

Total downloaded: 22,041,165 bytes, within the 30 MiB acquisition limit. These rows describe reference tables, not measured genes in any project dataset. Complete historical Ensembl split/merge mappings have not been acquired; source IDs absent from this target remain held.

The HGNC snapshots are dated by their saved acquisition/HTTP metadata and content hashes. Subsequent runs do not silently redownload a mutable latest table. Source annotation release and exon/intron counting definitions remain source-specific unknowns unless independently recovered.

## Map a native feature axis

```bash
python3 map_features.py \
  --features /path/to/native_features.tsv \
  --resource GSE135893 \
  --references references_20261116/reference_manifest.json \
  --output-dir /path/to/gene_mapping_output
```

Required TSV columns: `source_row_1based`, `source_gene_id`, `source_symbol`, `genome`, `feature_type`.
The source gene ID may be blank for symbol-only exports. Rows must preserve the full native axis in original contiguous order. Extra columns are not treated as expression.

Outputs are `gene_mapping.tsv` and `mapping_summary.json`. Each mapping retains source row, symbol, complete original ID, source feature-file hash, canonical candidates, target identity, method, confidence/hold reasons and frozen-reference hashes.

- **Strict identity eligibility:** unique source ENSG present in the fixed release, with acceptable genome evidence and no source-ID/symbol/row conflict. This does not establish equivalent quantification across references.
- **Symbol-only provisional:** unique HGNC approved/previous/alias correspondence to a target gene; retained for review and exploratory identity matching while the original reference crosswalk is unresolved. It is not silently promoted to the strict set.
- **Held:** ambiguous alias/history, absent target-release ID, conflicting ID/symbol, unsupported feature identity or multiple source rows sharing a canonical ID. No first-match selection, count aggregation or count duplication occurs.
- Valid Ensembl version suffixes are parsed; punctuation in symbols such as `RP11-34P13.3` is preserved.
- Source `MT-` flags and common-reference mitochondrial chromosome annotations are separate metadata. No QC denominator is changed.

## Describe cross-axis coverage

```bash
python3 coverage.py \
  --mapping GSE135893=/path/to/one/gene_mapping.tsv \
  --mapping GSE159354=/path/to/two/gene_mapping.tsv \
  --output-dir /path/to/coverage_output
```

Outputs distinguish source-ID unique, symbol-only provisional, mapping held and absence from the mapped export. They provide both shared strict features and shared unique candidates including provisional symbols. A zero strict intersection caused by symbol-only input is reported transparently; it does not mean the studies have no shared biological genes. Export absence remains coverage unknown, not zero expression or proof of an off-panel gene.

These are identity/coverage tables. Derived count projections, their included count-mass conservation and marker/program adequacy require separate verification in the full preparation workflow. No matrix harmonization or final cohort eligibility is implied by a successful metadata mapping.

## Full native intake execution

`map_all_axes.py --run RUN --references MANIFEST --output-dir NEW_OUTPUT` reads the three completed RNA input manifests, completed spatial summary and locked intake. It requires exact coverage of every source-qualified unit and validates the declared feature hash, length and donor crosswalk. Axes are deduplicated by the complete ordered native metadata content, preserving raw-file checksums and each unit's mapping path.

`identity_marker_coverage.py --batch NEW_OUTPUT/batch_mapping_summary.json --seed-code RUN/code/rna_review.py --output-dir NEW_MARKER_OUTPUT` audits the fixed identity-marker dictionary in the actual RNA review code. This describes source-label and mapping availability, not observed expression or annotation accuracy.

The actual execution covered all 134 native units in `preparation_001`; see [RESULTS.md](RESULTS.md) and the verified small local copies under `results_20261116/`. Large row-level mappings remain on the server. No source count matrix was read or changed by this gene-metadata execution.

## Narrow verification

Run `python3 check_mapping_fixture.py`. Twelve synthetic checks cover version parsing, dotted symbols, alias ambiguity, retired IDs, ID/symbol conflicts, provisional symbols, previous names, non-gene features, alternate genome review, duplicate-source/canonical rows and source-order failure. `mapping_fixture_result.json` records the result. These fixtures do not establish real-dataset duplicate counts.

The actual official references were separately loaded with SHA/size checks, complete axis parsing and duplicate-ID validation. Root execution must run the mapping CLI against every distinct native axis and preserve the axis-to-library/section crosswalk before asserting complete feature coverage.
