# LYSIS 4.0 — Next Sprint Changelog

## Added
- Six-frame translation endpoint and UI using the standard genetic code.
- Ambiguous codons represented as X, with per-frame counts of amino acids, ambiguous codons and stop codons.
- Export of six translated frames as a FASTA protein (`.faa`) file.
- IUPAC motif scan endpoint and UI, supporting degenerate nucleotide codes, overlapping hits and both orientations.
- 1-based hit coordinates and bounded output (first 500 matches, with truncation flag).
- Paper catalogue export formats documented: CSV, JSON, RIS and BibTeX.

## Validation
- 62 automated unit/integration tests passed in the build environment.
- Added tests for translation frame count, standard codon translation, ambiguous codons, invalid alphabets and empty inputs.
- Added tests for forward and degenerate motif matching, overlapping hits and invalid patterns.
- Existing 54 v3 tests retained.
- Browser JavaScript syntax checked separately; local Windows UI acceptance remains required.

## Scientific boundaries
- Translation is a six-frame sequence preview, not gene calling or protein annotation.
- Motif hits do not establish restriction enzyme activity, promoter function, host range, pathogenicity or treatment value.
- No new biological dataset, phage-host prediction model or clinical decision support was added.
