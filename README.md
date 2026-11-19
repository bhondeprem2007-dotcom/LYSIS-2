# LYSIS Research Desk — local workbench 4.0

## Run on Windows
1. Extract the ZIP to a normal folder (do not run from inside the ZIP preview).
2. Open the extracted `LYSIS_Prototype` folder.
3. Click the File Explorer address bar, type `cmd`, and press Enter.
4. Run `py -m pip install -r requirements.txt` and then `py app.py`.
5. Open **http://127.0.0.1:8765**. Keep the terminal open; press `Ctrl+C` to stop.
6. Alternatively, double-click `Run_LYSIS.bat` from the extracted folder.

## Modules in 4.0
- **Literature search:** local TF-IDF + NumPy truncated-SVD/LSA retrieval with transparent ranking components.
- **Paper explorer:** inspect the 19-record starter corpus; export CSV, JSON, RIS and BibTeX.
- **Compare papers:** side-by-side metadata and attributed corpus notes for 2–4 records.
- **Evidence brief:** create a downloadable source-linked Markdown compilation from retrieved records; this is not an LLM synthesis.
- **My reading list:** browser-local saved papers and personal notes.
- **FASTA checker:** IUPAC validation, lengths, GC%, base counts, ambiguity counts and N50.
- **Sequence lab:** reverse complement, k-mer counts and exploratory six-frame ORF scan.
- **Six-frame translation (new):** standard-code translation of all three forward and three reverse-complement reading frames; ambiguous codons render as X; export as `.faa`.
- **IUPAC motif scan (new):** overlapping motif matches on both orientations, including degenerate symbols, with 1-based coordinates.
- **Corpus audit and system scope:** structural metadata checks and explicit capability boundaries.

## Tests
Run `py -m unittest discover -s tests -v`. The 62-test suite covers API routes, validation, retrieval, FASTA, sequence utilities, six-frame translation, IUPAC motif scanning, corpus structure and evidence workflows. The build-environment run recorded 62 passing tests. Re-run on the target Windows PC for local acceptance.

## Scope and safety
Educational research-support software only. It does not implement MedCPT inference, LLM-generated answers, full-text PDF ingestion, validated gene annotation, sequence homology search, phage-host prediction, receptor matching, bacterial susceptibility testing, phage matching, antibiotic interaction prediction, candidate ranking or treatment recommendation. Translation and motif matching are descriptive computational utilities, not functional or clinical interpretations. The 19-record starter corpus is not the 103-paper catalogue described in project materials. Do not enter patient-identifiable information. Inspect original publications before citing or interpreting corpus notes.

## Privacy and provenance
The server binds only to `127.0.0.1`; reading-list data stays in browser local storage. The app makes no cloud upload or external API calls. DOI links open original publication sites only when clicked. The retrieval code is a reconstructed/revised implementation, not claimed to be the original historical file. `SHA256_MANIFEST.json` records hashes for this package build.
