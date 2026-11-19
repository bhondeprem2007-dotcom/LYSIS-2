# LYSIS 4.0 — Windows manual acceptance checklist

- [ ] Extract ZIP to a normal folder.
- [ ] Run `py -m pip install -r requirements.txt`.
- [ ] Run `py app.py`; homepage opens at `http://127.0.0.1:8765`.
- [ ] Literature search returns source records and DOI links.
- [ ] Paper explorer filter and CSV/JSON/RIS/BibTeX exports work.
- [ ] Compare 2–4 known paper IDs.
- [ ] Generate and download an evidence brief.
- [ ] Save a reading-list record and note; refresh and confirm persistence.
- [ ] FASTA checker validates a valid sequence and flags invalid symbols.
- [ ] Sequence lab reverse complement and k-mer output display.
- [ ] Six-frame translation: input `ATGGCC`; forward frame 1 begins `MA`; export `.faa`.
- [ ] IUPAC motif scan: sequence `AACCGGTT`, motif `CCG`; confirm forward hit at 1-based position 3.
- [ ] Motif scan: motif `AR` matches A/G at second position in suitable test input.
- [ ] Corpus audit displays structural status.
- [ ] System Scope accurately lists what is and is not implemented.
- [ ] Confirm no patient-identifiable data is entered.
