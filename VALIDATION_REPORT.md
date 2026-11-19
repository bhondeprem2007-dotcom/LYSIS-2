# LYSIS Research Desk 4.0 — Validation Report

**Build:** 4.0-sprint · local educational prototype

## Automated results

Command: `py -m unittest discover -s tests -v`

**Result in build environment: 62 tests passed, 0 failed.** This includes the previous 54 test cases and 8 new cases for six-frame translation and IUPAC motif scanning.

Coverage includes homepage and health routes, search ranking shape and bounds, source note presence, paper catalogue and detail routes, comparison input validation, evidence brief generation, FASTA syntax and composition, sequence-analysis utilities, translation behavior, motif scanning, malformed/oversized requests, unknown routes, security headers and corpus structural invariants.

## Manual Windows acceptance still required

- Extract package outside ZIP preview.
- Install dependencies with `py -m pip install -r requirements.txt`.
- Start with `py app.py`; confirm `http://127.0.0.1:8765`.
- Test one literature query, one FASTA input, one translation input and one motif scan.
- Confirm downloads for catalogue and protein FASTA work in the browser.

## Limitations

Automated tests establish expected software behavior for the covered cases; they do not establish biological validity, source-record correctness, model performance, usability on every Windows installation or clinical utility. The application has not undergone independent security review, laboratory validation or clinical validation.
