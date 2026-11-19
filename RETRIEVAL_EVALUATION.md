# Retrieval smoke evaluation — 19-record starter corpus

**Evaluation date:** 19 November 2026  
**Method:** five manually written topic queries; top-five results; project-team relevance labels based on the corpus titles and short research notes.

| Query | Top-five IDs (ranked) | Precision@5 |
|---|---|---:|
| Pseudomonas phage-antibiotic cocktails | P008, P003, P015, P011, P019 | 1.00 |
| Pseudomonas biofilm / phage resistance | P009, P015, P008, P013, P019 | 0.80 |
| Personalized therapy / clinical outcomes | P002, P001, P015, P018, P008 | 0.40 |
| MedCPT / biomedical retrieval | P007, P006, P008, P014, P005 | 0.40 |
| In-vitro synergy / clinical translation | P018, P015, P008, P003, P014 | 0.60 |

**Aggregate:** Hit@5 = 1.00; mean Precision@5 = 0.64; mean reciprocal rank at 5 = 1.00.

## Interpretation

The first relevant record appeared at rank 1 for all five queries, but only 64% of the top-five slots were labelled relevant on average. The personalized-therapy and MedCPT topics had lower precision, indicating that the small starter corpus and current lexical/LSA weighting need refinement. The MedCPT query is also a narrow topic represented by only two directly relevant records.

This is a smoke test, not an independently labelled benchmark. The relevance labels were created by the project team from the limited corpus itself; they are not blinded, independently adjudicated or representative of real user queries. The metrics should not be presented as a general performance claim. A proper evaluation would require a larger query set, independent domain-expert relevance labels, held-out queries and comparison against baseline retrieval methods.

Reproduce with:

```powershell
py evaluate_retrieval.py
```
