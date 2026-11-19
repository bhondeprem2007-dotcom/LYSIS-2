#!/usr/bin/env python3
"""Small, manually labelled retrieval smoke evaluation for the starter corpus.

The relevance labels are a project-team topic judgement based on each record's
bibliographic title and summary. This is not an independently curated IR gold
standard and must not be reported as a generalizable performance benchmark.
"""
from __future__ import annotations
import json
from pathlib import Path
import build_embeddings as retrieval

ROOT = Path(__file__).resolve().parent
retrieval.CORPUS = ROOT / "rag_corpus.jsonl"
retrieval.VECTORS_PATH = ROOT / "lsa_vectors.npy"
retrieval.VOCAB_PATH = ROOT / "lsa_vocab.json"
retrieval.META_PATH = ROOT / "chunk_meta.json"
retrieval.MANIFEST_PATH = ROOT / "lsa_index_manifest.json"

CASES = [
    {
        "query": "Pseudomonas aeruginosa phage antibiotic cocktails",
        "relevant": {"P003", "P008", "P009", "P010", "P011", "P012", "P014", "P015", "P016", "P017", "P018", "P019"},
    },
    {
        "query": "Pseudomonas aeruginosa biofilm phage resistance trade-offs",
        "relevant": {"P008", "P009", "P010", "P011", "P012", "P013", "P016", "P017", "P019"},
    },
    {
        "query": "personalized phage therapy clinical outcomes consecutive cases",
        "relevant": {"P001", "P002", "P004", "P005"},
    },
    {
        "query": "MedCPT biomedical information retrieval PubMed transformer",
        "relevant": {"P006", "P007"},
    },
    {
        "query": "in vitro phage antibiotic synergy clinical translation Pseudomonas",
        "relevant": {"P008", "P015", "P016", "P017", "P018"},
    },
]

def main() -> int:
    if not retrieval.index_is_current():
        retrieval.build_lsa(retrieval.load_corpus(ROOT / "rag_corpus.jsonl"))
    rows = []
    for case in CASES:
        results = retrieval.search(case["query"], 5)
        ranked = [r["paper_id"] for r in results]
        relevant = case["relevant"]
        hits = [i + 1 for i, paper_id in enumerate(ranked) if paper_id in relevant]
        rows.append({
            "query": case["query"],
            "top_5": ranked,
            "relevant_ids": sorted(relevant),
            "hits_at_5": len(hits),
            "precision_at_5": round(len(hits) / 5, 3),
            "hit_at_5": bool(hits),
            "reciprocal_rank": round(1 / hits[0], 3) if hits else 0.0,
        })
    summary = {
        "evaluation_type": "small manually labelled smoke set; not a formal benchmark",
        "queries": len(rows),
        "mean_precision_at_5": round(sum(r["precision_at_5"] for r in rows) / len(rows), 3),
        "hit_rate_at_5": round(sum(r["hit_at_5"] for r in rows) / len(rows), 3),
        "mean_reciprocal_rank_at_5": round(sum(r["reciprocal_rank"] for r in rows) / len(rows), 3),
        "cases": rows,
    }
    out = ROOT / "retrieval_evaluation.json"
    out.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
