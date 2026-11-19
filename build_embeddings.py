#!/usr/bin/env python3
"""Build and query a small LSA-based literature retrieval index for LYSIS.

This is a reconstructed and revised student-project component. It is not the
complete LYSIS application, a clinical decision-support system, or a validated
phage-selection tool.

Default retrieval pipeline:
  1. Read one JSON research record per line from rag_corpus.jsonl.
  2. Tokenize title, attributed summary, and subject areas.
  3. Build TF-IDF document vectors.
  4. Apply truncated SVD to obtain LSA document representations.
  5. Rank documents using a transparent weighted score:
        0.65 * LSA cosine similarity + 0.35 * query-term coverage.

Only NumPy is required. The MedCPT path mentioned in the earlier copyright
listing is not implemented by this script; it is treated as future work.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import re
import sys
from collections import Counter
from pathlib import Path
from typing import Any

import numpy as np

VERSION = "4.0-sprint"
HERE = Path(__file__).resolve().parent
# In the distributed project pack, data is stored beside this script in ../data.
DEFAULT_CORPUS = HERE.parent / "data" / "rag_corpus.jsonl"
# Also accept a corpus placed in the same directory as the script, as in the
# original user workflow.
LOCAL_CORPUS = HERE / "rag_corpus.jsonl"
CORPUS = LOCAL_CORPUS if LOCAL_CORPUS.exists() else DEFAULT_CORPUS
VECTORS_PATH = HERE / "lsa_vectors.npy"
VOCAB_PATH = HERE / "lsa_vocab.json"
META_PATH = HERE / "chunk_meta.json"
MANIFEST_PATH = HERE / "lsa_index_manifest.json"

STOP = set("""a about above after again against all am an and any are as at be because been before being
below between both but by can cannot could did do does doing down during each few for from further had has
have having he her here hers herself him himself his how i if in into is it its itself me more most my myself
no nor not of off on once only or other ought our ours ourselves out over own same she should so some such
than that the their theirs them themselves then there these they this those through to too under until up
very was we were what when where which while who whom why with would you your yours yourself yourselves
we our using used use also may can based show shows shown however thus therefore""".split())
TOKEN = re.compile(r"[a-z0-9][a-z0-9\-_]{1,}")
REQUIRED_FIELDS = ("chunk_id", "paper_id", "title", "text", "year", "venue", "doi")


def tokenize(text: str) -> list[str]:
    """Return normalized content tokens, excluding common stop words."""
    if not isinstance(text, str):
        return []
    return [t for t in TOKEN.findall(text.lower()) if t not in STOP and len(t) > 2]


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_corpus(corpus_path: Path | str | None = None) -> list[dict[str, Any]]:
    """Load JSONL records and fail clearly on malformed or incomplete rows."""
    path = Path(corpus_path) if corpus_path is not None else Path(CORPUS)
    if not path.is_file():
        raise FileNotFoundError(
            f"Corpus not found: {path}\n"
            "Place rag_corpus.jsonl beside this script or use the distributed data folder."
        )

    chunks: list[dict[str, Any]] = []
    seen_ids: set[str] = set()
    with path.open("r", encoding="utf-8") as handle:
        for line_no, raw in enumerate(handle, start=1):
            if not raw.strip():
                continue
            try:
                item = json.loads(raw)
            except json.JSONDecodeError as exc:
                raise ValueError(f"Invalid JSON on line {line_no}: {exc.msg}") from exc
            if not isinstance(item, dict):
                raise ValueError(f"Line {line_no} must contain a JSON object.")
            missing = [field for field in REQUIRED_FIELDS if field not in item]
            if missing:
                raise ValueError(f"Line {line_no} is missing required fields: {', '.join(missing)}")
            chunk_id = str(item["chunk_id"]).strip()
            if not chunk_id:
                raise ValueError(f"Line {line_no} has an empty chunk_id.")
            if chunk_id in seen_ids:
                raise ValueError(f"Duplicate chunk_id '{chunk_id}' on line {line_no}.")
            if not str(item["text"]).strip():
                raise ValueError(f"Line {line_no} has an empty text field.")
            if not str(item["title"]).strip() or not str(item["doi"]).strip():
                raise ValueError(f"Line {line_no} needs a non-empty title and DOI.")
            areas = item.get("areas", [])
            if isinstance(areas, str):
                areas = [areas]
            if not isinstance(areas, list):
                raise ValueError(f"Line {line_no}: areas must be a list of strings.")
            item["areas"] = [str(x) for x in areas]
            item["chunk_id"] = chunk_id
            item["paper_id"] = str(item["paper_id"])
            item["year"] = int(item["year"])
            item["source_url"] = str(item.get("source_url") or f"https://doi.org/{item['doi']}")
            item["authors"] = str(item.get("authors", ""))
            seen_ids.add(chunk_id)
            chunks.append(item)

    if len(chunks) < 2:
        raise ValueError(f"At least two corpus records are required for LSA; found {len(chunks)}.")
    return chunks


def searchable_text(chunk: dict[str, Any]) -> str:
    areas = chunk.get("areas", [])
    return " ".join([str(chunk.get("title", "")), str(chunk.get("text", "")), " ".join(areas)])


def _normalise_rows(matrix: np.ndarray) -> np.ndarray:
    norms = np.linalg.norm(matrix, axis=1, keepdims=True)
    return matrix / np.maximum(norms, 1e-12)


def build_lsa(chunks: list[dict[str, Any]], dim: int = 128) -> np.ndarray:
    """Build TF-IDF + truncated-SVD LSA vectors and metadata files."""
    if dim < 1:
        raise ValueError("--dim must be a positive integer.")
    if len(chunks) < 2:
        raise ValueError("At least two corpus records are required to build an LSA index.")

    docs = [tokenize(searchable_text(chunk)) for chunk in chunks]
    doc_freq = Counter()
    for doc in docs:
        doc_freq.update(set(doc))

    n_docs = len(docs)
    # Retain terms seen in at least one document. This is intentional for a
    # small starter corpus; dropping all hapax terms can erase niche biomedical
    # terms such as an organism or biofilm-specific phrase.
    vocab = sorted(term for term, count in doc_freq.items() if count < n_docs)
    if len(vocab) < 2:
        raise ValueError(
            "The corpus does not contain enough distinct indexable terms. "
            "Add more descriptive research text or subject-area tags."
        )

    vocab_index = {term: index for index, term in enumerate(vocab)}
    idf = np.array(
        [math.log((1 + n_docs) / (1 + doc_freq[term])) + 1.0 for term in vocab],
        dtype=np.float64,
    )

    tfidf = np.zeros((n_docs, len(vocab)), dtype=np.float64)
    for row, doc in enumerate(docs):
        tf = Counter(doc)
        if not tf:
            continue
        max_tf = max(tf.values())
        for term, count in tf.items():
            col = vocab_index.get(term)
            if col is not None:
                tfidf[row, col] = (0.5 + 0.5 * count / max_tf) * idf[col]

    tfidf = _normalise_rows(tfidf)
    if not np.any(tfidf):
        raise ValueError("No searchable terms were found in the corpus text.")

    # Direct SVD is numerically stable for this small project corpus.
    U, singular_values, Vt = np.linalg.svd(tfidf, full_matrices=False)
    max_rank = min(n_docs - 1, len(vocab) - 1, len(singular_values))
    rank = min(dim, max_rank)
    if rank < 1:
        raise ValueError("The corpus is too small to produce a non-zero LSA dimension.")

    components = Vt[:rank, :].T
    vectors = tfidf @ components
    vectors = _normalise_rows(vectors)
    if not np.isfinite(vectors).all() or not np.isfinite(components).all():
        raise FloatingPointError("SVD generated non-finite values; inspect the corpus and retry.")

    np.save(VECTORS_PATH, vectors.astype(np.float32))
    vocab_payload = {
        "version": VERSION,
        "vocab": vocab,
        "idf": idf.tolist(),
        "components": components.tolist(),
        "dimension": int(rank),
    }
    VOCAB_PATH.write_text(json.dumps(vocab_payload, ensure_ascii=False), encoding="utf-8")

    metadata: dict[str, dict[str, Any]] = {}
    for chunk in chunks:
        metadata[chunk["chunk_id"]] = {
            key: chunk.get(key, "")
            for key in ("paper_id", "title", "year", "venue", "doi", "authors", "areas", "source_url", "pmid", "record_type", "text")
        }
    META_PATH.write_text(json.dumps(metadata, ensure_ascii=False, indent=2), encoding="utf-8")

    manifest = {
        "version": VERSION,
        "algorithm": "TF-IDF + truncated SVD (LSA)",
        "dense_weight": 0.65,
        "lexical_weight": 0.35,
        "lexical_method": "query-term coverage (not BM25)",
        "requested_dimension": int(dim),
        "actual_dimension": int(rank),
        "record_count": len(chunks),
        "corpus_sha256": sha256_file(Path(CORPUS)),
        "corpus_filename": Path(CORPUS).name,
    }
    MANIFEST_PATH.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(
        f"built LSA index: {vectors.shape[0]} records x {vectors.shape[1]} dimensions, "
        f"vocabulary={len(vocab)}"
    )
    print(f"corpus: {CORPUS}")
    return vectors


def index_is_current() -> bool:
    paths = (VECTORS_PATH, VOCAB_PATH, META_PATH, MANIFEST_PATH)
    if not all(path.is_file() for path in paths):
        return False
    try:
        manifest = json.loads(MANIFEST_PATH.read_text(encoding="utf-8"))
        if manifest.get("version") != VERSION:
            return False
        if manifest.get("record_count") != len(load_corpus(Path(CORPUS))):
            return False
        if manifest.get("corpus_sha256") != sha256_file(Path(CORPUS)):
            return False
        vectors = np.load(VECTORS_PATH, mmap_mode="r")
        return vectors.ndim == 2 and vectors.shape[0] == manifest.get("record_count") and vectors.shape[1] == manifest.get("actual_dimension")
    except (OSError, ValueError, KeyError, json.JSONDecodeError):
        return False


def embed_query_lsa(query: str) -> tuple[np.ndarray, set[str]]:
    """Project a query into the indexed LSA space and return its known terms."""
    payload = json.loads(VOCAB_PATH.read_text(encoding="utf-8"))
    vocab = payload["vocab"]
    idf = np.asarray(payload["idf"], dtype=np.float64)
    components = np.asarray(payload["components"], dtype=np.float64)
    vocab_index = {term: index for index, term in enumerate(vocab)}
    query_terms = set(tokenize(query))
    known_terms = query_terms.intersection(vocab_index)
    if not known_terms:
        return np.zeros(components.shape[1], dtype=np.float64), set()

    counts = Counter(tokenize(query))
    vector = np.zeros(len(vocab), dtype=np.float64)
    max_tf = max(counts.values())
    for term, count in counts.items():
        col = vocab_index.get(term)
        if col is not None:
            vector[col] = (0.5 + 0.5 * count / max_tf) * idf[col]
    vector /= max(np.linalg.norm(vector), 1e-12)
    projected = vector @ components
    projected /= max(np.linalg.norm(projected), 1e-12)
    return projected, known_terms


def search(query: str, k: int = 5) -> list[dict[str, Any]]:
    if not query.strip():
        raise ValueError("Please provide a non-empty query.")
    if k < 1:
        raise ValueError("--k must be a positive integer.")

    vectors = np.load(VECTORS_PATH)
    metadata = json.loads(META_PATH.read_text(encoding="utf-8"))
    chunks = load_corpus()
    chunk_by_id = {chunk["chunk_id"]: chunk for chunk in chunks}
    ids = [chunk_id for chunk_id in metadata if chunk_id in chunk_by_id]
    if not ids:
        raise ValueError("The index has no matching records. Rebuild the index.")

    query_vector, known_terms = embed_query_lsa(query)
    if not known_terms:
        print("No query terms matched the indexed vocabulary; try different or more specific wording.")
        return []

    row_lookup = {chunk["chunk_id"]: index for index, chunk in enumerate(chunks)}
    ordered_ids = [chunk["chunk_id"] for chunk in chunks if chunk["chunk_id"] in metadata]
    row_indices = [row_lookup[chunk_id] for chunk_id in ordered_ids]
    if vectors.shape[0] != len(chunks):
        raise ValueError("Index/corpus size mismatch. Rebuild the index using --rebuild.")

    dense_scores = vectors[row_indices].astype(np.float64) @ query_vector
    lexical_scores = []
    for chunk_id in ordered_ids:
        doc_terms = set(tokenize(searchable_text(chunk_by_id[chunk_id])))
        lexical_scores.append(len(known_terms.intersection(doc_terms)) / max(len(known_terms), 1))
    lexical_array = np.asarray(lexical_scores, dtype=np.float64)
    combined_scores = 0.65 * dense_scores + 0.35 * lexical_array
    order = np.argsort(-combined_scores, kind="stable")[: min(k, len(ordered_ids))]

    results: list[dict[str, Any]] = []
    print(f"\nQuery: {query}")
    print("=" * 88)
    for rank, index in enumerate(order, start=1):
        chunk_id = ordered_ids[int(index)]
        meta = metadata[chunk_id]
        row = {
            "rank": rank,
            "chunk_id": chunk_id,
            "paper_id": meta["paper_id"],
            "title": meta["title"],
            "year": meta["year"],
            "venue": meta["venue"],
            "doi": meta["doi"],
            "source_url": meta.get("source_url", f"https://doi.org/{meta['doi']}"),
            "pmid": meta.get("pmid", ""),
            "record_type": meta.get("record_type", "not specified"),
            "snippet": (str(chunk_by_id[chunk_id].get("text", "")).strip()[:420] + ("…" if len(str(chunk_by_id[chunk_id].get("text", "")).strip()) > 420 else "")),
            "score": float(combined_scores[int(index)]),
            "lsa_similarity": float(dense_scores[int(index)]),
            "lexical_coverage": float(lexical_array[int(index)]),
        }
        results.append(row)
        print(
            f"{rank:2d}. [{row['score']:.3f}] {row['paper_id']} "
            f"({row['year']}, {row['venue']})"
        )
        print(f"     LSA={row['lsa_similarity']:.3f} | lexical={row['lexical_coverage']:.3f}")
        print(f"     {row['title']}")
        print(f"     {row['source_url']}\n")
    return results


def main(argv: list[str] | None = None) -> int:
    global CORPUS
    parser = argparse.ArgumentParser(
        description="Build or query the LYSIS starter literature-retrieval index."
    )
    parser.add_argument("--query", default=None, help="Research question or topic to retrieve.")
    parser.add_argument("--k", type=int, default=5, help="Number of results to display (default: 5).")
    parser.add_argument("--dim", type=int, default=128, help="Requested LSA dimension (default: 128).")
    parser.add_argument("--rebuild", action="store_true", help="Rebuild even if the index is current.")
    parser.add_argument("--corpus", type=Path, default=CORPUS, help="Path to JSONL corpus.")
    args = parser.parse_args(argv)

    CORPUS = args.corpus.resolve()
    try:
        if args.dim < 1:
            raise ValueError("--dim must be a positive integer.")
        if args.k < 1:
            raise ValueError("--k must be a positive integer.")
        if args.rebuild or not index_is_current():
            if args.query:
                print("Index missing or stale; building it from the current corpus first.")
            build_lsa(load_corpus(CORPUS), args.dim)
        if args.query:
            search(args.query, args.k)
        return 0
    except (OSError, ValueError, KeyError, json.JSONDecodeError, FloatingPointError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
