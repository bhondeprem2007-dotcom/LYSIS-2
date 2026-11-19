#!/usr/bin/env python3
"""LYSIS local research desk. Research support only; not a clinical tool."""
from __future__ import annotations
import csv, io, json, re, threading
from collections import Counter
from datetime import datetime, timezone
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from urllib.parse import urlparse, parse_qs
import build_embeddings as retrieval

ROOT = Path(__file__).resolve().parent
STATIC = ROOT / "static"
CORPUS_PATH = ROOT / "rag_corpus.jsonl"
HOST, PORT, MAX_BODY = "127.0.0.1", 8765, 2_000_000
retrieval.HERE = ROOT
retrieval.DEFAULT_CORPUS = CORPUS_PATH
retrieval.LOCAL_CORPUS = CORPUS_PATH
retrieval.CORPUS = CORPUS_PATH
retrieval.VECTORS_PATH = ROOT / "lsa_vectors.npy"
retrieval.VOCAB_PATH = ROOT / "lsa_vocab.json"
retrieval.META_PATH = ROOT / "chunk_meta.json"
retrieval.MANIFEST_PATH = ROOT / "lsa_index_manifest.json"
_INDEX_LOCK = threading.Lock()


def ensure_index():
    retrieval.CORPUS = CORPUS_PATH
    if not retrieval.index_is_current():
        with _INDEX_LOCK:
            if not retrieval.index_is_current():
                retrieval.build_lsa(retrieval.load_corpus(CORPUS_PATH))


def load_records():
    ensure_index()
    return retrieval.load_corpus(CORPUS_PATH)


def corpus_audit(records=None):
    records = records if records is not None else retrieval.load_corpus(CORPUS_PATH)
    ids, dois = set(), set(); duplicate_ids=[]; duplicate_dois=[]; missing=[]; bad_year=[]
    for r in records:
        pid=str(r.get('paper_id','')); doi=str(r.get('doi','')).strip().lower()
        if pid in ids: duplicate_ids.append(pid)
        if doi in dois: duplicate_dois.append(doi)
        ids.add(pid); dois.add(doi)
        absent=[k for k in ('title','year','venue','doi','text','source_url') if not str(r.get(k,'')).strip()]
        if absent: missing.append({'paper_id':pid,'fields':absent})
        if not (1900 <= int(r.get('year',0)) <= datetime.now().year+1): bad_year.append(pid)
    return {'records':len(records),'unique_paper_ids':len(ids),'unique_dois':len(dois),'duplicate_paper_ids':duplicate_ids,'duplicate_dois':duplicate_dois,'missing_metadata':missing,'year_outliers':bad_year,'status':'pass' if not (duplicate_ids or duplicate_dois or missing or bad_year) else 'review'}


def fasta_validate(text: str) -> dict:
    """Validate FASTA syntax and IUPAC nucleotide alphabet; no annotation."""
    if not isinstance(text, str) or not text.strip():
        return {'valid':False,'errors':['Paste a FASTA sequence first.'],'records':[],'scope':'Format validation only.'}
    if len(text)>1_500_000:
        return {'valid':False,'errors':['Input exceeds the 1.5 MB safety limit.'],'records':[],'scope':'Format validation only.'}
    lines=[x.strip() for x in text.replace('\r','').split('\n') if x.strip()]
    records,errors,current,seq_parts=[],[],None,[]
    allowed=set('ACGTRYSWKMBDHVN')
    def finish():
        nonlocal current,seq_parts
        if current is not None:
            seq=''.join(seq_parts).upper()
            invalid=sorted(set(seq)-allowed)
            if not seq: errors.append(f'{current}: sequence is empty.')
            if invalid: errors.append(f"{current}: unsupported nucleotide symbols: {', '.join(invalid)}.")
            if seq and not invalid:
                counts={b:seq.count(b) for b in 'ACGT'}
                gc=counts['G']+counts['C']
                records.append({'header':current,'length':len(seq),'gc_percent':round(100*gc/len(seq),2),'n_count':seq.count('N'),'ambiguous_count':sum(seq.count(x) for x in 'RYSWKMBDHVN'),'base_counts':counts,'at_percent':round(100*(counts['A']+counts['T'])/len(seq),2)})
        current,seq_parts=None,[]
    for n,line in enumerate(lines,1):
        if line.startswith('>'):
            finish(); current=line[1:].strip()
            if not current: errors.append(f'Line {n}: FASTA header is empty.'); current=f'unnamed_record_{len(records)+1}'
        elif current is None: errors.append(f'Line {n}: sequence data appears before the first FASTA header.')
        else: seq_parts.append(re.sub(r'\s+','',line))
    finish()
    if not records and not errors: errors.append('No valid FASTA records were found.')
    lengths=sorted((r['length'] for r in records),reverse=True)
    total=sum(lengths); running=0; n50=0
    for length in lengths:
        running+=length
        if running>=total/2: n50=length; break
    return {'valid':not errors and bool(records),'errors':errors,'records':records,'summary':{'record_count':len(records),'total_length':total,'longest':max(lengths,default=0),'shortest':min(lengths,default=0),'n50':n50},'scope':'FASTA syntax and nucleotide composition only. No genome annotation, gene calling, host prediction, phage matching, pathogenicity or treatment recommendation is performed.'}


def sequence_analysis(sequence: str, k: int = 3):
    """Descriptive nucleotide utilities; ORF scan is exploratory, not gene annotation."""
    seq=re.sub(r"\s+", "", str(sequence)).upper()
    allowed=set('ACGTRYSWKMBDHVN')
    if not seq: raise ValueError('Paste a nucleotide sequence first.')
    if len(seq)>100_000: raise ValueError('Sequence lab limit is 100,000 nt.')
    invalid=sorted(set(seq)-allowed)
    if invalid: raise ValueError('Unsupported nucleotide symbols: '+', '.join(invalid))
    if not 2<=int(k)<=6: raise ValueError('k-mer size must be between 2 and 6.')
    comp=str.maketrans('ACGTRYSWKMBDHVN','TGCAYRSWMKVHDBN')
    rc=seq.translate(comp)[::-1]
    kmers=Counter(seq[i:i+k] for i in range(max(0,len(seq)-k+1)))
    stops={'TAA','TAG','TGA'}; orfs=[]
    for strand, s in [('+',seq),('-',rc)]:
        for frame in range(3):
            start=None
            for pos in range(frame,len(s)-2,3):
                codon=s[pos:pos+3]
                if start is None and codon=='ATG': start=pos
                elif start is not None and codon in stops:
                    nt_len=pos+3-start
                    if nt_len>=90:
                        orfs.append({'strand':strand,'frame':frame+1,'start_1based':start+1,'end_1based':pos+3,'length_nt':nt_len})
                    start=None
    orfs.sort(key=lambda x:x['length_nt'],reverse=True)
    return {'length_nt':len(seq),'reverse_complement':rc,'gc_percent':round(100*(seq.count('G')+seq.count('C'))/len(seq),2),'kmer_size':k,'top_kmers':[{'kmer':a,'count':b} for a,b in kmers.most_common(20)],'orf_candidates':orfs[:30],'orf_scope':'Six-frame ATG-to-stop exploratory scan with a 90-nt minimum. These are ORF candidates only, not validated genes or phage functions.'}



CODON_TABLE = {
    'TTT':'F','TTC':'F','TTA':'L','TTG':'L','TCT':'S','TCC':'S','TCA':'S','TCG':'S',
    'TAT':'Y','TAC':'Y','TAA':'*','TAG':'*','TGT':'C','TGC':'C','TGA':'*','TGG':'W',
    'CTT':'L','CTC':'L','CTA':'L','CTG':'L','CCT':'P','CCC':'P','CCA':'P','CCG':'P',
    'CAT':'H','CAC':'H','CAA':'Q','CAG':'Q','CGT':'R','CGC':'R','CGA':'R','CGG':'R',
    'ATT':'I','ATC':'I','ATA':'I','ATG':'M','ACT':'T','ACC':'T','ACA':'T','ACG':'T',
    'AAT':'N','AAC':'N','AAA':'K','AAG':'K','AGT':'S','AGC':'S','AGA':'R','AGG':'R',
    'GTT':'V','GTC':'V','GTA':'V','GTG':'V','GCT':'A','GCC':'A','GCA':'A','GCG':'A',
    'GAT':'D','GAC':'D','GAA':'E','GAG':'E','GGT':'G','GGC':'G','GGA':'G','GGG':'G'
}
IUPAC = {'A':'A','C':'C','G':'G','T':'T','R':'AG','Y':'CT','S':'GC','W':'AT','K':'GT','M':'AC','B':'CGT','D':'AGT','H':'ACT','V':'ACG','N':'ACGT'}

def translate_sequence(sequence: str):
    seq=re.sub(r"\s+", "", str(sequence)).upper()
    if not seq: raise ValueError('Paste a nucleotide sequence first.')
    if len(seq)>100_000: raise ValueError('Translation limit is 100,000 nt.')
    invalid=sorted(set(seq)-set(IUPAC))
    if invalid: raise ValueError('Unsupported nucleotide symbols: '+', '.join(invalid))
    comp=str.maketrans('ACGTRYSWKMBDHVN','TGCAYRSWMKVHDBN')
    rc=seq.translate(comp)[::-1]
    frames=[]
    for strand,label in [(seq,'+'),(rc,'-')]:
        for frame in range(3):
            aa=''.join(CODON_TABLE.get(strand[i:i+3],'X') if len(strand[i:i+3])==3 else '' for i in range(frame,len(strand)-2,3))
            frames.append({'strand':label,'frame':frame+1,'amino_acids':aa,'aa_length':len(aa),'complete_codons':len(aa),'ambiguous_codons':aa.count('X'),'stop_count':aa.count('*')})
    return {'length_nt':len(seq),'frames':frames,'scope':'Standard genetic code, six-frame translation preview. Ambiguous codons are X. This is not gene annotation, protein function prediction or clinical interpretation.'}

def scan_motif(sequence: str, motif: str):
    seq=re.sub(r"\s+", "", str(sequence)).upper()
    pat=re.sub(r"\s+", "", str(motif)).upper()
    if not seq or not pat: raise ValueError('Enter both a nucleotide sequence and a motif.')
    if len(seq)>100_000: raise ValueError('Motif scan limit is 100,000 nt.')
    if len(pat)>100: raise ValueError('Motif length limit is 100 symbols.')
    bad_seq=sorted(set(seq)-set(IUPAC)); bad_pat=sorted(set(pat)-set(IUPAC))
    if bad_seq: raise ValueError('Unsupported sequence symbols: '+', '.join(bad_seq))
    if bad_pat: raise ValueError('Unsupported motif symbols: '+', '.join(bad_pat))
    import re as _re
    regex=''.join('['+IUPAC[c]+']' for c in pat)
    comp=str.maketrans('ACGTRYSWKMBDHVN','TGCAYRSWMKVHDBN')
    rc_pat=pat.translate(comp)[::-1]
    hits=[]
    for strand, query in [('+',pat),('-',rc_pat)]:
        rx=_re.compile('(?=('+''.join('['+IUPAC[c]+']' for c in query)+'))')
        for m in rx.finditer(seq):
            i=m.start(); start=i+1; end=i+len(pat)
            hits.append({'strand':strand,'start_1based':start,'end_1based':end,'matched_sequence':seq[i:i+len(pat)],'query_motif':pat})
    hits.sort(key=lambda h:(h['start_1based'],h['strand']))
    return {'sequence_length':len(seq),'motif':pat,'reverse_complement_motif':rc_pat,'hit_count':len(hits),'hits':hits[:500],'truncated':len(hits)>500,'scope':'IUPAC pattern matching on the supplied sequence in both orientations. Matches are sequence-level observations only; no restriction-enzyme activity, biological function or clinical meaning is inferred.'}

def citation(r):
    authors=r.get('authors') or 'Author metadata not supplied'
    return f"{authors} ({r.get('year','n.d.')}). {r.get('title','Untitled')}. {r.get('venue','Source not supplied')}. https://doi.org/{r.get('doi','')}"


def brief_for(query, results, records):
    by_id={r['paper_id']:r for r in records}
    rows=[]
    for hit in results:
        src=by_id.get(hit['paper_id'],{})
        rows.append({'paper_id':hit['paper_id'],'title':hit['title'],'year':hit['year'],'venue':hit['venue'],'doi':hit['doi'],'source_url':hit['source_url'],'ranking_score':round(hit['score'],4),'research_note':src.get('text',''),'citation':citation(src)})
    return {'query':query,'created_utc':datetime.now(timezone.utc).isoformat(),'records':rows,'caveat':'This is a retrieval-based evidence brief, not a generated systematic review, clinical recommendation, or assessment of treatment efficacy. Verify each original publication.'}


class Handler(BaseHTTPRequestHandler):
    server_version='LYSISResearchDesk/4.0'
    def send_payload(self,status,payload,ctype='application/json; charset=utf-8'):
        data=payload if isinstance(payload,bytes) else (payload.encode() if isinstance(payload,str) else json.dumps(payload,ensure_ascii=False).encode())
        self.send_response(status); self.send_header('Content-Type',ctype); self.send_header('Content-Length',str(len(data)))
        self.send_header('Cache-Control','no-store'); self.send_header('X-Content-Type-Options','nosniff'); self.send_header('X-Frame-Options','DENY'); self.send_header('Referrer-Policy','no-referrer'); self.send_header('Content-Security-Policy',"default-src 'self' 'unsafe-inline'; connect-src 'self'; img-src 'self' data:; base-uri 'none'; form-action 'self'"); self.end_headers(); self.wfile.write(data)
    def read_body(self):
        try: length=int(self.headers.get('Content-Length','0'))
        except ValueError: raise ValueError('Invalid Content-Length header.')
        if length<0 or length>MAX_BODY: raise OverflowError('Request too large (2 MB maximum).')
        raw=self.rfile.read(length)
        body=json.loads(raw or b'{}')
        if not isinstance(body,dict): raise ValueError('Request body must be a JSON object.')
        return body
    def do_GET(self):
        parsed=urlparse(self.path); path=parsed.path
        try:
            if path in ('/','/index.html'):
                return self.send_payload(200,(STATIC/'index.html').read_text(encoding='utf-8'),'text/html; charset=utf-8')
            if path=='/api/health':
                records=load_records(); manifest=json.loads(retrieval.MANIFEST_PATH.read_text())
                return self.send_payload(200,{'ok':True,'records':len(records),'index_dimension':manifest['actual_dimension'],'version':retrieval.VERSION,'mode':'local LSA + lexical retrieval','modules':['literature search','paper explorer','compare papers','evidence brief','reading list','FASTA analysis','sequence lab','six-frame translation','IUPAC motif scan','corpus audit','catalogue export']})
            if path=='/api/papers':
                records=load_records(); q=parse_qs(parsed.query).get('q',[''])[0].lower().strip()
                found=[r for r in records if not q or q in ' '.join([r.get('title',''),r.get('paper_id',''),r.get('venue',''),' '.join(r.get('areas',[])),r.get('text','')]).lower()]
                payload=[{k:r.get(k,'') for k in ('paper_id','title','year','venue','doi','source_url','authors','pmid','areas','record_type')} for r in found]
                return self.send_payload(200,{'count':len(payload),'papers':payload})
            if path=='/api/paper':
                pid=parse_qs(parsed.query).get('id',[''])[0]
                r=next((x for x in load_records() if x['paper_id']==pid),None)
                return self.send_payload(200 if r else 404,r or {'error':'Paper not found.'})
            if path=='/api/audit':
                return self.send_payload(200,corpus_audit(load_records()))
            return self.send_payload(404,{'error':'Not found'})
        except Exception as e:
            print(f'[LYSIS] GET error: {type(e).__name__}: {e}')
            return self.send_payload(500,{'error':'Could not complete request. Check terminal output.'})
    def do_POST(self):
        try:
            body=self.read_body(); path=urlparse(self.path).path
            if path=='/api/search':
                query=str(body.get('query','')).strip(); k=int(body.get('k',5))
                if not query: return self.send_payload(400,{'error':'Enter a research question.'})
                if len(query)>1000: return self.send_payload(400,{'error':'Query must be 1,000 characters or fewer.'})
                if not 1<=k<=20: return self.send_payload(400,{'error':'Choose between 1 and 20 results.'})
                records=load_records(); results=retrieval.search(query,k)
                return self.send_payload(200,{'query':query,'count':len(results),'results':results,'scoring':'0.65 × LSA cosine similarity + 0.35 × query-term coverage','warning':'Scores are retrieval ranking values, not probabilities, clinical confidence or treatment efficacy.','corpus_records':len(records)})
            if path=='/api/fasta': return self.send_payload(200,fasta_validate(str(body.get('fasta',''))))
            if path=='/api/sequence-analysis': return self.send_payload(200,sequence_analysis(str(body.get('sequence','')),int(body.get('k',3))))
            if path=='/api/translate': return self.send_payload(200,translate_sequence(str(body.get('sequence',''))))
            if path=='/api/motif-scan': return self.send_payload(200,scan_motif(str(body.get('sequence','')),str(body.get('motif',''))))
            if path=='/api/brief':
                query=str(body.get('query','')).strip(); k=max(1,min(10,int(body.get('k',5))))
                if not query:return self.send_payload(400,{'error':'Enter a topic for the evidence brief.'})
                records=load_records(); hits=retrieval.search(query,k)
                return self.send_payload(200,brief_for(query,hits,records))
            if path=='/api/compare':
                ids=body.get('paper_ids',[])
                if not isinstance(ids,list) or not 2<=len(ids)<=4:return self.send_payload(400,{'error':'Select 2 to 4 paper IDs to compare.'})
                records=load_records(); by={r['paper_id']:r for r in records}
                if len(set(ids))!=len(ids) or any(x not in by for x in ids):return self.send_payload(400,{'error':'Paper selection contains duplicates or unknown IDs.'})
                fields=('paper_id','title','year','venue','authors','doi','areas','record_type','text')
                return self.send_payload(200,{'papers':[{k:by[x].get(k,'') for k in fields} for x in ids],'note':'Side-by-side metadata and source-note comparison only; no clinical suitability ranking.'})
            return self.send_payload(404,{'error':'Not found'})
        except OverflowError as e:return self.send_payload(413,{'error':str(e)})
        except (ValueError,TypeError,json.JSONDecodeError) as e:return self.send_payload(400,{'error':str(e) or 'Invalid request.'})
        except Exception as e:
            print(f'[LYSIS] internal error: {type(e).__name__}: {e}')
            return self.send_payload(500,{'error':'The prototype could not complete this request. Check the terminal for details.'})
    def log_message(self,fmt,*args):print('[LYSIS]',fmt%args)


def main():
    ensure_index(); print(f'\nLYSIS Research Desk v4.0: http://{HOST}:{PORT}\nLocal-only prototype. Press Ctrl+C to stop.')
    server=ThreadingHTTPServer((HOST,PORT),Handler); server.daemon_threads=True
    try: server.serve_forever()
    except KeyboardInterrupt: print('\nLYSIS stopped.')
    finally: server.server_close()
if __name__=='__main__': main()
