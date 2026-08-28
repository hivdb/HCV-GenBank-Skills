#!/usr/bin/env python3
"""Assign NS3/NS5A/NS5B genotype, then genotype-matched subtype, from FASTAs."""
from __future__ import annotations
import argparse,csv,json,re,subprocess,tempfile
from collections import defaultdict
from pathlib import Path

GENES=("NS3","NS5A","NS5B"); EXT={".fa",".fasta",".fna"}
OUTFMT="6 qseqid sseqid length mismatch gaps pident evalue bitscore"
def fasta(path):
 h=None;c=[]
 for line in path.read_text().splitlines():
  if line.startswith('>'):
   if h is not None: yield h,''.join(c).upper()
   h=line[1:].strip();c=[]
  elif line.strip(): c.append(re.sub(r'\s+','',line))
 if h is not None: yield h,''.join(c).upper()
def write(path, rows):
 with path.open('w') as o:
  for h,s in rows:o.write(f'>{h}\n{s}\n')
def blast(query, db, output, threads):
 subprocess.run(['blastn','-query',str(query),'-db',str(db),'-dust','no','-task','blastn','-num_threads',str(threads),'-evalue','1e-6','-max_hsps','1','-max_target_seqs','1000','-outfmt',OUTFMT,'-out',str(output)],check=True,stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)
 return [x.split('\t') for x in output.read_text().splitlines() if x]
def db(path,prefix):
 subprocess.run(['makeblastdb','-in',str(path),'-dbtype','nucl','-out',str(prefix),'-parse_seqids'],check=True,stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)
def main():
 p=argparse.ArgumentParser();p.add_argument('--fasta-dir',required=True);p.add_argument('--output-dir',default='archived-skills/outputs/folder_assignments');p.add_argument('--gt-reference-fasta',default='HCVData/HCV_GT_RefSeqs.fasta');p.add_argument('--subtype-json',default='HCVData/HCV_Subtype_Refs_By_Genome_NA.json');p.add_argument('--min-aligned-nt',type=int,default=200);p.add_argument('--threads',type=int,default=1);a=p.parse_args()
 if a.threads < 1: raise SystemExit('--threads must be at least 1')
 out=Path(a.output_dir);out.mkdir(parents=True,exist_ok=True); records=[]
 for f in sorted(Path(a.fasta_dir).rglob('*')):
  if f.is_file() and f.suffix.lower() in EXT:
   for h,s in fasta(f): records.append((str(len(records)),h.split()[0],h,str(f),s))
 if not records: raise SystemExit('No FASTA records found')
 refs=[]
 for h,s in fasta(Path(a.gt_reference_fasta)):
  m=re.match(r'HCV([1-8])(NS3|NS5A|NS5B)',h.split()[0])
  if m: refs.append((f'{m.group(2)}|{m.group(1)}',s))
 with tempfile.TemporaryDirectory() as d:
  d=Path(d); q=d/'q.fa';write(q,[(r[0],r[4]) for r in records]);rf=d/'gt.fa';write(rf,refs);gdb=d/'gt';db(rf,gdb);hits=blast(q,gdb,d/'gt.tsv',a.threads)
  best={}
  for x in hits:
   gene,gt=x[1].rsplit('|',2)[-2:]; key=(x[0],gene); score=(float(x[7]),int(x[2]))
   if int(x[2])>=a.min_aligned_nt and (key not in best or score>best[key][0]):best[key]=(score,gt,x)
  subtype=json.loads(Path(a.subtype_json).read_text()); bygt=defaultdict(list)
  for r in subtype:
   m=re.search(r'Genotype\s*([1-8][A-Za-z0-9]*)',str(r.get('genotypeName','')))
   if m and r.get('sequence'):bygt[m.group(1)[0]].append((m.group(1),str(r.get('accession','')),str(r['sequence']).upper()))
  result=defaultdict(list)
  for gene in GENES:
   groups=defaultdict(list)
   for (qid,g),(_,gt,x) in best.items():
    if g==gene:groups[gt].append(qid)
   for gt,qids in groups.items():
    sf=d/f'{gene}_{gt}.fa'; meta={f'S{i}':v for i,v in enumerate(bygt[gt])};write(sf,[(k,v[2]) for k,v in meta.items()]);sdb=d/f'{gene}_{gt}';db(sf,sdb);qq=d/f'{gene}_{gt}_q.fa';write(qq,[(q,records[int(q)][4]) for q in qids]);sh=blast(qq,sdb,d/f'{gene}_{gt}.tsv',a.threads); sb={}
    for x in sh:
     score=(float(x[7]),int(x[2]));
     x[1]=x[1].rsplit('|',1)[-1]
     if int(x[2])>=a.min_aligned_nt and (x[0] not in sb or score>sb[x[0]][0]):sb[x[0]]=(score,x)
    for qid in qids:
     _,acc,hdr,src,_=records[int(qid)];_,_,gh=best[(qid,gene)]; row={'accession':acc,'header':hdr,'source_fasta':src,'gene':gene,'genotype':gt,'genotype_aligned_nt':gh[2],'genotype_pident':gh[5],'subtype':'','subtype_reference_accession':'','subtype_aligned_nt':'','subtype_pident':''}
     if qid in sb:
      x=sb[qid][1]; st,sa,_=meta[x[1]];row.update(subtype=st,subtype_reference_accession=sa,subtype_aligned_nt=x[2],subtype_pident=x[5])
     result[gene].append(row)
  fields=['accession','header','source_fasta','gene','genotype','genotype_aligned_nt','genotype_pident','subtype','subtype_reference_accession','subtype_aligned_nt','subtype_pident']
  for gene in GENES:
   with (out/f'{gene}_assignments.csv').open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(result[gene])
  (out/'assignment_summary.json').write_text(json.dumps({'input_records':len(records),'assignments_by_gene':{g:len(result[g]) for g in GENES}},indent=2)+'\n')
if __name__=='__main__':main()
