import argparse,json,time,copy,pathlib,Bio.SeqIO,hashlib,collections.Counter
from Bio import SeqIO
import matching

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('corpus_directory')
    parser.add_argument('database_json')
    parser.add_argument('--limit',default=20,type=int)
    parser.add_argument('--output',type=str,default=None)
    args=parser.parse_args()

    start_time=time.time()
    try:
        with open(args.database_json,'r') as f:db=json.load(f)
    except Exception as e:print(json.dumps({"error":str(e)}));return
    ref_sha256={e['source_sha256']:e for e in db}
    corpus_dir=pathlib.Path(args.corpus_directory)
    files=list(corpus_dir.rglob('*.dna'))
    groups={}
    for f in files:groups[f.parent]=f
    rounds=0;removed_count=0;errors=[];cases=0
    expected_counter=collections.Counter()
    actual_counter=collections.Counter()

    for parent in sorted(groups.keys()):
        if rounds>=args.limit:break
        group_files=groups[parent]
        for f in group_files:
            try:
                seq=SeqIO.read(f,'snapgene')
                canonical_map={e['canonical_name'].casefold():e['canonical_name'] for e in db}
                targets=[e for e in db if e.get('features') and any(e['features'][i]['label'].casefold()==canonical_map.get(k,k).casefold() for k in e['features'][i].keys())]
                if not targets:continue
                original=copy.deepcopy(seq)
                indices=[i for i,feat in enumerate(original.features) if feat.id in [t['id'] for t in targets]]
                new_features=[f for i,f in enumerate(original.features) if i not in indices]
                masked=original.replace(features=new_features)
                hits=matching.match_features(masked,db)
                status_new=[h for h in hits if h['status']=='new' and h.get('confidence')=='High']
                if not status_new:continue
                cases+=1
                removed_count+=len(status_new)
                for hit in status_new:
                    entry=next((e for e in db if e['canonical_name']==hit['canonical_name']),None)
                    if entry:
                        loc=tuple(sorted(set(hit['location'])))
                        actual_counter[(entry['canonical_name'].casefold(),hit['strand'],loc)]+=1
                expected_counter.update({(t['canonical_name'].casefold(),t['strand'],tuple(sorted(set(t['location'])))) for t in targets})
            except Exception as e:errors.append((str(f),str(e)));continue
        rounds+=1

    tp=sum((expected_counter & actual_counter).values())
    fp=sum((actual_counter - expected_counter).values())
    fn=sum((expected_counter - actual_counter).values())
    precision=tp/(tp+fp) if (tp+fp)>0 else None
    recall=tp/(tp+fn) if (tp+fn)>0 else None
    recovery=recall

    report={
        "scope":"Agreement with supplied annotations on source-excluded files; not independent biological validation",
        "reference_count":len(db),
        "excluded_reference_source_files_count":sum(1 for e in db if e['source_sha256'] in ref_sha256),
        "cases":cases,
        "removed":removed_count,
        "tp":int(tp),"fp":int(fp),"fn":int(fn),
        "precision":precision,"recall":recall,"recovery":recovery,
        "errors":errors,"time_seconds":round(time.time()-start_time,3)
    }

    if args.output:
        with open(args.output,'w',encoding='utf8') as f:json.dump(report,f,indent=2)
    else:print(json.dumps(report))

if __name__=='__main__':main()