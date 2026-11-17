import numpy as np

def estimate_delay(reference, gps, max_ms=2000, step_ms=10):
    """???? GPS ???????"""
    rt=np.array(sorted(reference),float); rv=np.array([[reference[t].get('vel_n',np.nan),reference[t].get('vel_e',np.nan)] for t in rt])
    gt=np.array(sorted(gps),float); gv=np.array([[gps[t].get('vel_n',np.nan),gps[t].get('vel_e',np.nan)] for t in gt])
    gvok=np.isfinite(gv).all(1); rvok=np.isfinite(rv).all(1)
    rt,rv=rt[rvok],rv[rvok]; gt,gv=gt[gvok],gv[gvok]
    scores=[]
    for delay in range(-max_ms,max_ms+1,step_ms):
        pred=np.column_stack([np.interp(rt+delay,gt,gv[:,k],left=np.nan,right=np.nan) for k in range(2)])
        ok=np.isfinite(pred).all(1)
        if ok.sum()<20: scores.append((np.inf,delay)); continue
        a=rv[ok]-rv[ok].mean(0); b=pred[ok]-pred[ok].mean(0)
        scores.append((-float(np.sum(a*b)/(np.linalg.norm(a)*np.linalg.norm(b)+1e-12)),delay))
    score,delay=min(scores); return float(delay),float(-score)

def estimate_clock_offset(reference, gps, max_s=120.0, step_s=0.25):
    """Use horizontal speed signatures to find the large FFC/GPS clock offset."""
    rt=np.array(sorted(reference),float); rs=np.array([np.hypot(reference[t].get('vel_n',np.nan),reference[t].get('vel_e',np.nan)) for t in rt])
    gt=np.array(sorted(gps),float); gs=np.array([np.hypot(gps[t].get('vel_n',np.nan),gps[t].get('vel_e',np.nan)) for t in gt])
    goodr=np.isfinite(rs); goodg=np.isfinite(gs)
    rt,rs=rt[goodr],rs[goodr]; gt,gs=gt[goodg],gs[goodg]
    if len(rt)<20 or len(gt)<20: return float(min(gt)-min(rt))
    grid=np.arange(max(min(rt),min(gt)-max_s*1000),min(max(rt),max(gt)+max_s*1000)+1,100.0)
    best=(-np.inf,0.0)
    for off in np.arange(-max_s,max_s+step_s/2,step_s):
        # FFC time + off seconds is compared with GPS time.
        q=np.interp(grid-off*1000,rt,rs,left=np.nan,right=np.nan); g=np.interp(grid,gt,gs,left=np.nan,right=np.nan)
        ok=np.isfinite(q)&np.isfinite(g)
        if ok.sum()<30: continue
        a=q[ok]-q[ok].mean(); b=g[ok]-g[ok].mean(); score=float(np.sum(a*b)/(np.linalg.norm(a)*np.linalg.norm(b)+1e-12))
        if score>best[0]: best=(score,off)
    return float(min(gt)-min(rt)+best[1]*1000)
