"""无真值局部去趋势残差比较；不把残差解释为可辨识的绝对噪声。"""
from pathlib import Path
from datetime import datetime
import csv
import numpy as np
from .accuracy_analysis import static_errors

HEADERS=('指标','单位','GPS1 残差 RMS','GPS2 残差 RMS','改善比例','结论','差分残差 RMS','有效点数')


def detrend(times, values, window):
    """完整居中时间窗二次拟合，至少11点，不跨超过0.5秒空洞。"""
    result=np.full_like(values,np.nan,dtype=float)
    for i,t in enumerate(times):
        lo=np.searchsorted(times,t-window/2);hi=np.searchsorted(times,t+window/2,side='right')
        if hi-lo<11 or t-times[lo]<window*.4 or times[hi-1]-t<window*.4:continue
        if np.max(np.diff(times[lo:hi]))>.5:continue
        x=(times[lo:hi]-t)/(window/2)
        design=np.column_stack((np.ones(len(x)),x,x*x))
        coef=np.linalg.lstsq(design,values[lo:hi],rcond=None)[0]
        result[i]=values[i]-coef[0]
    return result


def analyze(first,second,window=3.,tolerance_ms=60.):
    if not np.isfinite(window) or not 1<=window<=30:raise ValueError('滑窗长度必须在 1～30 秒之间')
    if not np.isfinite(tolerance_ms) or not 0<=tolerance_ms<=500:raise ValueError('匹配容差必须在 0～500 ms 之间')
    a=np.array(sorted(first));b=np.array(sorted(second));pairs=[];used=set()
    # 一对一最近邻，避免低频接收机的同一记录重复参与统计。
    for t in a:
        j=int(np.searchsorted(b,t));options=[k for k in (j-1,j) if 0<=k<len(b) and k not in used]
        if not options:continue
        k=min(options,key=lambda k:abs(int(b[k])-int(t)))
        if abs(int(b[k])-int(t))<=tolerance_ms:pairs.append((int(t),int(b[k])));used.add(k)
    if len(pairs)<30:raise ValueError('一对一共同匹配点不足 30 个，请检查日志或容差')
    ta=np.array([p[0] for p in pairs])/1000.;tb=np.array([p[1] for p in pairs])/1000.
    rows_a={i:first[t] for i,(t,_) in enumerate(pairs)};rows_b={i:second[t] for i,(_,t) in enumerate(pairs)}
    origin=rows_a[0]
    def values(rows):
        position=static_errors(rows,origin['lat'],origin['lon'],origin['alt'])
        return np.column_stack(([ [r[k] for r in rows.values()] for k in ('vel_n','vel_e','vel_d')]+[position[k] for k in ('位置北向','位置东向','位置地向')]))
    va,vb=values(rows_a),values(rows_b)
    ra,rb=detrend(ta,va,window),detrend(tb,vb,window)
    good=np.isfinite(ra).all(1)&np.isfinite(rb).all(1)
    if good.sum()<20:raise ValueError('完整滑窗内有效点不足 20 个，请缩短窗口或检查数据空洞')
    ra,rb=ra[good],rb[good];out=[]
    for k,title in enumerate(('速度 N','速度 E','速度 D','位置 N','位置 E','位置 D')):
        x,y=np.sqrt(np.mean(ra[:,k]**2)),np.sqrt(np.mean(rb[:,k]**2))
        equal=np.isclose(x,y,rtol=1e-6,atol=1e-9)
        tag='neutral' if equal else ('better' if y<x else 'worse')
        ratio=f'{(x-y)/x*100:+.2f}%' if x>1e-12 else '—'
        out.append(((title,'m/s' if k<3 else 'm',f'{x:.6g}',f'{y:.6g}',ratio,
                     '相同' if equal else ('GPS2 更平滑' if y<x else 'GPS2 波动更大'),
                     f'{np.sqrt(np.mean((rb[:,k]-ra[:,k])**2)):.6g}',int(good.sum())),tag))
    directory=Path(__file__).resolve().parent.parent/'analysis_output';directory.mkdir(exist_ok=True)
    stem='relative_accuracy_'+datetime.now().strftime('%Y%m%d_%H%M%S_%f')
    path=directory/(stem+'.csv')
    with path.open('w',newline='',encoding='utf-8-sig') as stream:
        writer=csv.writer(stream);writer.writerow(HEADERS);writer.writerows(row for row,_ in out)
    detail=directory/(stem+'_residuals.csv')
    with detail.open('w',newline='',encoding='utf-8-sig') as stream:
        writer=csv.writer(stream);writer.writerow(['gps1_itow_ms','gps2_itow_ms']+[f'g{g}_{k}' for g in (1,2) for k in ('vn','ve','vd','pn','pe','pd')])
        for i,index in enumerate(np.flatnonzero(good)):writer.writerow([*pairs[index],*ra[i],*rb[i]])
    summary=f'fix=3 唯一 PVT：GPS1 {len(first)}，GPS2 {len(second)}；一对一匹配 {len(pairs)}，完整窗口统计 {good.sum()} 点。\nGPS1 统计 iTOW：{ta[good][0]*1000:.0f}～{ta[good][-1]*1000:.0f} ms；GPS2：{tb[good][0]*1000:.0f}～{tb[good][-1]*1000:.0f} ms。\n窗口 {window:g} 秒，二次趋势，容差 {tolerance_ms:g} ms。未自动补偿相对时延。\n汇总：{path}\n逐点残差：{detail}\n'
    fields=[]
    for k,title in enumerate(('速度 N','速度 E','速度 D','位置 N','位置 E','位置 D')):
        unit='m/s' if k<3 else 'm'
        fields.extend([(f'raw{k}',title,unit),(f'residual{k}',title+' 去趋势残差',unit)])
    aligned=[]
    for times,values,residual in ((ta,va,detrend(ta,va,window)),(tb,vb,detrend(tb,vb,window))):
        rows=[]
        for index in np.flatnonzero(good):
            row={f'raw{k}':values[index,k] for k in range(6)}
            row.update({f'residual{k}':residual[index,k] for k in range(6)})
            rows.append((int(round(times[index]*1000)),row))
        aligned.append(rows)
    return out,summary,({t:{} for t,_ in aligned[0]},aligned,fields)
