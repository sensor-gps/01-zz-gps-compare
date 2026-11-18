"""Explicit fusion frame conversion and boot timestamp -> GPS iTOW mapping."""
import numpy as np

def load(path,anchor_boot='',anchor_tow='',frame='local',north_confirmed=False,gps1=None):
    from pyulog import ULog
    topics={t.name:t.data for t in ULog(path).data_list}
    if 'drone_fusion' not in topics: raise ValueError('ULog 缺少 drone_fusion')
    f=topics['drone_fusion']; t=np.asarray(f['timestamp'],dtype=float)
    gps=topics.get('sensor_gps',{})
    key=next((k for k in ('iTOW','itow','i_tow','gps_itow','time_of_week') if k in gps),None)
    if gps1 is not None:
        # Match receiver values, not truth motion: preserve residual GPS latency.
        required=('timestamp','veln','vele','veld','lat','lon')
        if any(k not in gps for k in required): raise ValueError('sensor_gps lacks clock-anchor fields')
        signature={}
        for tow,row in gps1.items():
            sig=tuple(round(row[k]*1000) for k in ('vel_n','vel_e','vel_d'))
            signature.setdefault(sig,[]).append((tow,row))
        offsets=[]
        for i,stamp in enumerate(gps['timestamp']):
            sig=tuple(round(float(gps[k][i])*1000) for k in ('veln','vele','veld'))
            matches=signature.get(sig,[])
            lat=float(gps['lat'][i]); lon=float(gps['lon'][i])
            if abs(lat)>90: lat/=1e7
            if abs(lon)>180: lon/=1e7
            matches=[(tow,row) for tow,row in matches if abs(row['lat']-lat)<2e-7 and abs(row['lon']-lon)<2e-7]
            if len(matches)==1: offsets.append(matches[0][0]-float(stamp)/1000)
        if len(offsets)<20: raise ValueError('sensor_gps / GPS1 identity matches insufficient; select the receiver recorded by this ULog')
        offsets=np.asarray(offsets); center=np.median(offsets)
        accepted=offsets[np.abs(offsets-center)<200]
        if len(accepted)<20 or len(accepted)<len(offsets)*.8: raise ValueError('sensor_gps / GPS1 clock mapping is inconsistent')
        mapped=t/1000+np.median(accepted)
        inside=np.ones(len(t),dtype=bool)
    elif key:
        gt=np.asarray(gps['timestamp'],dtype=float)
        tow=np.asarray(gps[key],dtype=float)
        tow=np.unwrap(tow*2*np.pi/604800000)*604800000/(2*np.pi)
        valid=np.isfinite(tow)&(tow>0)
        gt,tow=gt[valid],tow[valid]
        if len(gt)<2: raise ValueError('ULog GPS iTOW 有效锚点不足')
        mapped=np.interp(t,gt,tow)
        inside=(t>=gt[0])&(t<=gt[-1])
    else:
        if str(anchor_boot).strip() and str(anchor_tow).strip():
            boot=float(anchor_boot); tow=float(anchor_tow)
            if not np.isfinite([boot,tow]).all() or not 0<tow<604800000: raise ValueError('时间锚点不合法')
            mapped=tow+(t-boot)/1000
        else:
            # 参考 compare_dual_gps.py：先保留飞机相对时间，外层使用速度
            # 相关性把它与 UBX 的 iTOW 时钟对齐，不要求用户手填锚点。
            mapped=t/1000.0
        inside=np.ones(len(t),dtype=bool)
    g=topics.get('fusion_gps')
    if g is None or 'yaw' not in g: raise ValueError('缺少 fusion_gps.yaw，无法建立真北参考')
    gt=np.asarray(g['timestamp'],dtype=float)
    fy=np.unwrap(np.deg2rad(g['yaw']))
    yaw_gps=np.interp(t,gt,fy)
    yaw_ref=np.unwrap(np.deg2rad(f.get('yaw',np.zeros(len(t)))))
    inside &= (t>=gt[0])&(t<=gt[-1])
    v=np.column_stack([f['velocity_x'],f['velocity_y'],f['velocity_z']]).astype(float)
    p=np.column_stack([f.get('position_x',np.zeros(len(t))),f.get('position_y',np.zeros(len(t))),f.get('position_z',np.zeros(len(t))) ]).astype(float)
    # Reference implementation rotates by negative(drone_yaw - gps_yaw).
    yaw_diff=np.arctan2(np.sin(yaw_ref-yaw_gps),np.cos(yaw_ref-yaw_gps))
    angle=-yaw_diff
    c,s=np.cos(angle),np.sin(angle)
    pn=c*p[:,0]-s*p[:,1]; pe=s*p[:,0]+c*p[:,1]; pd=p[:,2]
    vn=c*v[:,0]-s*v[:,1]; ve=s*v[:,0]+c*v[:,1]; vd=v[:,2]
    valid=inside & np.isfinite(vn+ve+vd+mapped)
    if 'is_hori_valid' in f: valid &= np.asarray(f['is_hori_valid'])!=0
    if 'is_vert_valid' in f: valid &= np.asarray(f['is_vert_valid'])!=0
    if 'using_fake_fusion' in f: valid &= np.asarray(f['using_fake_fusion'])==0
    result={int(round(mapped[i])):{'pos_n':pn[i],'pos_e':pe[i],'pos_d':pd[i],
             'vel_n':vn[i],'vel_e':ve[i],'vel_d':vd[i],
             'yaw_diff_deg':float(np.rad2deg(yaw_diff[i])), 'speed':float(np.hypot(vn[i],ve[i]))} for i in np.flatnonzero(valid)}
    if not result: raise ValueError('转换后没有有效 fusion 数据')
    return result
