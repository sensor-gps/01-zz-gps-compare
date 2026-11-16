"""Explicit fusion frame conversion and boot timestamp -> GPS iTOW mapping."""
import numpy as np

def load(path,anchor_boot='',anchor_tow='',frame='body',north_confirmed=False):
    from pyulog import ULog
    topics={t.name:t.data for t in ULog(path).data_list}
    if 'drone_fusion' not in topics: raise ValueError('ULog 缺少 drone_fusion')
    f=topics['drone_fusion']; t=np.asarray(f['timestamp'],dtype=float)
    gps=topics.get('sensor_gps',{})
    key=next((k for k in ('iTOW','itow','i_tow','gps_itow','time_of_week') if k in gps),None)
    if key:
        gt=np.asarray(gps['timestamp'],dtype=float)
        tow=np.asarray(gps[key],dtype=float)
        tow=np.unwrap(tow*2*np.pi/604800000)*604800000/(2*np.pi)
        valid=np.isfinite(tow)&(tow>0)
        gt,tow=gt[valid],tow[valid]
        if len(gt)<2: raise ValueError('ULog GPS iTOW 有效锚点不足')
        mapped=np.interp(t,gt,tow)
        inside=(t>=gt[0])&(t<=gt[-1])
    else:
        if not str(anchor_boot).strip() or not str(anchor_tow).strip():
            raise ValueError('此 ULog 的 sensor_gps 不含 iTOW。请填写同一时刻的 ULog 时间戳（微秒）及 GPS iTOW（毫秒），不能直接把启动时间当作 iTOW。')
        boot=float(anchor_boot); tow=float(anchor_tow)
        if not np.isfinite([boot,tow]).all() or not 0<tow<604800000: raise ValueError('时间锚点不合法')
        mapped=tow+(t-boot)/1000; inside=np.ones(len(t),dtype=bool)
    if not north_confirmed: raise ValueError('请先确认 fusion_gps.yaw 的定义为真北航向，再勾选确认。若它是磁北或任意局部航向，需要先获取真北转换关系。')
    g=topics.get('fusion_gps')
    if g is None or 'yaw' not in g: raise ValueError('缺少 fusion_gps.yaw，无法建立真北参考')
    gt=np.asarray(g['timestamp'],dtype=float)
    yaw=np.interp(t,gt,np.unwrap(np.deg2rad(g['yaw'])))
    inside &= (t>=gt[0])&(t<=gt[-1])
    v=np.column_stack([f['velocity_x'],f['velocity_y'],f['velocity_z']]).astype(float)
    if frame=='body':
        roll=np.deg2rad(f['roll']); pitch=np.deg2rad(f['pitch'])
        cr,sr=np.cos(roll),np.sin(roll); cp,sp=np.cos(pitch),np.sin(pitch)
        cy,sy=np.cos(yaw),np.sin(yaw); x,y,z=v.T
        vn=cp*cy*x+(sr*sp*cy-cr*sy)*y+(cr*sp*cy+sr*sy)*z
        ve=cp*sy*x+(sr*sp*sy+cr*cy)*y+(cr*sp*sy-sr*cy)*z
        vd=-sp*x+sr*cp*y+cr*cp*z
    else:
        # 局部世界系与 NED 的同一机体航向之差给出坐标系旋转角。
        angle=yaw-np.unwrap(np.deg2rad(f['yaw']))
        vn=np.cos(angle)*v[:,0]-np.sin(angle)*v[:,1]
        ve=np.sin(angle)*v[:,0]+np.cos(angle)*v[:,1]; vd=v[:,2]
    valid=inside & np.isfinite(vn+ve+vd+mapped)
    if 'is_hori_valid' in f: valid &= np.asarray(f['is_hori_valid'])!=0
    if 'is_vert_valid' in f: valid &= np.asarray(f['is_vert_valid'])!=0
    if 'using_fake_fusion' in f: valid &= np.asarray(f['using_fake_fusion'])==0
    result={int(round(mapped[i])):{'vel_n':vn[i],'vel_e':ve[i],'vel_d':vd[i],
             'speed':float(np.hypot(vn[i],ve[i]))} for i in np.flatnonzero(valid)}
    if not result: raise ValueError('转换后没有有效 fusion 数据')
    return result
