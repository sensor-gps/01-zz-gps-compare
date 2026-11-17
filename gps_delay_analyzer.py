"""GPS iTOW/PVT 对齐与延迟分析 GUI。

依赖: Python 3, matplotlib；读取 ULog 真值还需要 pyulog (pip install pyulog)。
运行: python gps_delay_analyzer.py
"""
from __future__ import annotations
import os, re, struct, threading, csv, json
import tkinter as tk
from tkinter import ttk, filedialog, messagebox
import numpy as np

try:
    import matplotlib.pyplot as plt
except ImportError:
    plt = None

FIELDS = ("lat", "lon", "alt", "vel_n", "vel_e", "vel_d", "speed", "heading", "hacc", "vacc", "num_sv")

def _add(d, tow, **kw):
    if tow is not None and tow >= 0:
        d.setdefault(int(tow), {}).update(kw)

def parse_ubx(path):
    data = open(path, "rb").read(); out = {}
    i = 0
    while i + 8 <= len(data):
        j = data.find(b"\xb5\x62", i)
        if j < 0 or j + 8 > len(data): break
        cls, mid, n = data[j+2], data[j+3], struct.unpack_from("<H", data, j+4)[0]
        end = j + 8 + n
        if end > len(data):
            i=j+2
            continue
        p = data[j+6:j+6+n]
        ca=cb=0
        for byte in data[j+2:j+6+n]:
            ca=(ca+byte)&255; cb=(cb+ca)&255
        if data[j+6+n:j+8+n] != bytes((ca,cb)):
            i=j+2
            continue
        if cls == 1 and mid == 7 and n >= 92: # NAV-PVT
            tow = struct.unpack_from("<I", p, 0)[0]
            lon, lat = struct.unpack_from("<ii", p, 24)
            hmsl = struct.unpack_from("<i", p, 36)[0]
            hacc, vacc = struct.unpack_from("<II", p, 40)
            vel_n, vel_e, vel_d, speed, head = struct.unpack_from("<iiiii", p, 48)
            num_sv = p[23]
            # 排除尚无有效定位的启动数据；保留原始 GPS 周内毫秒。
            if not (0 < tow < 604800000 and p[20]>=3 and p[21]&1
                    and -1800000000<=lon<=1800000000 and -900000000<=lat<=900000000):
                i=end
                continue
            _add(out, tow, lat=lat/1e7, lon=lon/1e7, alt=hmsl/1000, vel_n=vel_n/1000,
                 vel_e=vel_e/1000, vel_d=vel_d/1000, speed=speed/1000, heading=head/1e5,
                 hacc=hacc/1000, vacc=vacc/1000, num_sv=num_sv)
        i = end
    return out

def parse_text(path):
    out = {}; raw = open(path, "rb").read()
    for m in re.finditer(rb"\$(G[PN]RMC|G[PN]GGA),([^*\r\n]*)", raw):
        typ, s = m.group(1).decode(), m.group(2).decode(errors="ignore"); a = s.split(',')
        try:
            if typ.endswith("RMC") and len(a) >= 9 and a[1] and a[2] == 'A':
                tow = int(float(a[1])*1000); lat = float(a[2][:2])+float(a[2][2:])/60; lon = float(a[4][:3])+float(a[4][3:])/60
                if a[3]=='S': lat=-lat
                if a[5]=='W': lon=-lon
                _add(out,tow,lat=lat,lon=lon,speed=float(a[6])*0.514444,heading=float(a[7] or 0))
            elif typ.endswith("GGA") and len(a) >= 9 and a[0]:
                tow=int(float(a[0])*1000); lat=float(a[1][:2])+float(a[1][2:])/60; lon=float(a[3][:3])+float(a[3][3:])/60
                if a[2]=='S': lat=-lat
                if a[4]=='W': lon=-lon
                _add(out,tow,lat=lat,lon=lon,alt=float(a[8] or 0),num_sv=int(a[6] or 0))
        except (ValueError, IndexError): pass
    return out

def parse_gps(path):
    # .log 也可能是混合 UBX 二进制日志，不能仅按扩展名判断。
    records=parse_ubx(path)
    if not records:
        raise RuntimeError(f'未提取到校验通过且定位有效的 UBX NAV-PVT：{path}')
    export_dir=os.path.join(os.path.dirname(os.path.abspath(__file__)), 'analysis_output')
    os.makedirs(export_dir,exist_ok=True)
    export_path=os.path.join(export_dir,os.path.basename(path)+'.pvt.csv')
    with open(export_path,'w',newline='',encoding='utf-8-sig') as f:
        writer=csv.DictWriter(f,fieldnames=('itow_ms',)+FIELDS)
        writer.writeheader()
        for t,row in sorted(records.items()): writer.writerow(dict(itow_ms=t,**row))
    return records

def parse_ulog(path):
    try:
        from pyulog import ULog
    except ImportError as e: raise RuntimeError('读取 ULog 需要安装 pyulog：pip install pyulog') from e
    u=ULog(path); topics={x.name:x for x in u.data_list}; out={}
    # drone_fusion 是飞控融合输出：position/velocity 为机体系，yaw 为相对真北航向。
    # 先将机体系 x/y 旋转到 N/E，z 保持 Down 方向。
    for name in ('drone_fusion','vehicle_gps_position','gps_position','vehicle_global_position'):
        if name not in topics: continue
        d=topics[name].data; keys=d.keys(); tk=next((k for k in ('time','timestamp','timestamp_sample') if k in keys),None)
        if not tk: continue
        for i,t in enumerate(d[tk]):
            # 若主题带 iTOW/itow，直接使用；否则用时间戳毫秒作为统一时间轴。
            if 'iTOW' in keys: tow = d['iTOW'][i]
            elif 'itow' in keys: tow = d['itow'][i]
            else: tow = d.get('time_utc_usec', d.get('timestamp', d[tk]))[i] / 1000.0
            kw={}
            if name == 'drone_fusion':
                def get(*ks):
                    k0=next((x for x in ks if x in keys),None)
                    return float(d[k0][i]) if k0 else None
                px,py,pz=get('position_x','pos_x'),get('position_y','pos_y'),get('position_z','pos_z')
                vx,vy,vz=get('velocity_x','vel_x'),get('velocity_y','vel_y'),get('velocity_z','vel_z')
                yaw=get('yaw','heading')
                if yaw is not None:
                    # 常见日志 yaw 为度；弧度值通常绝对值不超过 2*pi。
                    import math
                    yr = yaw if abs(yaw) <= 2*math.pi+0.1 else math.radians(yaw)
                    c,s=math.cos(yr),math.sin(yr)
                    if px is not None: kw['pos_n'],kw['pos_e'],kw['pos_d'] = c*px-s*py, s*px+c*py, pz if pz is not None else 0
                    if vx is not None: kw['vel_n'],kw['vel_e'],kw['vel_d'] = c*vx-s*vy, s*vx+c*vy, vz if vz is not None else 0
                    kw['heading']=math.degrees(yr)%360
                if not kw: continue
                _add(out,tow,**kw)
                continue
            for k, aliases, scale in [('lat',('lat','latitude'),1e-7),('lon',('lon','longitude'),1e-7),('alt',('alt','altitude_msl'),1e-3),('vel_n',('vel_n','vel_n_m_s'),1),('vel_e',('vel_e','vel_e_m_s'),1),('vel_d',('vel_d','vel_d_m_s'),1),('hacc',('eph','hacc'),1e-2),('vacc',('epv','vacc'),1e-2),('num_sv',('satellites_used','num_sv'),1)]:
                k0=next((x for x in aliases if x in keys),None)
                if k0: kw[k]=float(d[k0][i])*scale
            _add(out, tow, **kw)
        if out: return out
    raise RuntimeError('ULog 中未找到 vehicle_gps_position/gps_position 主题')

def align(ref, series, tolerance=200):
    rk=np.array(sorted(ref)); sk=np.array(sorted(series)); rows=[]
    for t in rk:
        j=np.searchsorted(sk,t); cand=[x for x in (j-1,j) if 0<=x<len(sk)]
        if not cand: continue
        k=min(cand,key=lambda q:abs(sk[q]-t))
        if abs(sk[k]-t)<=tolerance:
            row=dict(series[int(sk[k])])
            row['_source_itow_ms']=int(sk[k])
            rows.append((t,row))
    return rows

def show_comparison(parent, reference, ref, selected, aligned):
    from matplotlib.figure import Figure
    from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
    window=tk.Toplevel(parent)
    window.title('GPS PVT 对比 — 分页查看')
    window.geometry('1200x800')
    ttk.Label(window,text='滚轮：以鼠标位置缩放｜左键：放置垂直光标｜←/→：移动光标｜+/-：缩放｜Home：恢复｜工具栏：平移、框选、保存').pack(fill='x',padx=10,pady=6)
    notebook=ttk.Notebook(window); notebook.pack(fill='both',expand=True)
    fields=[('lat','纬度','°'),('lon','经度','°'),('alt','海拔高度','m'),
            ('vel_n','北向速度','m/s'),('vel_e','东向速度','m/s'),('vel_d','地向速度','m/s'),
            ('speed','地速','m/s'),('heading','航向角','°'),('hacc','水平精度','m'),
            ('vacc','垂直精度','m'),('num_sv','卫星数','颗')]
    origin=min(reference)
    for field,title,unit in fields:
        page=ttk.Frame(notebook); notebook.add(page,text=title)
        fig=Figure(figsize=(10,6)); ax=fig.add_subplot(111)
        canvas=FigureCanvasTkAgg(fig,master=page)
        toolbar=NavigationToolbar2Tk(canvas,page,pack_toolbar=False); toolbar.pack(side='bottom',fill='x')
        readout=tk.StringVar(value='点击曲线区域放置垂直对比光标')
        ttk.Label(page,textvariable=readout,anchor='w',justify='left').pack(side='bottom',fill='x',padx=8,pady=6)
        canvas.get_tk_widget().pack(fill='both',expand=True)
        curves=[]
        source_times=[]
        if ref is not None:
            times=np.array(sorted(ref)); curves.append(('FFC 真值',(times-origin)/1000,np.array([ref[t].get(field,np.nan) for t in times])))
            source_times.append(times)
        for (_,path),rows in zip(selected,aligned):
            curves.append((os.path.splitext(os.path.basename(path))[0],
                           np.array([(t-origin)/1000 for t,_ in rows]),np.array([row.get(field,np.nan) for _,row in rows])))
            source_times.append(np.array([row.get('_source_itow_ms',t) for t,row in rows]))
        plot_lines=[]
        for name,x,y in curves:
            line,=ax.plot(x,y,label=name,linewidth=1.1,marker='o',markersize=5,
                          markeredgecolor='white',markeredgewidth=0.6)
            plot_lines.append(line)
        def adapt_markers(event=None,ax=ax,fig=fig,canvas=canvas,curves=curves,plot_lines=plot_lines):
            lo,hi=sorted(ax.get_xlim())
            for line,(_,x,y) in zip(plot_lines,curves):
                left=int(np.searchsorted(x,lo)); right=int(np.searchsorted(x,hi,side='right'))
                count=right-left
                spacing=ax.bbox.width/max(count,1) # 可见采样点的平均像素间距
                size=float(np.clip(spacing*.55*72/fig.dpi,2.2,6))
                line.set_markersize(size)
                line.set_markeredgewidth(.6 if spacing>=7 else 0)
                # 每个真实采样点始终显示，只调整尺寸，避免缩放后标记被跳过。
                line.set_markevery(None)
            canvas.draw_idle()
        ax.callbacks.connect('xlim_changed',adapt_markers)
        canvas.mpl_connect('resize_event',adapt_markers)
        ax.set(title=f'{title}对比',xlabel='相对时间 (s)',ylabel=f'{title} ({unit})'); ax.grid(True,alpha=.35); ax.legend()
        fig.tight_layout(); adapt_markers(); canvas.draw()
        initial=(ax.get_xlim(),ax.get_ylim())
        cursors=[ax.axvline(0,color=color,ls='--',lw=1,visible=False) for color in ('purple','brown')]
        from matplotlib.patches import Rectangle
        box=Rectangle((0,0),0,0,fill=False,edgecolor='gray',linestyle='--',visible=False)
        ax.add_patch(box)
        state={'times':[None,None],'active':0,'annotations':[],'drag':None}
        def place(t,slot=None,ax=ax,canvas=canvas,cursors=cursors,state=state,curves=curves,readout=readout,unit=unit,source_times=source_times,has_ref=ref is not None):
            if slot is None: slot=state['active']
            state['times'][slot]=t
            for annotation in state['annotations']: annotation.remove()
            state['annotations'].clear()
            for k,time in enumerate(state['times']):
                cursors[k].set_visible(time is not None)
                if time is not None: cursors[k].set_xdata([time,time])
            ta,tb=state['times']; lines=[]
            if ta is not None and tb is not None: lines.append(f'A={ta:.6f} s    B={tb:.6f} s    Δt (B−A)={tb-ta:.6f} s')
            else: lines.append('右键依次放置 A、B 光标；左键按住光标线可拖动；按 A/B 选择光标')
            for row,(name,x,y) in enumerate(curves):
                vals=[]
                for k,time in enumerate(state['times']):
                    value=None
                    if time is not None and len(x) and x[0]<=time<=x[-1]:
                        j=int(np.searchsorted(x,time))
                        if j<len(x) and x[j]==time: value=y[j]
                        elif 0<j<len(x) and np.isfinite(y[j-1:j+1]).all(): value=float(np.interp(time,x[j-1:j+1],y[j-1:j+1]))
                    vals.append(value)
                    if value is not None and np.isfinite(value):
                        label=ax.annotate(f'{"AB"[k]} {name}: {value:.8g}',(time,value),xytext=(6 if k==0 else -6,10+row*14),textcoords='offset points',ha='left' if k==0 else 'right',fontsize=8,bbox=dict(boxstyle='round',fc='white',alpha=.8))
                        state['annotations'].append(label)
                va,vb=vals
                fmt=lambda v: '无数据' if v is None else f'{v:.8g} {unit}'
                line=f'{name}: A={fmt(va)}    B={fmt(vb)}'
                if va is not None and vb is not None: line+=f'    ΔY (B−A)={vb-va:.8g} {unit}'
                lines.append(line)
            # 插值交点没有唯一的原始 iTOW；单独列出最近的对齐样本及其真实采样时间。
            for slot,time in enumerate(state['times']):
                if time is None: continue
                samples=[]
                for (name,x,y),original in zip(curves,source_times):
                    if not len(x) or not x[0]<=time<=x[-1]: continue
                    j=int(np.argmin(np.abs(x-time)))
                    samples.append((name,int(original[j])))
                lines.append(f'{"AB"[slot]} 最近采样点原始时间（ms，非插值交点时间）：'+ '； '.join(f'{name}={stamp}' for name,stamp in samples))
                gps_samples=samples[1:] if has_ref and samples and samples[0][0]=='FFC 真值' else samples
                for j,(name,stamp) in enumerate(gps_samples):
                    for other,other_stamp in gps_samples[j+1:]:
                        lines.append(f'{"AB"[slot]} ΔiTOW [{other} − {name}] = {other_stamp-stamp:+d} ms')
            readout.set('\n'.join(lines)); canvas.draw_idle()
        def zoom(factor,cx=None,cy=None,ax=ax,canvas=canvas):
            a,b=ax.get_xlim(); c,d=ax.get_ylim()
            cx=(a+b)/2 if cx is None else cx; cy=(c+d)/2 if cy is None else cy
            ax.set_xlim(cx+(a-cx)*factor,cx+(b-cx)*factor)
            ax.set_ylim(cy+(c-cy)*factor,cy+(d-cy)*factor); canvas.draw_idle()
        def click(event,ax=ax,canvas=canvas,toolbar=toolbar,place=place,state=state,box=box):
            canvas.get_tk_widget().focus_set()
            if event.inaxes!=ax or toolbar.mode: return
            if event.button==3:
                place(event.xdata); state['active']=1-state['active']
            elif event.button in (1,2):
                if event.button==1:
                    # 使用屏幕像素判定命中，缩放后仍保持一致的选择宽度。
                    candidates=[(abs(ax.transData.transform((t,0))[0]-event.x),k)
                                for k,t in enumerate(state['times']) if t is not None]
                    if candidates and min(candidates)[0]<=9:
                        state['active']=min(candidates)[1]
                        state['drag']=('cursor',)
                        place(event.xdata)
                        return
                state['drag']=(event.button,event.x,event.y,event.xdata,event.ydata,ax.get_xlim(),ax.get_ylim())
                if event.button==1:
                    box.set_xy((event.xdata,event.ydata)); box.set_width(0); box.set_height(0); box.set_visible(True)
        def motion(event,ax=ax,canvas=canvas,state=state,box=box,place=place):
            drag=state['drag']
            if drag is None: return
            if drag[0]=='cursor':
                if event.inaxes==ax: place(event.xdata)
                return
            button,px,py,x,y,xlim,ylim=drag
            if button==2:
                dx=(event.x-px)/ax.bbox.width*(xlim[1]-xlim[0]); dy=(event.y-py)/ax.bbox.height*(ylim[1]-ylim[0])
                ax.set_xlim(xlim[0]-dx,xlim[1]-dx); ax.set_ylim(ylim[0]-dy,ylim[1]-dy)
            elif event.inaxes==ax:
                box.set_width(event.xdata-x); box.set_height(event.ydata-y)
            canvas.draw_idle()
        def release(event,ax=ax,canvas=canvas,state=state,box=box):
            drag=state['drag']; state['drag']=None; box.set_visible(False)
            if drag and drag[0]==1 and event.inaxes==ax and abs(event.x-drag[1])>5 and abs(event.y-drag[2])>5:
                ax.set_xlim(sorted((drag[3],event.xdata))); ax.set_ylim(sorted((drag[4],event.ydata)))
            canvas.draw_idle()
        def scroll(event,ax=ax,zoom=zoom):
            if event.inaxes==ax: zoom(.8 if event.button=='up' else 1.25,event.xdata,event.ydata)
        def key(event,ax=ax,canvas=canvas,initial=initial,state=state,place=place,zoom=zoom):
            if event.key in ('+','=','-'): zoom(1.25 if event.key=='-' else .8)
            elif event.key in ('home',' ','space'):
                ax.set_xlim(initial[0]); ax.set_ylim(initial[1]); canvas.draw_idle()
            elif event.key in ('a','b'): state['active']=0 if event.key=='a' else 1
            elif event.key in ('left','right'):
                a,b=ax.get_xlim(); t=state['times'][state['active']]
                place(((a+b)/2 if t is None else t)+(b-a)/100*(1 if event.key=='right' else -1))
        for event,handler in [('button_press_event',click),('motion_notify_event',motion),('button_release_event',release),('scroll_event',scroll),('key_press_event',key)]: canvas.mpl_connect(event,handler)
        page.plot_objects=(fig,canvas,toolbar,click,motion,release,scroll,key,adapt_markers)

class App(tk.Tk):
    def __init__(self):
        super().__init__(); self.title('GPS 延迟分析（iTOW/PVT）'); self.geometry('1000x620'); self.items=[]
        self.config_path=os.path.join(os.path.dirname(os.path.abspath(__file__)),'.gps_delay_history.json')
        left=ttk.Frame(self); left.pack(side='left',fill='y',padx=8,pady=8); right=ttk.Frame(self); right.pack(side='right',fill='both',expand=True,padx=8,pady=8)
        ttk.Label(left,text='分析项目').pack(anchor='w'); self.list=tk.Listbox(left,width=24,height=12); self.list.pack(fill='y',expand=True); self.list.insert(0,'GPS 延迟分析（iTOW/PVT）'); self.list.selection_set(0)
        ttk.Label(right,text='VIO 真值 ffc.csv').grid(row=0,column=0,sticky='w'); self.ref=tk.StringVar(); ttk.Entry(right,textvariable=self.ref,width=75).grid(row=0,column=1); ttk.Button(right,text='选择',command=lambda:self.pick(self.ref,'ulog')).grid(row=0,column=2)
        self.gps=[]
        for i in range(3):
            v=tk.StringVar(); self.gps.append(v); ttk.Label(right,text=f'GPS {i+1} 原始文件').grid(row=i+1,column=0,sticky='w'); ttk.Entry(right,textvariable=v,width=75).grid(row=i+1,column=1); ttk.Button(right,text='选择',command=lambda x=v:self.pick(x)).grid(row=i+1,column=2)
        ttk.Label(right,text='对齐容差(ms)').grid(row=4,column=0,sticky='w'); self.tol=tk.StringVar(value='200'); ttk.Entry(right,textvariable=self.tol,width=10).grid(row=4,column=1,sticky='w')
        self.log=tk.Text(right,height=18); self.log.grid(row=6,column=0,columnspan=3,sticky='nsew'); right.rowconfigure(6,weight=1); right.columnconfigure(1,weight=1); ttk.Button(right,text='开始分析并绘图',command=self.run).grid(row=5,column=1,sticky='e',pady=8)
        self.list.insert('end','RTCM 星历检测')
        scrollbar=ttk.Scrollbar(left,orient='vertical',command=self.list.yview)
        scrollbar.pack(side='right',fill='y'); self.list.configure(yscrollcommand=scrollbar.set)
        self.project_area=right
        self.eph_area=ttk.Frame(self)
        ttk.Label(self.eph_area,text='RTCM 星历检测文件').grid(row=0,column=0,padx=8,pady=12)
        self.eph_path=tk.StringVar()
        ttk.Entry(self.eph_area,textvariable=self.eph_path,width=70).grid(row=0,column=1,sticky='ew')
        ttk.Button(self.eph_area,text='选择',command=lambda:self.pick(self.eph_path)).grid(row=0,column=2,padx=8)
        ttk.Button(self.eph_area,text='解析并绘制星历缺失比例',command=self.run_ephemeris).grid(row=1,column=1,pady=12)
        ttk.Label(self.eph_area,text='支持含 RTCM3 二进制的日志；按 MSM 时间每 10 秒统计。\n左轴：星历接收缺失比例；右轴：该时间窗观测到的卫星数。\n缺失指本文件中尚未收到该卫星的星历，不代表接收机实际无星历。').grid(row=2,column=0,columnspan=3,sticky='w',padx=8)
        self.eph_area.columnconfigure(1,weight=1)
        self.list.bind('<<ListboxSelect>>',self.select_project)
        self.fusion_frame=tk.StringVar(value='local')
        self.north_confirmed=tk.BooleanVar(value=True)
        self.anchor_boot=tk.StringVar(); self.anchor_tow=tk.StringVar()
        self.load_history(); self.protocol('WM_DELETE_WINDOW',self.on_close)
    def select_project(self,event=None):
        selection=self.list.curselection()
        if not selection: return
        self.project_area.pack_forget(); self.eph_area.pack_forget()
        area=self.eph_area if selection[0]==1 else self.project_area
        area.pack(side='right',fill='both',expand=True,padx=8,pady=8)
    def run_ephemeris(self):
        try:
            path=self.eph_path.get().strip()
            if not path: raise ValueError('请选择需要分析的 RTCM 日志文件')
            import rtcm_ephemeris
            if plt is not None:
                plt.rcParams['font.sans-serif']=['Microsoft YaHei','SimHei','DejaVu Sans']
                plt.rcParams['axes.unicode_minus']=False
            rtcm_ephemeris.show(self,path)
        except Exception as e: messagebox.showerror('星历检测失败',str(e))
    def pick(self,v,kind='gps'):
        filetypes=[('所有文件','*.*')] if kind=='ffc' else [('GPS files','*.ubx *.log *.bin *.*'),('所有文件','*.*')]
        f=filedialog.askopenfilename(filetypes=filetypes);
        if f:v.set(f)
    def load_history(self):
        try:
            with open(self.config_path,encoding='utf-8') as f: d=json.load(f)
            saved_ref=d.get('ref','')
            self.ref.set(saved_ref if str(saved_ref).lower().endswith(('.ulg','.ulog')) else '')
            for v,p in zip(self.gps,d.get('gps',[])): v.set(p)
            self.tol.set(str(d.get('tol','200'))); self.fusion_frame.set(d.get('fusion_frame','local'))
            self.north_confirmed.set(bool(d.get('north_confirmed',False)))
            self.anchor_boot.set(str(d.get('anchor_boot',''))); self.anchor_tow.set(str(d.get('anchor_tow','')))
            self.eph_path.set(d.get('eph_path',''))
        except (OSError,ValueError,TypeError,json.JSONDecodeError): pass
    def save_history(self):
        d={'ref':self.ref.get(),'gps':[v.get() for v in self.gps],'tol':self.tol.get(),
           'fusion_frame':self.fusion_frame.get(),'north_confirmed':self.north_confirmed.get(),
           'anchor_boot':self.anchor_boot.get(),'anchor_tow':self.anchor_tow.get(),'eph_path':self.eph_path.get()}
        try:
            with open(self.config_path,'w',encoding='utf-8') as f: json.dump(d,f,ensure_ascii=False,indent=2)
        except OSError: pass
    def on_close(self): self.save_history(); self.destroy()
    def run(self):
        if getattr(self,'_running',False): return
        self.save_history()
        self._running=True
        self.log.insert('end','正在后台解析文件，请稍候，界面仍可响应……\n')
        threading.Thread(target=self._run,daemon=True).start()
    def _run(self):
        try:
            if plt is None: raise RuntimeError('请先安装 matplotlib：pip install matplotlib')
            ref_path=self.ref.get().strip()
            if ref_path and not ref_path.lower().endswith(('.ulg','.ulog')):
                raise ValueError('?? GPS ??????? ULog ?????? .ulg ????? ffc.csv ??????')
            from fusion_reference import load
            ref=None
            selected=[(i+1,v.get().strip()) for i,v in enumerate(self.gps) if v.get().strip()]
            gps=[parse_gps(path) for _,path in selected]
            if not gps: raise RuntimeError('至少选择一个 GPS 文件')
            for (number,path),data in zip(selected,gps):
                if not data: raise RuntimeError(f'GPS {number} 未解析到有效数据：{path}')
            if ref_path:
                if not self.gps[0].get().strip(): raise ValueError('GPS1 is required as the sensor_gps clock anchor')
                ref=load(ref_path,gps1=gps[0])
                self.after(0,lambda:self.log.insert('end','ULog sensor_gps -> GPS1 PVT identity clock mapping; GPS2/3 share iTOW.\n'))
            reference=ref if ref is not None else gps[0]
            if not reference: raise RuntimeError('参考文件中没有可分析的数据')
            tol=float(self.tol.get())
            if not np.isfinite(tol) or tol<0: raise ValueError('对齐容差必须是非负有限数值')
            aligned=[align(reference,x,tol) for x in gps]
            if ref is not None:
                from gps_timing import estimate_delay
                for (number,path),data in zip(selected,gps):
                    delay,score=estimate_delay(ref,data)
                    text=f'GPS {number}: NE lag={delay:.0f} ms, correlation={score:.5f}\n'
                    self.after(0,lambda text=text:self.log.insert('end',text))
            if ref is not None and not any(aligned):
                raise ValueError(f'GPS 与 fusion 的 iTOW 无重叠，请检查时间锚点。fusion 范围 {min(ref)}～{max(ref)} ms')
            mode='ULog drone_fusion NED' if ref is not None else 'GPS only'
            self.after(0,lambda: self.log.insert('end',f'参考：{mode}\n'))
            for (number,_),data,rows in zip(selected,gps,aligned):
                self.after(0,lambda number=number,data=data,rows=rows: self.log.insert('end',f'GPS {number}：{len(data)} 点，匹配 {len(rows)} 点\n'))
            plt.rcParams['font.sans-serif']=['Microsoft YaHei','SimHei','DejaVu Sans']
            plt.rcParams['axes.unicode_minus']=False
            self.after(0,lambda: show_comparison(self,reference,ref,selected,aligned))
        except Exception as e:
            # Python 会在离开 except 块时清理异常变量 e，不能在稍后执行的
            # Tkinter lambda 中直接引用它；先复制成普通字符串。
            error_text = str(e)
            self.after(0, lambda msg=error_text: messagebox.showerror('分析失败', msg))
        finally:
            self._running=False

if __name__=='__main__': App().mainloop()
