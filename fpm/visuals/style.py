"""Shared scientific identities, published style profiles and export helpers."""
from pathlib import Path
import hashlib,json
import numpy as np,pandas as pd,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

ORDER=['truncated_square_488','honeycomb','kagome','square','triangular','rr_k3','rr_k4','rr_k6']
STYLE={'truncated_square_488':('4.8.8','#6F3FA0','8'),'honeycomb':('Honeycomb','#E69F00','h'),'kagome':('Kagome','#CC79A7','D'),'square':('Square','#0072B2','s'),'triangular':('Triangular','#D55E00','^'),'rr_k3':('RRG3','#9E9E9E','o'),'rr_k4':('RRG4','#666666','o'),'rr_k6':('RRG6','#222222','o')}
STATE_COLORS={'infected':'#4DAF4A','protected':'#377EB8'}
LANDMARKS={'establishment':('Establishment midpoint','#222222','o'),'reach_valley':('Reach minimum','#0072B2','*'),'q_time_max':('Time-averaged protection maximum','#6F3FA0','p'),'q_final_max':('Terminal protection maximum','#009E73','D'),'rebound_ci_separated':('Rebound (cell-CI rule)','#D55E00','^'),'net_flux_max':('Net strategy-balance maximum','#4D4D4D','v')}

# Profile differences reproduce the approved printed compositions. Network and
# state identities above are shared by every profile.
PROFILES={'standard': {'font.family': 'STIXGeneral', 'mathtext.fontset': 'stix', 'font.size': 8, 'axes.labelsize': 9, 'axes.linewidth': 0.65, 'lines.linewidth': 0.85, 'pdf.fonttype': 42, 'xtick.labelsize': 7, 'ytick.labelsize': 8, 'axes.spines.top': False, 'axes.spines.right': False}, 'establishment': {'font.family': 'STIXGeneral', 'mathtext.fontset': 'stix', 'font.size': 9, 'axes.labelsize': 10, 'xtick.labelsize': 8.2, 'ytick.labelsize': 8.2, 'axes.linewidth': 0.7, 'axes.spines.top': False, 'axes.spines.right': False, 'xtick.major.width': 0.65, 'ytick.major.width': 0.65, 'xtick.major.size': 3.1, 'ytick.major.size': 3.1, 'xtick.direction': 'out', 'ytick.direction': 'out', 'pdf.fonttype': 42, 'ps.fonttype': 42, 'svg.fonttype': 'none', 'svg.hashsalt': 'FPM_FIG3_ALL_CI_20261022', 'savefig.facecolor': 'white', 'axes.unicode_minus': False}, 'crossovers': {'font.family': 'STIXGeneral', 'mathtext.fontset': 'stix', 'font.size': 8, 'axes.labelsize': 9, 'axes.linewidth': 0.65, 'lines.linewidth': 0.85, 'pdf.fonttype': 42, 'xtick.labelsize': 7, 'ytick.labelsize': 8, 'axes.spines.top': False, 'axes.spines.right': False}, 'landmarks': {'font.family': 'STIXGeneral', 'mathtext.fontset': 'stix', 'font.size': 8, 'axes.labelsize': 9, 'axes.linewidth': 0.65, 'pdf.fonttype': 42, 'xtick.labelsize': 7, 'ytick.labelsize': 8, 'axes.spines.top': False, 'axes.spines.right': False}, 'supplement': {'font.family': 'STIXGeneral', 'mathtext.fontset': 'stix', 'font.size': 8, 'axes.labelsize': 9, 'axes.linewidth': 0.65, 'lines.linewidth': 0.8, 'pdf.fonttype': 42, 'xtick.labelsize': 7, 'ytick.labelsize': 7, 'axes.spines.top': False, 'axes.spines.right': False}, 'paired': {'font.family': 'STIXGeneral', 'mathtext.fontset': 'stix', 'font.size': 8, 'axes.labelsize': 9, 'axes.linewidth': 0.65, 'pdf.fonttype': 42, 'xtick.labelsize': 8, 'ytick.labelsize': 8, 'axes.spines.top': False, 'axes.spines.right': False}}

class PlotContext:
    def __init__(self,root,profile,namespace='main'):
        self.root=Path(root);self.profile=profile;self.namespace=namespace
        self.order=list(ORDER);self.style=STYLE;self.qa=[];self.coordinates={};self.exports=[]
        if profile=='establishment':self.order=['honeycomb','truncated_square_488']+ORDER[2:]
        self.output=self.root/('manuscript/figures' if namespace=='main' else 'supplement/figures')
        if profile=='establishment' or profile=='crossovers':self.output=self.root/'figures'
        self.preview=self.root/'qa/figure_previews'
        self.output.mkdir(parents=True,exist_ok=True);self.preview.mkdir(parents=True,exist_ok=True)
        (self.root/'qa').mkdir(exist_ok=True)
    def __enter__(self):
        self.previous_rc=matplotlib.rcParams.copy()
        plt.rcdefaults();plt.rcParams.update(PROFILES[self.profile])
        return self
    def __exit__(self,*args):
        plt.close('all');matplotlib.rcParams.update(self.previous_rc)
    def export(self,fig,path,**kwargs):
        path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
        fig.savefig(path,**kwargs)
        self.exports.append(path)
    @staticmethod
    def clean(ax):ax.grid(alpha=.13,lw=.4);ax.set_axisbelow(True)
    @staticmethod
    def panel(ax,letter):ax.text(0,1.04,letter,transform=ax.transAxes,ha='left',va='bottom',fontweight='bold',fontsize=10)
    def handles(self,keys=None):
        return [Line2D([],[],color=STYLE[t][1],marker=STYLE[t][2],ms=4,label=STYLE[t][0]) for t in (self.order if keys is None else keys)]
    @staticmethod
    def legend(fig,items,ncol,**kw):
        seq=[items[i] for col in range(ncol) for i in range(col,len(items),ncol)]
        return fig.legend(handles=seq,ncol=ncol,frameon=False,fontsize=7,**kw)
    def topology_legend(self,fig,extra=None):
        self.legend(fig,self.handles(self.order[:5]),5,loc='lower center',bbox_to_anchor=(.53,.105))
        self.legend(fig,self.handles(self.order[5:]),3,loc='lower center',bbox_to_anchor=(.53,.055))
        if extra:self.legend(fig,extra,len(extra),loc='lower center',bbox_to_anchor=(.53,.005))
    @staticmethod
    def text_bounds(fig,tolerance):
        fig.canvas.draw();renderer=fig.canvas.get_renderer();outside=[]
        for text in fig.findobj(matplotlib.text.Text):
            if not text.get_visible() or not text.get_text():continue
            box=text.get_window_extent(renderer)
            if box.x0 < -tolerance or box.y0 < -tolerance or box.x1>fig.bbox.width+tolerance or box.y1>fig.bbox.height+tolerance:outside.append(text.get_text())
        return outside
    def save(self,fig,name,records,folder='supplement/figures'):
        if self.namespace=='main':
            for ax in fig.axes:
                xl,yl=ax.get_xlim(),ax.get_ylim()
                for axis,lim,setter in [(ax.xaxis,xl,ax.set_xticks),(ax.yaxis,yl,ax.set_yticks)]:
                    ticks=axis.get_majorticklocs();inside=[v for v in ticks if min(lim)<=v<=max(lim)]
                    if len(inside)!=len(ticks):setter(inside)
                ax.set_xlim(xl);ax.set_ylim(yl)
            outside=self.text_bounds(fig,1);assert not outside,(name,outside)
            self.export(fig,self.output/(name+'.pdf'),metadata={'CreationDate':None,'ModDate':None,'Creator':'FPM audited candidate renderer'})
            self.export(fig,self.preview/(name+'.png'),dpi=170)
            self.qa.append(dict(file=name+'.pdf',width_in=fig.get_figwidth(),height_in=fig.get_figheight(),plotted_records=len(records),text_outside_page=outside))
            self.coordinates[name]=records
        else:
            outside=self.text_bounds(fig,.5);assert not outside,(name,outside)
            pdf=self.root/folder/(name+'.pdf')
            self.export(fig,pdf,dpi=300,metadata={'CreationDate':None,'ModDate':None,'Creator':'FPM SI native-data visual standardization'})
            self.export(fig,self.root/'qa/previews'/(name+'.png'),dpi=180)
            (self.root/'data/rendered').mkdir(parents=True,exist_ok=True)
            pd.DataFrame(records).to_csv(self.root/'data/rendered'/(name+'.csv'),index=False,float_format='%.17g')
            qa=dict(status='PASS',font='STIXGeneral',math_font='stix',records=len(records),text_outside_page=outside,pdf_sha256=hashlib.sha256(pdf.read_bytes()).hexdigest(),raster_backgrounds=0)
            (self.root/'qa'/('RENDER_'+name+'.json')).write_text(json.dumps(qa,indent=2)+'\n',encoding='utf-8')
        plt.close(fig)
