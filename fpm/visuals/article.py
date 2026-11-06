"""Current article views composed from shared measured-data plot recipes.

This module changes presentation only. Dynamics, observations, confidence
endpoints and measured coordinates are inherited without new simulations.
"""
from pathlib import Path
import hashlib,json,shutil
import numpy as np,pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.container import ErrorbarContainer
from matplotlib.colors import to_hex
from . import style
from .style import PlotContext, LANDMARKS, ORDER, STYLE
from .response import parameter_response
from .topology import response_landmarks, degree_matched_reach

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def dump(p,x):p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(x,indent=2)+'\n',encoding='utf-8')

def trim_drawn_coordinates(ax,maximum):
    """Limit presentation objects to measured coordinates; retain source tables."""
    for line in ax.lines:
        x=np.asarray(line.get_xdata(),dtype=float);y=np.asarray(line.get_ydata(),dtype=float)
        if len(x)==len(y):line.set_data(x[x<=maximum],y[x<=maximum])
    for collection in ax.collections:
        if hasattr(collection,'get_segments'):
            collection.set_segments([s for s in collection.get_segments() if len(s) and float(s[0][0])<=maximum])

def _compose(outer,view):
    R=outer.root
    work=R/'tmp/render_context'
    for d in ['data','figures','qa']:(work/d).mkdir(parents=True,exist_ok=True)
    previews=R/'qa/figure_previews';previews.mkdir(parents=True,exist_ok=True)
    coords=R/'qa/coordinates';coords.mkdir(exist_ok=True)
    for name in ['fig4_values.csv','fig4_beta_dagger.txt','fig8_landmarks.csv','fig11_contrasts.csv']:
        shutil.copy2(R/'inputs'/name,work/'data'/name)
    V=pd.read_csv(R/'inputs/fig4_values.csv');V=V[V.within_historical_domain]
    bd=float((R/'inputs/fig4_beta_dagger.txt').read_text())
    qa={}
    metadata={'CreationDate':None,'ModDate':None,'Creator':'FPM candidate composition from shared renderer'}
    if view=='peak_duration':
        views=[('Fig4',['peak','duration'],[1,2],R/'manuscript/figures/FIG4_COMPACT_PEAK_DURATION.pdf')]
        for label,observables,kept,pdf in views:
            with PlotContext(work,'standard') as ctx:
                def save_selected(fig,name,records,folder=None):
                    axes=list(fig.axes)
                    for j,ax in enumerate(axes):
                        if j not in kept:ax.remove()
                    fig.set_size_inches(7.1,3.6)
                    for ax,pos,letter in [(axes[kept[0]],[.10,.29,.37,.62],'A'),(axes[kept[1]],[.61,.29,.37,.62],'B')]:
                        ax.set_position(pos);ax.set_title(letter,loc='left',fontweight='bold')
                    outside=ctx.text_bounds(fig,.8);assert not outside,outside
                    selected=[r for r in records if r['observable'] in observables]
                    frame=pd.DataFrame(selected)
                    for ob in observables:
                        for system in ctx.order+['well_mixed','ode']:
                            actual=frame[(frame.observable==ob)&(frame.system==system)].sort_values('beta_n')
                            expected=V[(V.observable==ob)&(V.system==system)].sort_values('beta_n')
                            if system=='ode' and ob!='peak':expected=expected[(expected.beta_n<=bd-.05)|(expected.beta_n>=bd+.05)]
                            assert len(actual)==len(expected),(label,ob,system)
                            columns=['beta_n','mean']+([] if system=='ode' else ['ci_low','ci_high'])
                            np.testing.assert_allclose(actual[columns].to_numpy(float),expected[columns].to_numpy(float),rtol=0,atol=0)
                    # Verify the actual drawn error-bar segments, not just the ledger.
                    segments_checked=0
                    for ax,ob in zip([axes[j] for j in kept],observables):
                        containers=[c for c in ax.containers if isinstance(c,ErrorbarContainer)]
                        assert len(containers)==9
                        for container,system in zip(containers,ctx.order+['well_mixed']):
                            expected=V[(V.observable==ob)&(V.system==system)].sort_values('beta_n')
                            segments=container.lines[2][0].get_segments()
                            assert len(segments)==len(expected)
                            for segment,row in zip(segments,expected.itertuples()):
                                np.testing.assert_allclose(segment,[[row.beta_n,row.ci_low],[row.beta_n,row.ci_high]],rtol=0,atol=1e-12)
                                segments_checked+=1
                    if label=='Fig4':
                        trim_drawn_coordinates(axes[2],15)
                        axes[2].set_xlim(.8,15.25);axes[2].set_xticks([1,5,10,15])
                        for container,system in zip([c for c in axes[2].containers if isinstance(c,ErrorbarContainer)],ctx.order+['well_mixed']):
                            expected=V[(V.observable=='duration')&(V.system==system)&(V.beta_n<=15)].sort_values('beta_n')
                            segments=container.lines[2][0].get_segments()
                            assert len(segments)==len(expected)
                            for segment,row in zip(segments,expected.itertuples()):
                                np.testing.assert_allclose(segment,[[row.beta_n,row.ci_low],[row.beta_n,row.ci_high]],rtol=0,atol=1e-12)
                    ctx.export(fig,pdf,metadata=metadata)
                    ctx.export(fig,previews/(pdf.stem+'.png'),dpi=185)
                    frame['within_visible_x_domain']=frame.beta_n.between(.8,20)
                    if label=='Fig4':frame.loc[frame.observable=='duration','within_visible_x_domain']=frame.loc[frame.observable=='duration','beta_n']<=15
                    frame.to_csv(coords/(label+'_PLOTTED_COORDINATES.csv'),index=False,float_format='%.17g')
                    qa[label]={'status':'PASS','source_sha256':sha(R/'inputs/fig4_values.csv'),
                        'observables':observables,'source_records':len(V),'selected_records':len(frame),
                        'visible_records':int(frame.within_visible_x_domain.sum()),
                        'systems':ctx.order+['well_mixed','ode'],'actual_errorbar_segments_checked':segments_checked,
                        'coordinate_and_interval_equality':True,
                        'displayed_measured_x_max':{'peak':20,'duration':15} if label=='Fig4' else {'reach':20,'q_final':20},
                        'x_limits':{'A':[.8,20.25],'B':[.8,15.25]} if label=='Fig4' else {'A':[.8,20.25],'B':[.8,20.25]},
                        'axis_padding_is_not_extra_measurements':True,
                        'ODE_branch_separator':bd,'ODE_branches_joined':False,
                        'text_outside_canvas':outside,'new_intervals':0,'pdf':pdf.relative_to(R).as_posix()}
                    plt.close(fig)
                ctx.save=save_selected
                parameter_response(ctx)
        qa['Fig4']['author_requested_axis_override']='Panel B measured coordinates through 15; A unchanged through 20'


    elif view=='selected_landmarks':
        keys=['establishment','reach_valley','q_time_max'];actual=[]
        with PlotContext(work,'landmarks') as ctx:
            export=ctx.export
            def export_filtered(fig,path,**kwargs):
                if not getattr(fig,'_filtered_review',False):
                    axes=list(fig.axes);right=axes[1]
                    right.set_xlim(8,20);right.set_xticks([8,10,15,20])
                    for ax in axes:
                        for artist in list(ax.lines):
                            x=np.asarray(artist.get_xdata());y=np.asarray(artist.get_ydata())
                            if len(x)!=1 or len(y)!=1:continue
                            marker=artist.get_marker();matches=[k for k,v in LANDMARKS.items() if v[2]==marker]
                            if len(matches)!=1:continue
                            key=matches[0]
                            if key not in keys:artist.remove();continue
                            new_y=float(y[0])+.15;artist.set_ydata([new_y])
                            actual.append({'topology':ctx.order[round(new_y)],'landmark':key,'beta_n':float(x[0]),'y':new_y,
                                           'marker':marker,'color':LANDMARKS[key][1],'errorbar_drawn':False})
                    right.get_legend().remove()
                    items=[Line2D([],[],color=LANDMARKS[k][1],marker=LANDMARKS[k][2],linestyle='none',
                           markersize=6.1 if LANDMARKS[k][2]=='*' else 4.4,label=LANDMARKS[k][0]) for k in keys]
                    legend=right.legend(handles=items,loc='lower right',frameon=True,fancybox=False,
                           framealpha=1,facecolor='white',edgecolor='#777777',fontsize=7,
                           handletextpad=.5,labelspacing=.4,borderpad=.65,borderaxespad=.7)
                    legend.get_frame().set_linewidth(.6);fig._filtered_review=True
                    outside=ctx.text_bounds(fig,.8);assert not outside,outside
                    fig.canvas.draw();box=legend.get_window_extent(fig.canvas.get_renderer())
                    for row in actual:
                        if row['beta_n']>=8:
                            x,y=right.transData.transform((row['beta_n'],row['y']))
                            assert not box.expanded(1.02,1.15).contains(x,y)
                    qa['Fig8']={'status':'PASS','selected_landmarks':keys,'displayed_records':len(actual),
                        'source_records':48,'x_segments':[[1,2.5],[8,20]],'displayed_errorbars':0,
                        'legend_inside':True,'legend_point_overlap':False,'text_outside_canvas':outside,
                        'source_sha256':sha(R/'inputs/fig8_landmarks.csv'),'new_intervals':0}
                target=(R/'manuscript/figures/FIG8_THREE_LANDMARKS.pdf' if Path(path).suffix=='.pdf'
                        else previews/'FIG8_THREE_LANDMARKS.png')
                export(fig,target,**kwargs)
            ctx.export=export_filtered
            response_landmarks(ctx)
        source=pd.read_csv(R/'inputs/fig8_landmarks.csv');expected=source[source.landmark.isin(keys)]
        plotted=pd.DataFrame(actual)
        joined=plotted.merge(expected,on=['topology','landmark'],suffixes=('_drawn','_source'),validate='one_to_one')
        assert len(joined)==len(actual)==len(expected)==24
        np.testing.assert_array_equal(joined.beta_n_drawn,joined.beta_n_source)
        plotted.to_csv(coords/'Fig8_PLOTTED_COORDINATES.csv',index=False,float_format='%.17g')
        source[~source.landmark.isin(keys)].to_csv(coords/'Fig8_ARCHIVED_NOT_IN_MAIN.csv',index=False,float_format='%.17g')
        qa['Fig8']['author_requested_axis_override']='Right linear segment ends at 20; all 24 selected markers retained'

    elif view=='degree_window':
        with PlotContext(work,'standard') as ctx:
            def save_geometry(fig,name,records,folder=None):
                ax=fig.axes[0];trim_drawn_coordinates(ax,15)
                ax.set_xlim(.5,15.3);ax.set_xticks([1,5,10,15])
                kept_records=[r for r in records if float(r['beta_n'])<=15]
                data=pd.DataFrame(kept_records)
                source=pd.read_csv(R/'inputs/fig11_contrasts.csv')
                expected=source[source.beta_n<=15]
                assert len(data)==len(expected)
                for line,topology in zip(ax.lines[:5],ctx.order[:5]):
                    rows=expected[expected.spatial==topology].sort_values('beta_n')
                    np.testing.assert_array_equal(np.asarray(line.get_xdata()),rows.beta_n.to_numpy())
                    np.testing.assert_array_equal(np.asarray(line.get_ydata()),rows.difference.to_numpy())
                assert len(ax.containers)==0 and len(ax.collections)==0
                assert ax.get_ylim()==(-.55,.3)
                outside=ctx.text_bounds(fig,.8);assert not outside,outside
                pdf=R/'manuscript/figures/FIG11_AUDITED_FIVE_PAIRS.pdf'
                ctx.export(fig,pdf,metadata=metadata);ctx.export(fig,previews/(pdf.stem+'.png'),dpi=185)
                data.to_csv(coords/'Fig11_PLOTTED_COORDINATES.csv',index=False,float_format='%.17g')
                source[source.beta_n>15].to_csv(coords/'Fig11_ARCHIVED_NOT_DISPLAYED.csv',index=False,float_format='%.17g')
                qa['Fig11']={'status':'PASS','source_sha256':sha(R/'inputs/fig11_contrasts.csv'),
                    'source_records':len(source),'displayed_records':len(data),'curves':5,
                    'displayed_measured_x_max':15,'x_limits':[.5,15.3],'x_ticks':[1,5,10,15],
                    'y_limits':[-.55,.3],'errorbar_containers':0,'errorbar_collections':0,
                    'coordinate_equality_to_archived_differences':True,'no_new_differences_computed':True,
                    'text_outside_canvas':outside,'axis_padding_is_not_extra_measurements':True}
                plt.close(fig)
            ctx.save=save_geometry;degree_matched_reach(ctx)

    else:raise ValueError('Unknown presentation view')
    filenames={'peak_duration':'FIG4_COMPACT_PEAK_DURATION.pdf',
               'selected_landmarks':'FIG8_THREE_LANDMARKS.pdf',
               'degree_window':'FIG11_AUDITED_FIVE_PAIRS.pdf'}
    outer.exports.append(R/'manuscript/figures'/filenames[view])
    dump(R/'qa'/('CURRENT_'+view+'.json'),{'status':'PASS','figures':qa,
         'new_simulations':0,'new_intervals':0})
    return qa

def peak_duration(ctx):return _compose(ctx,'peak_duration')
def selected_landmarks(ctx):return _compose(ctx,'selected_landmarks')
def degree_window(ctx):return _compose(ctx,'degree_window')

KEYS=['establishment','reach_valley','rebound_ci_separated','q_time_max','q_final_max','net_flux_max']
LABELS=[r'$\beta_{50}$',r'$\beta_{R,\min}$',r'$\beta_{\uparrow R}$',
        r'$\beta_{\bar q,\max}$',r'$\beta_{q_f,\max}$',r'$\beta_{\mathrm{NSB}}$']
MEANINGS=[
 'Establishment midpoint: P_macro=0.5 at the 5% cumulative-reach criterion.',
 'Measured coordinate of the post-establishment minimum mean reach.',
 'First measured coordinate after the minimum whose lower reach bound exceeds the upper bound at the minimum; archived pointwise interval-separation rule.',
 'Measured coordinate of maximum mean time-averaged active protection.',
 'Measured coordinate of maximum mean terminal active protection.',
 'Measured coordinate of maximum cumulative net strategy balance: (n_N_to_Q - n_Q_to_N)/N.',
]


def _relative_render(outer,mode,raw,relative):
    R=outer.root
    plt.rcParams.update(style.PROFILES['standard'])
    plt.rcParams.update({'axes.labelsize':10,'xtick.labelsize':10,'ytick.labelsize':8,
                         'savefig.facecolor':'white','axes.unicode_minus':False})
    fig,ax=plt.subplots(figsize=(7.1,3.65))
    fig.subplots_adjust(left=.10,right=.985,bottom=.165,top=.97)
    ax.set_xlim(-.53,5.53);ax.set_ylim(0,1.30)
    ax.set_xticks(np.arange(6),LABELS)
    ax.set_yticks(np.arange(0,1.21,.2))
    ax.set_ylabel(r'$\beta_j/\beta_{j,\,4.8.8}$')
    ax.grid(axis='y',color='.91',lw=.45);ax.set_axisbelow(True)
    reference=ax.axhline(1,color='.50',lw=.75,ls=':',zorder=1)
    reference.set_gid('normalization_reference')
    offsets=(np.arange(8)-3.5)*.089
    records=[];data_lines=[];bars=[]
    for i,topology in enumerate(ORDER):
        label,color,marker=STYLE[topology]
        x=np.arange(6)+offsets[i]
        y=relative.loc[topology].to_numpy(float)
        if mode=='bars':
            container=ax.bar(x,y,width=.079,bottom=0,color=color,edgecolor=color,linewidth=.35,zorder=2)
            bars.extend(container.patches)
            for patch,height in zip(container.patches,y):
                assert patch.get_y()==0 and patch.get_height()==height
                assert to_hex(patch.get_facecolor()).lower()==color.lower()
        artist=ax.plot(x,y,linestyle='none',marker=marker,ms=4.15 if mode=='bars' else 4.9,
            color=color,mfc=color,mec=color,mew=.45,zorder=3)[0]
        data_lines.append(artist)
        assert artist.get_linestyle()=='None'
        assert to_hex(artist.get_markerfacecolor()).lower()==color.lower()
        np.testing.assert_array_equal(artist.get_ydata(),y)
        for j,key in enumerate(KEYS):
            records.append({'topology':topology,'landmark':key,'tick_label':LABELS[j],
                'beta_n':float(raw.loc[topology,key]),
                'reference_beta_n':float(raw.loc['truncated_square_488',key]),
                'ratio':float(y[j]),'categorical_x':j,'display_x':float(x[j]),
                'display_offset':float(offsets[i]),'color':color,'marker':marker})
    handles=[Line2D([],[],linestyle='none',marker=STYLE[t][2],ms=4.4,mfc=STYLE[t][1],
                   mec=STYLE[t][1],mew=.4,label=STYLE[t][0]) for t in ORDER]
    ix=[0,4,1,5,2,6,3,7]
    legend=ax.legend(handles=[handles[i] for i in ix],ncol=4,loc='upper center',
                     frameon=True,facecolor='white',edgecolor='.65',framealpha=1,
                     fontsize=7.6,handlelength=1,columnspacing=1.25,handletextpad=.4,
                     borderpad=.5,borderaxespad=.55,labelspacing=.45)
    fig.canvas.draw();renderer=fig.canvas.get_renderer()
    outside=style.PlotContext.text_bounds(fig,1)
    assert not outside,outside
    legend_floor=ax.transData.inverted().transform(legend.get_window_extent(renderer).corners())[:,1].min()
    assert legend_floor>float(relative.to_numpy().max())+.015,legend_floor
    assert len(records)==48 and len(data_lines)==8 and len(bars)==(48 if mode=='bars' else 0)
    stem='LANDMARK_RELATIVE_'+('BARS' if mode=='bars' else 'SYMBOLS')+'_CANDIDATE'
    pdf=(R/'supplement/figures/S6_RELATIVE_TRANSMISSION_LANDMARKS_BARS.pdf' if mode=='bars' else R/'output/alternatives'/(stem+'.pdf'));png=R/'previews'/(stem+'.png');pdf.parent.mkdir(parents=True,exist_ok=True)
    fig.savefig(pdf,metadata={'CreationDate':None,'ModDate':None,
        'Creator':'FPM shared-style categorical landmark renderer',
        'Subject':'Candidate visual alternative; unchanged measured coordinates'})
    fig.savefig(png,dpi=200)
    pd.DataFrame(records).to_csv(R/'data/panel_B'/(stem+'.csv'),index=False,float_format='%.17g')
    receipt={'status':'NUMERIC_PASS_VISUAL_REVIEW_PENDING','mode':mode,'pdf':pdf.relative_to(R).as_posix(),
        'pdf_sha256':sha(pdf),'preview':png.relative_to(R).as_posix(),'records':48,
        'size_inches':[7.1,3.65],'y_limits':[0,1.30],'bar_baseline':0 if mode=='bars' else None,
        'bar_count':len(bars),'markers':48,'data_connecting_segments':0,
        'reference_line_y':1,'reference_equals_one_by_construction':True,
        'x_labels':LABELS,'x_offsets_are_display_only':True,'outside_text':outside,
        'legend_floor_data_y':float(legend_floor),'legend_overlap_with_data':False,
        'network_order':ORDER,'style':STYLE,'source_sha256':sha(R/'inputs/fig8_landmarks.csv'),
        'shared_style_sha256':sha(Path(style.__file__)),'new_simulations':0,
        'new_intervals':0,'new_fits':0,'PNG_inputs':0,'manuscript_integrated':mode=='bars'}
    (R/'qa'/(stem+'.json')).write_text(json.dumps(receipt,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    outer.exports.append(pdf)
    plt.close(fig)
    return receipt

def relative_landmark_bars(ctx):
    R=ctx.root
    for name in ['previews','data/panel_B','qa']:(R/name).mkdir(parents=True,exist_ok=True)
    frame=pd.read_csv(R/'inputs/fig8_landmarks.csv')
    assert len(frame)==48 and not frame.duplicated(['topology','landmark']).any()
    raw=frame.pivot(index='topology',columns='landmark',values='beta_n').loc[ORDER,KEYS]
    assert (raw.loc['truncated_square_488']>0).all()
    relative=raw.div(raw.loc['truncated_square_488'],axis=1)
    prototype=pd.read_csv(R/'data/panel_B/RELATIVE_COORDINATES.csv',index_col=0).loc[ORDER,KEYS]
    np.testing.assert_allclose(relative.to_numpy(),prototype.to_numpy(),rtol=0,atol=5e-16)
    assert np.all(relative.loc['truncated_square_488']==1)
    for t,c in [('rr_k3','#9E9E9E'),('rr_k4','#666666'),('rr_k6','#222222')]:assert STYLE[t][1:]==(c,'o')
    return _relative_render(ctx,'bars',raw,relative)
