# -*- coding: utf-8 -*-
"""
Created on Mon Jun 25 22:27:20 2018

@author: lemitri
"""

import matplotlib.pyplot as plt
import seaborn as sb
sb.set_style('ticks')

size_pp = 20

font = {'family': 'times new roman',
        'color':  'black',
        'weight': 'normal',
        'size': size_pp,
        }
        
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset




#%%

#%%

for g in CHP:
    plt.plot(range(len(day)),[sum(loss_hours[g,d,t] for t in time) for d in day])
    #plt.plot(range(len(day)),[sum(Q_set['HUC',g,d,t] for t in time) for d in day])
        
#%%
for g in HP:
    plt.plot(range(len(day)),[sum(loss_hours[g,d,t] for t in time) for d in day])

#%%
    
for g in incinerator:
    plt.plot(range(len(day)),[sum(loss_hours[g,d,t] for t in time) for d in day])

#%%
    
    
color_list = [sb.xkcd_rgb['navy'],sb.xkcd_rgb['peach'],sb.xkcd_rgb['salmon'],sb.xkcd_rgb['sage'],sb.xkcd_rgb['poop green'],sb.xkcd_rgb['french blue'],sb.xkcd_rgb['blood'],sb.xkcd_rgb['french blue'],sb.xkcd_rgb['cherry'],sb.xkcd_rgb['dark orange'],sb.xkcd_rgb['silver'],sb.xkcd_rgb['charcoal']]    


f,ax=plt.subplots(figsize=(15,5))

sb.despine(offset=1)

plt.subplots_adjust(top=0.915)
plt.subplots_adjust(bottom=0.125)
plt.subplots_adjust(left=0.11)
plt.subplots_adjust(right=0.98)
plt.subplots_adjust(wspace=0.27)



plt.savefig("/Users/lmitridati3/Documents/Postdoc projects/Electricity-aware unit commitment/Python/Electricity Aware Unit Commitment copy/Q_selected_CHPs.eps",dpi=600)
plt.savefig("/Users/lmitridati3/Documents/Postdoc projects/Electricity-aware unit commitment/Python/Electricity Aware Unit Commitment copy/Q_selected_CHPs.png",dpi=600)


#%%

f,(ax1,ax2)=plt.subplots(2,1,figsize=(15,10))

sb.despine(offset=1)

plt.subplots_adjust(top=0.915)
plt.subplots_adjust(bottom=0.125)
plt.subplots_adjust(left=0.11)
plt.subplots_adjust(right=0.98)
#plt.subplots_adjust(wspace=0.27)

for g in [4,5]:

    ax1.plot(range(len(day)),[sum(Q_set['HUC',CHP[g],d,t] for t in time) for d in day],c=color_list[g],label='CHP{0} (decoupled)'.format(g+1))
    ax1.plot(range(len(day)),[sum(Q_set['EAHUC',CHP[g],d,t] for t in time) for d in day],c=color_list[g],ls='--',label='CHP{0} (electricity-aware)'.format(g+1))


for g in [5]:
    ax2.plot(range(len(day)),[sum(loss_hours[CHP[g],d,t] for t in time) for d in day],c=color_list[g],label='CHP{0}'.format(g+1))


ax1.set_ylabel('Heat production (MWh)',fontdict=font)
ax1.set_xlabel('(a)',fontdict=font)

ax1.set_xlim(-1,365.5)
ax1.set_ylim(-0.5,8000.5)

ax1.set_yticks((0,1000,2000,3000,4000,5000,6000,7000,8000))
ax1.set_xticks((0,31,31+28,31+28+31,31+28+31+30,31+28+31+30+31,31+28+31+30+31+30,31+28+31+30+31+30+31,31+28+31+30+31+30+31+31,31+28+31+30+31+30+31+31+30,31+28+31+30+31+30+31+31+30+31,31+28+31+30+31+30+31+31+30+31+30))

ax1.set_yticklabels(('0','1000','2000','3000','4000','5000','6000','7000','8000'), fontdict=font)
ax1.set_xticklabels((('Jan. 1','Feb. 1','Mar. 1','Apr. 1','May 1','Jun. 1','Jul. 1','Aug. 1','Sep. 1','Oct. 1','Nov. 1','Dec. 1')),fontdict=font)


ax2.set_ylabel('Financial losses (Euros)',fontdict=font)
ax2.set_xlabel('Days \n (b)',fontdict=font)

ax2.set_xlim(-1,365.5)
ax2.set_ylim(-500000,0.1)

ax2.set_yticks((-500000,-400000,-300000,-200000,-100000,0))
ax2.set_xticks((0,31,31+28,31+28+31,31+28+31+30,31+28+31+30+31,31+28+31+30+31+30,31+28+31+30+31+30+31,31+28+31+30+31+30+31+31,31+28+31+30+31+30+31+31+30,31+28+31+30+31+30+31+31+30+31,31+28+31+30+31+30+31+31+30+31+30))

ax2.set_yticklabels(('-500,000','-400,000','-300,000','-200,000','-100,000','0'), fontdict=font)
ax2.set_xticklabels((('Jan. 1','Feb. 1','Mar. 1','Apr. 1','May 1','Jun. 1','Jul. 1','Aug. 1','Sep. 1','Oct. 1','Nov. 1','Dec. 1')),fontdict=font)

ax1.legend(bbox_to_anchor=(0.99,0.99),
           bbox_transform=plt.gcf().transFigure,ncol=4,fontsize=17)

plt.savefig("/Users/lmitridati3/Documents/Postdoc projects/Electricity-aware unit commitment/Python/Electricity Aware Unit Commitment copy/losses_selected_CHPs.eps",dpi=600)
plt.savefig("/Users/lmitridati3/Documents/Postdoc projects/Electricity-aware unit commitment/Python/Electricity Aware Unit Commitment copy/losses_selected_CHPs.png",dpi=600)





