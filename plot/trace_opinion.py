import numpy as np
import matplotlib.pyplot as plt
import csv
import fun_opinion as op
from matplotlib.gridspec import GridSpec

# .dat -> .csv
datContent = [i.strip().split() for i in open("./opinion_data.dat").readlines()]
with open("./opinion_data.csv","wb") as csv_file:
    writer = csv.writer(csv_file)
    writer.writerows(datContent)

# read csv
data = np.genfromtxt('opinion_data.csv',delimiter=',')

# number of agents
n = (data.shape[1]-3)/2

t = data[:,0] # time
x0 = data[:,1] # leader opinion
yi = data[:,2:n+1] # agents relative opinions
y1 = yi[:,0]
yn = yi[:,-1]
xi = yi+x0.reshape(-1,1) # agents absolute opinion
u = data[:,2*n+1] # control
p0 = data[:,n+1] # leader costate
pi = data[:,n+2:2*n+1] # agents costates
p1 = pi[:,0]
pn = pi[:,-1]


""" parameters """
eta = 0.5
sigma = 5.0


u_th = op.ddh(y1)*op.h(y1) - op.ddh(yn)*op.h(yn)
u_th = u_th / (op.ddh(yn) - op.ddh(y1))

#u_th = u.copy()
#u_th2 = u.copy()

# singular control computation
#for k in range(t.shape[0]):
#    u_th[k] = op.u_sing(yi[k,:],pi[k,:])
#    u_th2[k] = op.u_sing2(yi[k,:],pi[k,:])

# plot opinion
fig = plt.figure(constrained_layout=True,figsize=(20,10))
gs = GridSpec(3, 2, figure=fig)

ax1 = fig.add_subplot(gs[0, 0])

ax1.plot(t,x0,'r',label='leader')
ax1.plot(t,xi[:,1],'b',label='agents')
ax1.plot(t,xi,'b')

ax1.set(xlabel='time')
ax1.set(ylabel='opinion')
ax1.grid()
ax1.legend(framealpha=1)

ax2 = fig.add_subplot(gs[1, 0])

ax2.plot(t,u,'r.',ms=1,label='applied control')
ax2.plot(t,u_th,'g--',label='singular control')
ax2.plot((0,t[-1]),(sigma,sigma),'k--')
ax2.plot((0,t[-1]),(-sigma,-sigma),'k--',label='saturation')
#ax2.plot(t,u_th2,'m:',label='singular control')

ax2.set(xlabel='time')
ax2.set(ylabel='control')
ax2.set_ylim(-sigma-0.5, sigma+0.5)
ax2.grid()
ax2.legend(framealpha=1,loc="upper right")

ax3 = fig.add_subplot(gs[2, 0])

ax3.plot(t,p1,'c',label='p1')
ax3.plot(t,-p1,'c--',label='-p1')
ax3.plot(t,pn,'b',label='pn')
ax3.plot(t,-pn,'b--',label='-pn')
ax3.plot(t,p1+pn,'r',label='p1+pn')

ax3.set(xlabel='time')
ax3.set(ylabel='costate')
ax3.grid()
ax3.legend(framealpha=1)

ax4 = fig.add_subplot(gs[:, 1])


ax4.plot(y1,yn,'c',label='extreme agents trajectory')
ax4.plot(y1[0],yn[0],'c*',label='initial opinions')

ymax = np.max(np.hstack((y1,yn)))
ymin = np.min(np.hstack((y1,yn)))
ax4.plot((ymax,-ymax),(ymax,-ymax),'k--')
ax4.plot((ymax,-ymax),(-ymax,ymax),'k--')
ax4.plot((0,0),(-ymax,ymax),'k-',lw=1)
ax4.plot((ymax,-ymax),(0,0),'k-',lw=1)
ax4.plot((eta,-eta,-eta),(eta,eta,-eta),'m-',label='final constraint')

a = np.linspace(-14,-1.1,100)
b = np.sqrt(1.0+4.0/(a**2-1.0))
ax4.plot(-b,-a,'k--')
#ax4.plot(b,-a,'k--')
ax4.plot(a,b,'k--')
#ax4.plot(a,-b,'k--')

ax4.set(xlabel='y1')
ax4.set(ylabel='yn')
ax4.grid()
ax4.legend(framealpha=1)

plt.savefig('opinion_state_control.png')
plt.savefig('opinion_state_control.pdf')

plt.show()
