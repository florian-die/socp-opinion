import numpy as np
import matplotlib.pyplot as plt
import csv
import fun_opinion as op

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
xi = yi+x0.reshape(-1,1) # agents absolute opinion 
u = data[:,2*n+1] # control
p0 = data[:,n+1] # leader costate
pi = data[:,n+2:2*n+1] # agents costates

u_th = u.copy()
#u_th2 = u.copy()

# singular control computation
for k in range(t.shape[0]):
    u_th[k] = op.u_sing(yi[k,:],pi[k,:])
#    u_th2[k] = op.u_sing2(yi[k,:],pi[k,:])

# plot opinion
fig, ax = plt.subplots(2,1,figsize=(10,10))

ax[0].plot(t,x0,'r',label='leader')
ax[0].plot(t,xi,'b',label='agents')

ax[0].set(xlabel='time')
ax[0].set(ylabel='opinion')

ax[0].grid()
ax[0].legend(framealpha=1)

ax[1].plot((0,t[-1]),(5,5),'k--')
ax[1].plot((0,t[-1]),(-5,-5),'k--',label='saturation')
ax[1].plot(t,u,'r',label='applied control')
ax[1].plot(t,u_th,'g--',label='singular control')
#ax[1].plot(t,u_th2,'m:',label='singular control')

ax[1].set(xlabel='time')
ax[1].set(ylabel='control')

ax[1].set_ylim(-5.5, 5.5)

ax[1].grid()
ax[1].legend(framealpha=1,loc="upper right")

plt.savefig('opinion_state_control.png')
plt.savefig('opinion_state_control.pdf')

plt.show()


