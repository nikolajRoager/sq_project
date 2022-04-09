import numpy as np
import matplotlib.pyplot as plt
import csv
import sys
from mpl_toolkits.mplot3d import Axes3D
from os.path import exists


if len(sys.argv) < 2:
        raise ValueError('Usage ',sys.argv[0],' data_folder ')

data_folder = sys.argv[1];

plot_min=0;
plot_max=0;
plot_lim=True
plot_steps=1;
particle_number = 0;
timesteps = [];
max_time = 0;
dt = 0;
names=[]
particle_colors = ['-r','-g','-k']#The default color map includes colors which are hard to read or tell apart
# open the CSV data setup file for setting up the
with open(data_folder+'/setup.csv', mode ='r') as setup:
    # reading the CSV file
    csvFile = csv.DictReader(setup,['name','val'])

    for item in csvFile :
        if (item['name'].strip()=='field_min'):
            plot_min=float(item['val'])#Python being weakly typed, somehow makes it worse at casting things correctly, this would course all mammer of chaos later if I allowed python to interpret this as a 64 bit float
        elif (item['name'].strip()=='field_max'):
            plot_max=float(item['val'])
        elif (item['name'].strip()=='field_steps'):
            plot_steps=int(item['val'])
        elif (item['name'].strip()=='particles'):
            particle_number=int(item['val'])
        elif (item['name'].strip()=='dt'):
            dt=float(item['val'])
        elif (item['name'].strip()=='T'):
            max_time=float(item['val'])
        elif (item['name'].strip()=='name'):
            names.append(item['val'])
        elif (item['name'].strip()=='timesteps'):
            timesteps.append(int(item['val']))
        else:
            print('Entry not understood ',item['name'])


# open a file with extra information to be displayed, this file is optional and does not need to exists
if exists(data_folder+'/extra.txt'):
    with open(data_folder+'/extra.txt', mode ='r') as extra:
        Lines = extra.readlines()

        for item in Lines:
            words = item.split(',')
            if (len(words)!=0):
                if (words[0].strip()=='nofield'):
                    plotE=False
                    plotB=False
                elif (words[0].strip()=='lengthunit'):
                    lengthunit = '/'+words[1];
                else:
                    print(item)
                    print('Option not valid')

fig = plt.figure()
ax = fig.add_subplot(111)
#ax.set_aspect(0.5)


X,Y = np.meshgrid(np.arange(-10,10,0.1) ,np.arange(-10,10,0.1) )

factor =  15;
r = np.sqrt(X*X+Y*Y);

ny = Y/r;
nx = X/r;

Bx = (3.0*((ny*nx)))
By = (((3.0*ny*ny-1)))



plt.streamplot(X,Y,Bx,By, density=1.4, linewidth=None, color='#A23BEC')

for i in range(0,particle_number):
    particle_data = np.fromfile(data_folder+"/particles"+str(i)+".bin",  dtype=np.float64).reshape(timesteps[i],3);
    if (len(names)>0 and False):
        ax.plot(particle_data[:,0], particle_data[:,2],particle_colors[i], label=r'$'+(names[i])+'$',lw=2)
        ax.plot(particle_data[-1,0], particle_data[-1,2],particle_colors[i]+'x', label=r'$'+(names[i])+'$',lw=2)
    else:
        ax.plot(particle_data[:,0], particle_data[:,2],particle_colors[i],lw=2)
        ax.plot(particle_data[-1,0], particle_data[-1,2],particle_colors[i]+'x',lw=2)






if (plot_min!=plot_max and plot_lim):#If a set range was given, use it
    ax.set_xlim([plot_min, plot_max])
    ax.set_ylim([plot_min, plot_max ])
plt.legend();
plt.show()

