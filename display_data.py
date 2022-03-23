import numpy as np
import matplotlib.pyplot as plt
import csv
import sys
from mpl_toolkits.mplot3d import Axes3D
from os.path import exists


#Get the coordinates for a circle centered and rotated here, with radius R
def get_circle(x,y,z,R,theta):
    phis = np.linspace(0,2*np.pi,64);
    return R*np.cos(phis)*np.sin(theta)+x,R*np.cos(phis)*np.cos(theta)+y,R*np.sin(phis)+z;

def get_circleZ(x,y,z,R):
    phis = np.linspace(0,2*np.pi,64);
    return R*np.cos(phis)+x,R*np.sin(phis)+y,z;


if len(sys.argv) < 2:
        raise ValueError('Usage ',sys.argv[0],' data_folder [optionally --3D (default) --x x_val --y y_val --z z_val]')


data_folder = sys.argv[1];

plotB=False;
plotE=False;
plot_min=0;
plot_max=0;
plot_lim=True
plot_steps=1;
particle_number = 0;
timesteps = [];
max_time = 0;
dt = 0;
names=[]
particle_colors = [':b','--k','-.r','-g','-m','-y','-c']#The default color map includes colors which are hard to read or tell apart
# open the CSV data setup file for setting up the
with open(data_folder+'/setup.csv', mode ='r') as setup:
    # reading the CSV file
    csvFile = csv.DictReader(setup,['name','val'])

    for item in csvFile :
        if (item['name'].strip()=='display_E_field'):
            plotE=bool(item['val'].strip()=='True')
        elif (item['name'].strip()=='display_B_field'):
            plotB=bool(item['val'].strip()=='True')
        elif (item['name'].strip()=='field_min'):
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

#We might want to draw a circle at position x,y,z rotated theta around the z axis and with radius R


circles = 0
c_x = []
c_y = []
c_z = []
c_theta = []
c_R = []
lengthunit = '';

circlesZ = 0
cz_x = []
cz_y = []
cz_z = []
cz_R = []


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
                elif (words[0].strip()=='circle'):
                    c_x.append(float(words[1]))
                    c_y .append(float(words[2]))
                    c_z .append(float(words[3]))
                    c_R .append(float(words[4]))
                    c_theta.append(float(words[5]))
                    circles=circles+1;
                elif (words[0].strip()=='circleZ'):
                    cz_x.append(float(words[1]))
                    cz_y .append(float(words[2]))
                    cz_z .append(float(words[3]))
                    cz_R .append(float(words[4]))
                    circlesZ=circlesZ+1;
                elif (words[0].strip()=='lengthunit'):
                    lengthunit = '/'+words[1];
                else:
                    print(item)
                    print('Option not valid')

if (plotB):
    Bdata = np.fromfile(data_folder+"/B.bin",  dtype=np.float64)
    Bxdata = Bdata[0:plot_steps**3].reshape(plot_steps,plot_steps,plot_steps);
    Bydata = Bdata[plot_steps**3:2*plot_steps**3].reshape(plot_steps,plot_steps,plot_steps);
    Bzdata = Bdata[2*plot_steps**3:3*plot_steps**3].reshape(plot_steps,plot_steps,plot_steps);

if (plotE):
    Edata = np.fromfile(data_folder+"/E.bin",  dtype=np.float64)
    Exdata = Edata[0:plot_steps**3].reshape(plot_steps,plot_steps,plot_steps);
    Eydata = Edata[plot_steps**3:2*plot_steps**3].reshape(plot_steps,plot_steps,plot_steps);
    Ezdata = Edata[2*plot_steps**3:3*plot_steps**3].reshape(plot_steps,plot_steps,plot_steps);

if (plotE or plotB):
    Edata = np.fromfile(data_folder+"/pos_ref.bin",  dtype=np.float64)
    xdata= Edata[0:plot_steps**3].reshape(plot_steps,plot_steps,plot_steps);
    ydata= Edata[plot_steps**3:2*plot_steps**3].reshape(plot_steps,plot_steps,plot_steps);
    zdata= Edata[2*plot_steps**3:3*plot_steps**3].reshape(plot_steps,plot_steps,plot_steps);




#See if we should plot 3D or 2D
if (len(sys.argv)==2):
    plot3D=True;
else: #Loop the arguments and see what we should plot
    plot3D=False;
    skip_arg = False;# I will need to ignore the second argument of --x x_val when looping, python style for-loops don't simply allow me to do ++i, so I need a flag
    for i in range(2,len(sys.argv)):
        if skip_arg:
            skip_arg=False;
        else:
            if (sys.argv[i]=='--3D'):
                plot3D=True;
            else:#2D plotting functions
                skip_arg=True
                if (len(sys.argv)==i+1):
                    print('Missing argument to ',sys.argv[i])
                else:
                    plane = float(sys.argv[i+1]);
                    #We need to find the index closest to this plane, such that
                    #plot_min+(plot_max-plot_min)/plot_steps*i=plane
                    #i.e.
                    plane_i = int( plot_steps*(plane-plot_min)/(plot_max-plot_min))

                    #Make a 3D plot of the fields
                    fig = plt.figure()
                    ax = fig.add_subplot(111)
                    ax.axis('equal')

                    if (sys.argv[i]=='--x'):
                        print('Plot x ',plane_i)
                        if (plotB or plotB):
                            x = ydata[plane_i,:,:]
                            y = zdata[plane_i,:,:]


                        if (plotB):
                            u = Bydata[plane_i,:,:]
                            v = Bzdata[plane_i,:,:]

                            N = 1/np.sqrt(u*u+v*v);
                            #plot normalized
                            ax.quiver(x, y, u*N, v*N,  color = 'red',label='B')

                        if (plotE):
                            u = Eydata[plane_i,:,:]
                            v = Ezdata[plane_i,:,:]

                            N = 1/np.sqrt(u*u+v*v);
                            ax.quiver(x, y, u*N, v*N,  color = 'green',label='E')
                        for i in range(0,particle_number):

                            particle_data = np.fromfile(data_folder+"/particles"+str(i)+".bin",  dtype=np.float64).reshape(timesteps[i],3);
                            if (len(names)>0):
                                ax.plot(particle_data[:,1], particle_data[:,2],particle_colors[i],label=r'$'+(names[i])+'$')
                            else:
                                ax.plot(particle_data[:,1], particle_data[:,2],particle_colors[i])
                        if (circles != 0):
                            for  i in range(0,circles):
                                cx,cy,cz = get_circle(c_x[i],c_y[i],c_z[i],c_R[i],c_theta[i]);
                                ax.plot(cy,cz, color = 'grey')
                            for  i in range(0,circlesZ):
                                cx,cy,cz = get_circleZ(cz_x[i],cz_y[i],cz_z[i],cz_R[i]);
                                ax.plot(cy,cz, color = 'grey')


                        ax.set_xlabel('y'+lengthunit)
                        ax.set_ylabel('z'+lengthunit)
                    elif (sys.argv[i]=='--y'):
                        print('Plot y',plane_i)
                        if (plotB or plotB):
                            x = xdata[:,plane_i,:]
                            y = zdata[:,plane_i,:]
                        if (plotB):
                            u = Bxdata[:,plane_i,:]
                            v = Bzdata[:,plane_i,:]

                            N = 1/np.sqrt(u*u+v*v);
                            ax.quiver(x, y, u*N, v*N,  color = 'red',label='B')

                        if (plotE):
                            u = Exdata[:,plane_i,:]
                            v = Ezdata[:,plane_i,:]

                            N = 1/np.sqrt(u*u+v*v);
                            ax.quiver(x, y, u*N, v*N,  color = 'green',label='E')
                        for i in range(0,particle_number):
                            particle_data = np.fromfile(data_folder+"/particles"+str(i)+".bin",  dtype=np.float64).reshape(timesteps[i],3);
                            if (len(names)>0):
                                ax.plot(particle_data[:,0], particle_data[:,2],particle_colors[i], label=r'$'+(names[i])+'$')
                            else:
                                ax.plot(particle_data[:,0], particle_data[:,2],particle_colors[i])
                        if (circles != 0):
                            for  i in range(0,circles):
                                cx,cy,cz = get_circle(c_x[i],c_y[i],c_z[i],c_R[i],c_theta[i]);
                                ax.plot(cx,cz, color = 'grey')
                            for  i in range(0,circlesZ):
                                cx,cy,cz = get_circleZ(cz_x[i],cz_y[i],cz_z[i],cz_R[i]);
                                ax.plot(cx,cz, color = 'grey')
                        ax.set_xlabel('x')
                        ax.set_ylabel('z')
                    elif (sys.argv[i]=='--z'):
                        print('Plot z',plane_i)
                        if (plotB or plotB):
                            x = xdata[:,:,plane_i]
                            y = ydata[:,:,plane_i]
                        if (plotB):
                            u = Bxdata[:,:,plane_i]
                            v = Bydata[:,:,plane_i]

                            N = 1/np.sqrt(u*u+v*v);
                            ax.quiver(x, y, u*N, v*N, color = 'red',label='B')

                        if (plotE):
                            u = Exdata[:,:,plane_i]
                            v = Eydata[:,:,plane_i]

                            N = 1/np.sqrt(u*u+v*v);
                            ax.quiver(x, y, u*N, v*N, color = 'green',label='E')
                        for i in range(0,particle_number):
                            particle_data = np.fromfile(data_folder+"/particles"+str(i)+".bin",  dtype=np.float64).reshape(timesteps[i],3);
                            if (len(names)>0):
                                ax.plot(particle_data[:,0], particle_data[:,1],particle_colors[i], label=r'$'+(names[i])+'$')
                            else:
                                ax.plot(particle_data[:,0], particle_data[:,1],particle_colors[i])
                        if (circles != 0):
                            for  i in range(0,circles):
                                cx,cy,cz = get_circle(c_x[i],c_y[i],c_z[i],c_R[i],c_theta[i]);
                                ax.plot(cx,cy, color = 'grey')
                            for  i in range(0,circlesZ):
                                cx,cy,cz = get_circleZ(cz_x[i],cz_y[i],cz_z[i],cz_R[i]);
                                ax.plot(cx,cy, color = 'grey')
                        ax.set_xlabel('x')
                        ax.set_ylabel('y')
                    else:
                        print('Argument not recognized')

                    ax.grid()






                    if (plot_min!=plot_max and plot_lim):#If a set range was given, use it
                        ax.set_xlim([plot_min, plot_max])
                        ax.set_ylim([plot_min, plot_max ])
                    plt.legend();
                    plt.show()


if (plot3D):
    #Make a 3D plot of the fields
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    fig.subplots_adjust(left=-0.2, right=1.2, bottom=-0.0, top=1.2)
    ax.view_init(10, -80)
    if (plotB):
        u = Bxdata
        v = Bydata
        w = Bzdata
        N = 1/np.sqrt(u*u+v*v+w*w);
        ax.quiver(xdata, ydata, zdata, N*u, N*v, N*w, length=(plot_max-plot_min)/plot_steps, color = 'red',label='B', alpha=0.5)

    if (plotE):
        u = Exdata
        v = Eydata
        w = Ezdata
        N = 1/np.sqrt(u*u+v*v+w*w);
        ax.quiver(xdata, ydata, zdata, N*u, N*v, N*w, length=(plot_max-plot_min)/plot_steps, color = 'green',label='E', alpha=0.5)


    for i in range(0,particle_number):
        particle_data = np.fromfile(data_folder+"/particles"+str(i)+".bin",  dtype=np.float64).reshape(timesteps[i],3);
        print(i);
        if (len(names)>0):
            print(i)
            print(names[i])
            ax.plot3D(particle_data[ : ,0 ], particle_data[ : ,1 ],particle_data[ : ,2 ],particle_colors[i],  label=r'$'+(names[i])+'$')
        else:
            ax.plot3D(particle_data[ : ,0 ], particle_data[ : ,1 ],particle_data[ : ,2 ],particle_colors[i])

    if (circles != 0):
        for  i in range(0,circles):
            cx,cy,cz = get_circle(c_x[i],c_y[i],c_z[i],c_R[i],c_theta[i]);
            ax.plot3D(cx,cy,cz, color = 'grey')
        for  i in range(0,circlesZ):
            cx,cy,cz = get_circleZ(cz_x[i],cz_y[i],cz_z[i],cz_R[i]);
            ax.plot3D(cx,cy,cz, color = 'grey')


    if (plot_min!=plot_max and plot_lim):#If a set range was given, use it
        ax.set_xlim([plot_min, plot_max])
        ax.set_ylim([plot_min, plot_max ])
        ax.set_zlim([plot_min, plot_max ])
    ax.set_xlabel('x'+lengthunit)
    ax.set_ylabel('y'+lengthunit)
    ax.set_zlabel('z'+lengthunit)
    print(names)
    plt.legend(loc='lower right');
    plt.show()
