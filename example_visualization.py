#%% 
import numpy as np
import matplotlib as plt
from Functions import rk3
from Example import bvector_system_1
###############################################################
##                  The Moderately Stiff Case                ##
###############################################################
#%% 

## 3.1 Input the system_1
a1=1000
a2=1
A=np.array([
    [-a1,  0],
    [ a1,-a2]
])
y0=[1,0]# Rigoriously, y0 should be a col vector, but according to the definition of our y matrix, it is a row vector here
interval=[0,0.1]
bvector=bvector_system_1.bvector_system_1

#%% 

## 3.2 Calculate erroer
k_size=10
Erro_rk3_list=[]
Erro_dirk3_list=[]
h_List=[]

for k in range(1,k_size+1):
    N=40*k

    x_rk3,y_rk3=rk3.rk3(A,bvector,y0,interval,N)
    #y1_rk3=y_rk3[0,:]
    y2_rk3=y_rk3[1,:]

    x_dirk3,y_dirk3=dirk3.dirk3(A,bvector,y0,interval,N)
    #y1_dirk3=y_dirk3[0,:]
    y2_dirk3=y_dirk3[1,:]

    x=x_rk3 #Since x_rk3=x_dirk3
    #y1_exact=np.exp(-a1*x)
    y2_exact=(a1/(a1-a2))*(np.exp(-a2*x)-np.exp(-a1*x))

    Erro_1_norm_rk3_y2=0
    Erro_1_norm_dirk3_y2=0
    h=x[1]-x[0] #stepsize
    #start from j=2 in python start from j=1 since y2==y2_exact
    for j in range(1,len(y2_exact)):
        Erro_1_norm_rk3_y2+=abs((y2_rk3[j]-y2_exact[j])/y2_exact[j])
        Erro_1_norm_dirk3_y2+=abs((y2_dirk3[j]-y2_exact[j])/y2_exact[j])

    Erro_rk3_list.append(h*Erro_1_norm_rk3_y2)
    Erro_dirk3_list.append(h*Erro_1_norm_dirk3_y2)
    h_List.append(h)

    Erro_rk3=np.array(Erro_rk3_list)
    Erro_dirk3=np.array(Erro_dirk3_list)
    h_array=np.array(h_List)

#%% 

## 3.3 Figure1:the convergence rate of rk3 in the moderately stiff case
# Find the best fit to the data: if err ~ A h^s then the best fit straight line of log(err) vs log(h) has slope s.
#p = np.polyfit(np.log(h_List), np.log(Erro_rk3_list), 1)
p = np.polyfit(np.log(h_List[1:]), np.log(Erro_rk3_list[1:]), 1) # exclude the outlier
plt.loglog(h_array,Erro_rk3,'rx',label='Numerical data')
plt.loglog(h_array,np.exp(p[1])*h_array**(p[0]),'b-',label="Line slope {:.2f}".format(p[0]))
plt.xlabel('h')
plt.ylabel('$\|$ Error $\|$')
plt.title("Figure1:the convergence rate of rk3 in the moderately stiff case")
plt.legend(loc = 2)
plt.grid()
plt.show()

#%% 

## 3.4 Figure2:the highest resolution (N=400) of RK3 in the moderately stiff case
# rk3 highest resolution
N=400
x_rk3,y_rk3=rk3(A,bvector,y0,interval,N)
y1_rk3=y_rk3[0,:]
y2_rk3=y_rk3[1,:]

x=x_rk3
y1_exact=np.exp(-a1*x)
y2_exact=(a1/(a1-a2))*(np.exp(-a2*x)-np.exp(-a1*x))

fig,axs = plt.subplots(1,2)
fig.suptitle('Figure2:the highest resolution (N=400) of RK3 in the moderately stiff case')

axs[0].semilogy(x,y1_exact,label='Exact') # To show the rapid decay of y1
axs[0].scatter(x,y1_rk3,s=5,color='red',label='RK3(N=400)')
axs[0].set_xlabel('x')
axs[0].set_ylabel('$y_{1}$')
axs[0].legend()
axs[0].grid()


axs[1].plot(x,y2_exact,label='Exact')
axs[1].scatter(x,y2_rk3,s=5,color='red',label='RK3(N=400)')
axs[1].set_xlabel('x')
axs[1].set_ylabel('$y_{2}$')
axs[1].legend()
axs[1].grid()

fig.tight_layout()
plt.show()

#%% 

## 3.5 Figure3:the convergence rate of dirk3 in the moderately stiff case
p = np.polyfit(np.log(h_List), np.log(Erro_dirk3_list), 1)
plt.loglog(h_array,Erro_dirk3,'rx',label='Numerical data')
plt.loglog(h_array,np.exp(p[1])*h_array**(p[0]),'b-',label="Line slope {:.2f}".format(p[0]))
plt.xlabel('h')
plt.ylabel('$\|$ Error $\|$')
plt.title("Figure3:the convergence rate of dirk3 in the moderately stiff case")
plt.legend(loc = 2)
plt.show()

#%% 

## 3.6 Figure4:the highest resolution (N=400) of DIRK3 in the moderately stiff case
N=400
x_dirk3,y_dirk3=rk3(A,bvector,y0,interval,N)
y1_dirk3=y_dirk3[0,:]
y2_dirk3=y_dirk3[1,:]

x=x_dirk3
y1_exact=np.exp(-a1*x)
y2_exact=(a1/(a1-a2))*(np.exp(-a2*x)-np.exp(-a1*x))

fig,axs = plt.subplots(1,2)
fig.suptitle('Figure4:the highest resolution (N=400) of DIRK3 in the moderately stiff case')

axs[0].semilogy(x,y1_exact,label='Exact') # To show the rapid decay of y1
axs[0].scatter(x,y1_dirk3,s=5,color='red',label='RK3(N=400)')
axs[0].set_xlabel('x')
axs[0].set_ylabel('$y_{1}$')
axs[0].legend()
axs[0].grid()


axs[1].plot(x,y2_exact,label='Exact')
axs[1].scatter(x,y2_dirk3,s=5,color='red',label='RK3(N=400)')
axs[1].set_xlabel('x')
axs[1].set_ylabel('$y_{2}$')
axs[1].legend()
axs[1].grid()

fig.tight_layout()
plt.show()

#%% 

###############################################################
##                        The  Stiff Case                   ##
###############################################################


def bvector_system_2(x):
    """
    b vector definition in Task 4.

    Parameters
    ----------
    x : float
        Coordinate

    Returns
    -------
    b : array of float
        b as given by equation (13)
    """

    # Define vector b according to equation (13)
    b=np.array([
        np.cos(10*x)-10*np.sin(10*x),
        199*np.cos(10*x)-10*np.sin(10*x),
        208*np.cos(10*x)+10000*np.sin(10*x)
    ])
    # Return vector b
    return b

# %% 
# 4.1 Input system_2
A=np.array([
    [-1,        0,      0],
    [-99,    -100,      0],
    [-10098, 9900, -10000]
])

bvector=bvector_system_2

y0=[0,1,0]

interval=[0,1]

#%% 

## 4.2 Calculate erroer
k_size=13
#Erro_rk3_list=[]
Erro_dirk3_list=[]
h_List=[]

for k in range(4,k_size+4):
    N=200*k

    #x_rk3,y_rk3=rk3(A,bvector,y0,interval,N)
    #y1_rk3=y_rk3[0,:]
    #y2_rk3=y_rk3[1,:]
    #y3_rk3=y_rk3[2,:]

    x_dirk3,y_dirk3=dirk3(A,bvector,y0,interval,N)
    #y1_dirk3=y_dirk3[0,:]
    #y2_dirk3=y_dirk3[1,:]
    y3_dirk3=y_dirk3[2,:]


    x=x_dirk3 #Since x_rk3=x_dirk3
    #y1_exact=np.exp(-a1*x)
    y3_exact=np.sin(10*x)+2*np.exp(-x)-np.exp(-100*x)-np.exp(-10000*x)

    #Erro_1_norm_rk3=0
    Erro_1_norm_dirk3_y3=0
    h=x[1]-x[0] #stepsize
    #start from j=2 in python start from j=1 since y2==y2_exact
    for j in range(1,len(y3_exact)):
        #Erro_1_norm_rk3+=abs((y2_rk3[j]-y2_exact[j])/y2_exact[j])
        Erro_1_norm_dirk3_y3+=abs((y3_dirk3[j]-y3_exact[j])/y3_exact[j])

    #Erro_rk3_list.append(h*Erro_1_norm_rk3)
    Erro_dirk3_list.append(h*Erro_1_norm_dirk3_y3)
    h_List.append(h)

    #Erro_rk3=np.array(Erro_rk3_list)
    Erro_dirk3=np.array(Erro_dirk3_list)
    h_array=np.array(h_List)

#%% 

## 4.3 Figure5:the convergence rate of dirk3 in the stiff case
p = np.polyfit(np.log(h_List[1:]), np.log(Erro_dirk3_list[1:]), 1)# Exclude the outlier
plt.loglog(h_array,Erro_dirk3,'rx',label='Numerical data')
plt.loglog(h_array,np.exp(p[1])*h_array**(p[0]),'b-',label="Line slope {:.2f}".format(p[0]))
plt.xlabel('h')
plt.ylabel('$\|$ Error $\|$')
plt.title("Figure5:the convergence rate of dirk3 in the stiff case")
plt.legend(loc = 2)
fig.tight_layout()
plt.show()

#%% 

## 4.4 Figure6:the highest resolution (N=3200) of DIRK3 in the stiff case
N=16*200
x_dirk3,y_dirk3=dirk3(A,bvector,y0,interval,N)
y1_dirk3=y_dirk3[0,:]
y2_dirk3=y_dirk3[1,:]
y3_dirk3=y_dirk3[2,:]

x=x_dirk3

y1_exact = np.cos(10*x) - np.exp(-x)
y2_exact = np.cos(10*x) + np.exp(-x) - np.exp(-100*x)
y3_exact = np.sin(10*x) + 2*np.exp(-x) - np.exp(-100*x) - np.exp(-10000*x)

fig,axs=plt.subplots(3,1)

fig.suptitle('Figure6:the highest resolution (N=3200) of DIRK3 in the stiff case')

axs[0].plot(x,y1_exact,label='Exact') # To show the rapid decay of y1
axs[0].scatter(x,y1_dirk3,s=1,color='red',label='DIRK3(N=3200)')
axs[0].set_xlabel('x')
axs[0].set_ylabel('$y_{1}$')
axs[0].legend()
axs[0].grid()


axs[1].plot(x,y2_exact,label='Exact')
axs[1].scatter(x,y2_dirk3,s=1,color='red',label='DIRK3(N=3200)')
axs[1].set_xlabel('x')
axs[1].set_ylabel('$y_{2}$')
axs[1].legend()
axs[1].grid()

axs[2].plot(x,y3_exact,label='Exact')
axs[2].scatter(x,y3_dirk3,s=1,color='red',label='DIRK3(N=3200)')
axs[2].set_xlabel('x')
axs[2].set_ylabel('$y_{3}$')
axs[2].legend()
axs[2].grid()

fig.tight_layout()
plt.show()

#%% 

## 4.5 Figure7:the highest resolution (N=3200) of RK3 in the stiff case
N=16*200
x_rk3,y_rk3=rk3(A,bvector,y0,interval,N)
y1_rk3=y_rk3[0,:]
y2_rk3=y_rk3[1,:]
y3_rk3=y_rk3[2,:]

x=x_rk3

y1_exact = np.cos(10*x) - np.exp(-x)
y2_exact = np.cos(10*x) + np.exp(-x) - np.exp(-100*x)
y3_exact = np.sin(10*x) + 2*np.exp(-x) - np.exp(-100*x) - np.exp(-10000*x)

fig,axs=plt.subplots(3,1)

fig.suptitle('Figure7:the highest resolution (N=3200) of RK3 in the stiff case')

axs[0].plot(x,y1_exact,label='Exact')
axs[0].scatter(x,y1_rk3,s=1,color='red',label='RK3(N=3200)')
axs[0].set_xlabel('x')
axs[0].set_ylabel('$y_{1}$')
axs[0].legend()
axs[0].grid()


axs[1].plot(x,y2_exact,label='Exact')
axs[1].scatter(x,y2_rk3,s=1,color='red',label='RK3(N=3200)')
axs[1].set_xlabel('x')
axs[1].set_ylabel('$y_{2}$')
axs[1].legend()
axs[1].grid()

axs[2].plot(x,y3_exact,label='Exact')
axs[2].scatter(x[:-2374],y2_rk3[:-2374],s=1,color='red',label='RK3(N=3200)')#after y2_rk3[-2374],val=NAN
axs[2].set_xlabel('x')
axs[2].set_ylabel('$y_{3}$')
axs[2].legend()
axs[2].grid()

fig.tight_layout()
plt.show()

###############################################################
##                           End                             ##
###############################################################

