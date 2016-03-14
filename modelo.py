import matplotlib
import numpy

import numpy
input2014=numpy.loadtxt('Input2014.txt')
print input2014[0]

#First,we define the value of h,the corresponding range in x and the constants
#If you need more information about the parameters of our model please visit our Wiki. 
h=0.200 #step value
a=90.000   #important, adjusted value 
bXR=300.000 #No va 
kXR=3000.000 #typical value 
gXR=0.060 #real, checked 
dXRa=0.1998 #real, important not checked 
aX=0.004691138 #real, basal checked, important
bX=300.000 #typical value 
kc=0.000001406334639 #real, important 
nX=1.000 #Assumed ...why?
kx=1409.148 #not important 
gCY=0.00071666667 #important checked, strange value  
kXRa=3000.000 #esto no va 
gXRa=0.0231 #not important not checked  


n_points = int((4000.0+h)/h)
t =[0]*(n_points)
XR =[0]*(n_points)
CY =[0]*(n_points) 
#A =[0]*(n_points)
XRa =[0]*(n_points)
#A no cambia 

#we define all the functions 
def func_XRa_prime(t,XRa,CY,XR,input2014):
    return kc*XR*input2014-(gXRa*XRa)-(dXRa*XRa)

def func_XR_prime(t,XRa,CY,XR,input2014):
    return a-(kc*input2014*kXR)-(gXR*XR)+(dXRa*XRa)

def func_CY_prime(t,XRa,CY,XR,input2014):
    return aX+(kc*input2014*XR)-(gCY*CY)
    
#First we initialize the arrays
t[0] = 0.0
XRa[0] = 0.0
XR[0] = 0.0
CY[0] = 0.0
#A[0] = 0.0

######EL INPUT2014 NO EMPIEZA EN CERO ####

#then we determine the input function

#for i in (0,n_points):
	#A[i]= input2014[i]
	
#for i in range(400):
   # A[i]=0.0
#for i in range(401,n_points):
   # A[i]=5000*sin(t[i]/100)+5000
    
    
#And finally we do 4th order Runge Kutta. 
for i in range(1,n_points):
    
    k1_XRa = func_XRa_prime(t[i-1],XRa[i-1],CY[i-1],XR[i-1],input2014[i-1])
    k1_CY = func_CY_prime(t[i-1],XRa[i-1],CY[i-1],XR[i-1],input2014[i-1])
    k1_XR = func_XR_prime(t[i-1],XRa[i-1],CY[i-1],XR[i-1],input2014[i-1])
    
    #first step
    t1 = t[i-1] + (h/2.0)
    XRa1 = XRa[i-1] + (h/2.0) * k1_XRa
    CY1 = CY[i-1] + (h/2.0) * k1_CY
    XR1 = XR[i-1] + (h/2.0) * k1_XR

    k2_XRa = func_XRa_prime(t1,XRa1,CY1,XR1,input2014[i-1])
    k2_CY = func_CY_prime(t1,XRa1,CY1,XR1,input2014[i-1])
    k2_XR = func_XR_prime(t1,XRa1,CY1,XR1,input2014[i-1])
    
    
    #second step
    t2 = t[i-1] + (h/2.0)
    XRa2 = XRa[i-1] + (h/2.0) * k2_XRa
    CY2 = CY[i-1] + (h/2.0) * k2_CY
    XR2 = XR[i-1] + (h/2.0) * k2_XR
     
    k3_XRa = func_XRa_prime(t2,XRa2,CY2,XR2,input2014[i-1])
    k3_CY = func_CY_prime(t2,XRa2,CY2,XR2,input2014[i-1])
    k3_XR = func_XR_prime(t2,XRa2,CY2,XR2,input2014[i-1])
    
    
    #third step
    t3 = t[i-1] + h
    XRa3 = XRa[i-1] + (h/2.0) * k3_XRa
    CY3 = CY[i-1] + (h/2.0) * k3_CY
    XR3 = XR[i-1] + (h/2.0) * k3_XR
    
    k4_XRa = func_XRa_prime(t3,XRa3,CY3,XR3,input2014[i-1])
    k4_CY = func_CY_prime(t3,XRa3,CY3,XR3,input2014[i-1])
    k4_XR = func_XR_prime(t3,XRa3,CY3,XR3,input2014[i-1])
    
    
     #fourth step
    average_k_XRa = (1.0/6.0)*(k1_XRa + 2.0*k2_XRa + 2.0*k3_XRa + k4_XRa)
    average_k_CY = (1.0/6.0)*(k1_CY + 2.0*k2_CY + 2.0*k3_CY + k4_CY)
    average_k_XR = (1.0/6.0)*(k1_XR + 2.0*k2_XR + 2.0*k3_XR + k4_XR)

    t[i] = t[i-1] + h
    XRa[i] = XRa[i-1] + h * average_k_XRa
    CY[i] = CY[i-1] + h * average_k_CY
    XR[i] = XR[i-1] + h * average_k_XR
    
  #  if (i < 400):
   #     A[i]=0.0
 #   if (401<i<n_points):
  #      A[i]=5000*sin(t[i]/100)+5000

#Plot

plot(t,XRa,label='XRa')
plot(t,CY,label='CY')
plot(t,XR,label='XR')
plot(t,Input2014,label='Input2014')
plt.xlabel('t(min)',size=20)
plt.ylabel('proteins',size=20)
legend()
figure(figsize(15,10))
plt.savefig("Shewanelladet.png", format='png',bbox_inches='tight',transparent=False)


print max(CY)