import pandas as pd #pandas is only used for extracting data from csv file provided. 
import cmath #used as calculations in fourier transfrom involve complex domains.
import matplotlib.pyplot as plt #used only for plotting purpose

data = pd.read_csv (r'C:\Users\yashn\Desktop\EE\Final SNS\data.csv') #mention the path for csv file here
x_n = data['x[n]'].tolist()
y_n = data['y[n]'].tolist()
h = [1/16, 4/16, 6/16, 4/16, 1/16 ]
l = [0.1,0.1,0.2,0.2,0.2,0.1,0.1]

#Please note that all values in frequency domain are stored as a dictionary with keys as angle starting from -π to +π.
#The difference between adjacent value of frequencies(keys/angles) is 0.01.
#Please note all values in n domain are stored in lists with 0 value at 0th index and kth value at kth index of list.

def plots(lst,lst2,lst3,lst4): # plotting 4 lists values v/s index of list
    x = range(len(lst))   
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    ax.plot(x, lst, color ='tab:blue',label='y_n')
    ax.plot(x, lst2, color ='tab:orange',label='x_n')
    ax.plot(x, lst3, color ='tab:green',label='x1[n] ')
    ax.plot(x, lst4, color ='tab:cyan',label='x2[n]')
    ax.legend()
    plt.xlabel("Value of n")
    plt.ylabel("Corresponding Values of Signals")
    plt.show()

def average(lst): # getting average value of list
    sum=0
    for i in lst:
        sum+=i
    return sum/len(lst)

def convolution(x,h): #Convolution of 2 signals x and h Assuming len(h)<len(x)
    y1_n = []
    temp = [0]*(len(h)//2) + x + [0]*(len(h)//2) #Since len(h)=k, adding k//2 zeros at right and left, so that convolution sum is evaluated correctly at end points
    temp2 = h[::-1] # reversing h, ie h[k]=h[-k]
    for i in range(len(h)//2,len(x)+len(h)//2):  #so as to exclude 0s added at end and start 
        y_i = 0
        for j in range(i - len(h)//2, i + len(h)//2+1): #shifting h[0] to h[i] and considering only values of non-zero h
            y_i+=temp[j]*temp2[j-i+len(h)//2]  #as temp2 has k values of h, we scale it down to [0,k-1] 
        y1_n.append(y_i)    
    return y1_n       #returning a list of convolution values
   
def fouriertrans(x_n,num): #num corrwsponds to 0th value in index of x_n ie x[0]
    freq={}
    angle = -cmath.pi 
    while angle <= cmath.pi:  
        X_angle = 0
        for k in range(len(x_n)): #calculating sum for a particluar angle ie frequency
            
            z = complex(0, angle*(k-num)*(-1)); # complex(x,y)=x+jy
            
            X_angle += x_n[k]*cmath.exp(z)  
        freq[angle] = X_angle
        angle+=0.01
    return freq # its a dictionary with angle:value of DTFT at that angle

def ratio_of_FT(n,d): #to calculate ratios of fourier, Y/X, Y/H and so on
    value={}
    angle=-cmath.pi
    while angle <=cmath.pi:
        value[angle] = n[angle]/d[angle] 
        angle+=0.01
    return value

def inverseF(dct, number): #To calculate inverse Fourier Transform of dct, with sequence going from 0 to number(num) specified
    lst=[]
    for n in range(number):
        angle=-cmath.pi
        f=0
        while angle <=cmath.pi:
            z=complex(0,angle*n)   
            f+=dct[angle]*cmath.exp(z)*0.01            
            angle+=0.01
        lst.append(0.5*f/cmath.pi)
    return lst    #returns a signal in n domain in the form of list


def ReciprocalFT(dct, freq): #To calculte reciprocal of a given Fourier Transform dct
    # freq corresponds to frequency after which 1/dct may tend to overshoot to high values, which may lead to errors and noise addtion.
    Z1={}
    angle = -cmath.pi
    while angle <= cmath.pi:
        if abs(angle) <= freq:
            Z1[angle]= 1/dct[angle]      
        else:
            Z1[angle]= 0.84
        angle+=0.01
    return Z1

def Products_of_FT(dct1, dct2): #product of 2 Fourier Tranforms dct1 and dct2
    Z1={}
    angle = -cmath.pi
    while angle <= cmath.pi:
        
        Z1[angle]= dct1[angle]*dct2[angle]
        angle+=0.01
    return Z1

def error(x_n, m_n): #Calculating RMS value of error
    #x_n is base of comparision, m_n is new array
    e_n=0
    for i in range(0,len(x_n)):
        e_n+=(x_n[i]-m_n[i])**2
    return abs((e_n/len(x_n)))

def remove_term(x_n, m_n): #calculate RMS error function by excluding 3 terms at start and 4 terms at end. 
    e_n=0
    for i in range(3,len(x_n)-4):
        e_n+=(x_n[i]-m_n[i])**2
    return abs((e_n/(len(x_n)-7)))

def analysis(lst): #in order to curb any unwanted high rise/fall in values of lst, which may help to diminish noise
    lst= [lst[0]]*2+lst+ [lst[-1]]*2
    a=0.5
    for i in range(2,len(lst)-2):    
        if abs(lst[i])>a + (abs(0.5*lst[i+2]+0.5*lst[i-2])):
            lst[i] = min((0.5*lst[i-2]+0.5*lst[i+2]),(0.5*lst[i-1]+0.5*lst[i+1])) +a
        elif abs(lst[i]) <  (abs(0.5*lst[i-2]+0.5*lst[i+2])) - a : #abs(0.5*lst[i-1]+0.5*lst[i+1])
            lst[i] = max((0.5*lst[i-2]+0.5*lst[i+2]),(0.5*lst[i-1]+0.5*lst[i+1])) -a
    return lst[2:len(lst)-2]

def Rounding(lst): #Limit all the values of list to 4 decimal places, as in given data. 
    lst2=[]
    for k in lst:
        lst2.append(round(k,4))
    return lst2

    
Y, H, Denoi = fouriertrans(y_n,0), fouriertrans(h,len(h)//2), fouriertrans(l,len(l)//2)
#Y is Fourier Transform of y_n, H is fourier transform of h_n, Denoi is Fourier Transform of Low pass filter used for Denoising l
Y2= fouriertrans([average(y_n)]*(len(l)//2)+y_n+[average(y_n)]*(len(l)//2),0)
for k in Denoi: #Since denoi is a low pass filter, making denoise frequencies 0 above a threshold frequency
    if abs(k)>2:
        Denoi[k]=0.0 
G = ReciprocalFT(H , 1.55) #G=1/H


#Temp corresponds to step of x_2[n], ie first deblur then denoise
Temp1 = Products_of_FT(Y,G)
Temp1 = inverseF(Temp1, len(y_n))
Temp1 = [average(Temp1)] * (len(l)//2) + Temp1 + [average(Temp1)]*(len(l)//2)
Temp2 = fouriertrans(Temp1,0)
Temp2 = Products_of_FT(Temp2, Denoi)
Temp3 = inverseF(Temp2,len(y_n)+len(l)//2*2)
Temp3 = Temp3[len(l)//2:len(y_n)+len(l)//2]
Temp4=[]
for k in Temp3:
    Temp4.append(abs(k))
Temp4 = analysis(Temp4) 
Xn_2 = Rounding(Temp4)


#Zemp corresponds to steps of x_1[n], ie first denoise then deblur
Zemp1 = Products_of_FT(Y2, Denoi)
Zemp1 = inverseF(Zemp1, len(y_n)+len(l)//2*2)
Zemp1 = Zemp1[len(l)//2:len(y_n)+len(l)//2]
Zemp2 = fouriertrans(Zemp1, 0)
Zemp2 = Products_of_FT(Zemp2, G)
Zemp3 = inverseF(Zemp2,len(y_n))
Zemp4=[]
for k in Zemp3:
    Zemp4.append(abs(k))
Zemp4 = analysis(Zemp4)
Xn_1 = Rounding(Zemp4)
  

e1= error(x_n,Xn_1)
e2 = error(x_n, Xn_2)
e3= error(x_n,y_n)
e4= remove_term(x_n,Xn_2)
e5= remove_term(x_n, Xn_1)

print('x_1[n]= ',Xn_1)
print('x_2[n]= ',Xn_2)
print('The original MSE error was',e3)
print('The MSE error when first denoisied i.e.error in X_1[n]',e1)
print('The MSE error excluding 3 start terms and 4 end terms when first denoisied i.e. error in X_1[n]',e5)
print('The MSE error when first deblur i.e. error in X_2[n]',e2)
print('The MSE error excluding 3 start and 4 end terms temrs when first deblur i.e. error in X_2[n]',e4)


plots(y_n,x_n,Xn_1,Xn_2) #plotting 4 lists 


    
    
    

     

            
            
            
 