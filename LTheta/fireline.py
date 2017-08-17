import numpy as np
import scipy.integrate as integrate
import scipy.interpolate as interp
from matplotlib import pyplot as plt
import os, sys


def star(r,points,mag):
    x = np.zeros(N)
    y = np.zeros(N)
    d = r + mag * np.cos(points*alpha)
    x = np.cos(alpha)*d
    y = np.sin(alpha)*d
    return x,y

def findNeighbors(xs,ys):
    #iterate over all x,y pairs
    zI = xs[lag,:]+1j*ys[lag,:]
    zn = np.fft.fft(zI)
    # z0 = xs[0,:]+1j*ys[0,:]
    neighborsX = np.zeros(len(xs[0,:]))
    neighborsY = np.zeros(len(xs[0, :]))
    modes = np.linspace(-N / 2, N / 2 - 1, N)
    modes = np.fft.ifftshift(modes)
    for i in range(len(xs[0,:])):
        guess = 2*np.pi/N*np.argmin(((xs[0,i]-xs[lag,:])**2+(ys[0,i]-ys[lag,:])**2)**(0.5))
        g = 1
        z0 = xs[0,i]+1j*ys[0,i]
        exponentials = np.exp(1j * modes * guess)
        z = np.sum(zn * exponentials) / N
        k=0
        while(abs(g)>1E-6):
            k+=1
            exponentials = np.exp(1j*modes*guess)
            z = np.sum(zn*exponentials)/N
            zP = np.sum(1j*modes*zn*exponentials)/N
            zPP = -np.sum((modes**2)*zn*exponentials)/N
            g = np.real((z-z0)*np.conj(zP))
            if(k>5):
                numberFile=open("NumberLog.dat","r")
                logNumber = int(numberFile.readline())
                print("NUMBER",logNumber)
                numberFile.close()
                #Now increase number and rewrite numberlog
                numberFile=open("NumberLog.dat","w")
                numberFile.write(str(logNumber+1))
                numberFile.close()
                """
                   File Structure:
                         N Points
                         Initial Theta Guess
                         Outer Point
                         Inner Geometry Points (N)
                  """
                logFile = open("FN"+str(k)+"Threshold"+str(logNumber+1),"wb")
                logFile.write(str.encode(str(N)+"\n"))
                logFile.write(str.encode(str(i*2*np.pi/N)+"\n"))
                logFile.write(str.encode(str(xs[0,i])+","+str(ys[0,i])+"\n"))
                for j in range(N):
                    logFile.write(str.encode(str(xs[lag,j])+","+str(ys[lag,j])+"\n"))
                logFile.close()
                """PLOTTING ISSUE
                  plt.figure(7)
                  # plt.plot(np.real(zI),np.imag(zI))
                  # plt.scatter(np.real(z0),np.imag(z0))
                  # plt.scatter(np.real(z),np.imag(z))
                  zed = xs[lag, :] + 1j * ys[lag, :]
                  plt.plot(alpha, np.real((z0 - zed) * (Spectral_Derivative(zed, False))))
                  zeroes = np.zeros(len(alpha))
                  plt.plot(alpha,zeroes)
                  plt.show()
                  """
                print("*********NEIGHBOR ISSUE***********")
                break
            gP = np.real((z-z0)*np.conj(zPP))+abs(zP)**2
            guess-=g/gP
        exponentials = np.exp(1j*modes*guess)
        zNearest = np.sum(zn*exponentials)/N
        neighborsX[i]=np.real(zNearest)
        neighborsY[i]=np.imag(zNearest)
    return neighborsX, neighborsY

########################
def bgFlow(x,y):
    # background flow
    # u = 1*np.cos(7*alpha)
    u = 0.0*np.ones(len(x))
    v = 8.0*np.ones(len(x))
    # u = 0*np.ones(len(x))#+2*np.random.uniform(-1,1)
    # u = np.cos(3*alpha)
    # v = 0*np.ones(len(x))#+4.*np.random.uniform(-1,1)
    # u=np.ones(len(x))
    # v=np.ones(len(x))

    return u,v
########################
def DrawSink(midPointX,midPointY,xWidth,yWidth,density):
    [xp,yp]=np.meshgrid(np.linspace(-xWidth/2,xWidth/2,density),np.linspace(-yWidth/2,yWidth/2,density))
    sinkx = np.zeros((len(xp[:,0]),len(xp[0,:])))
    sinky = np.zeros((len(yp[:,0]),len(yp[0,:])))
    for i in range(len(midPointX)):
        dist2 = (xp-midPointX[i])**2 + (yp-midPointY[i])**2
        sinkx+= (midPointX[i]-xp)*sa[i]/(dist2+eps2/2)
        sinky+= (midPointY[i]-yp)*sa[i]/(dist2+eps2/2)
        #SINK PLOTTING
        # plt.plot(xp[density/2-2,density/2-2],yp[density/2-2,density/2-2],'r*')
        # plt.quiver(xp,yp,sinkx,sinky)
        # plt.scatter(midPointX,midPointY)

def ThetaAdjust(theta):
    if (theta[0] < 0):
        theta[0] = 2 * np.pi - abs(theta[0])
    for k in range(1,N):
        if (theta[k] < 0):
            theta[k] = 2 * np.pi - abs(theta[k])
        while (abs(theta[k-1]-theta[k])>np.pi):
            theta[k]+=2*np.pi

def betaFunction(adjTheta):
    beta = np.zeros(N)
    beta[0]=adjTheta[0].copy()
    dtheta = (adjTheta[N-1]-beta[0])/(N-1)
    for i in range(1,N):
        beta[i]=beta[0]+i*dtheta
    return beta

def CalcBetaCubic(theta):
    A = np.zeros((3,3))
    b = np.zeros(3)
    indexes = [N/4-1,3*N/4-1,N-1]
    for i in range(3):
        index = indexes[i]
        b[i]=theta[index]
        if(i==2):
            b[i]=b[i]/2
        for j in range(3):
            A[i,j]=theta[index]**(3-j)
    coeff = np.linalg.solve(A,b)
    return coeff

def Calc_Vel_Tan(theta,vel_normal,L):
    # plt.figure(4)
    # plt.plot(alpha,theta-alpha)
    # thetab = theta-alpha
    # print(thetab[0],thetab[1],thetab[-1])
    # plt.figure(1)
    thetaPrime = np.real(Spectral_Derivative(theta-alpha,True))
    modes = np.linspace(-N / 2, N / 2 - 1, N)
    modes = np.fft.ifftshift(modes)
    thetaPH = np.fft.fft(thetaPrime)
    # thetaPH = np.fft.fftshift(thetaPH)
    # plt.figure(4)
    # plt.clf()
    # plt.plot(alpha,thetaPrime,'b-o')
    # plt.figure(1)
    #Add derivative of linear beta function
    thetaPrime += 1.
    dTheta = 2*np.pi/N
    avgVnThetaPrime = np.sum(thetaPrime * vel_normal) * 2 * np.pi / N
    avgVnThetaPrime = avgVnThetaPrime / (2 * np.pi )
    vel_tan = np.real(Spectral_Integration(thetaPrime*vel_normal))
    vnPrime = np.real(Spectral_Derivative(vel_normal,False))
    vel_tan -=alpha*avgVnThetaPrime
    theta -= dt / L * (vnPrime + thetaPrime * vel_tan)
    if (VERBOSE):
        print("<Vn*ThetaPrime>:", avgVnThetaPrime)
        print("DL",2*np.pi*dt*avgVnThetaPrime)
    L+= 2*np.pi*dt* avgVnThetaPrime
    return theta, vel_tan, L

def Spectral_Integration(function):
    N = len(function)
    modes = np.linspace(-N / 2, N / 2 - 1, N)
    modes = np.fft.ifftshift(modes)
    xn = np.fft.fft(function)
    xn[1:] = xn[1:] / (1j * modes[1:])
    x0 = xn[0].copy()
    xn[0] = (-1)* np.sum(xn[1:])
    return np.fft.ifft(xn)+x0*np.linspace(0,2*np.pi,N,endpoint=False)/N

def GaussianFilter(modes,c):
    return np.exp(-modes**2/(2*c**2))

def Spectral_Derivative(function,PLOT):
    N=len(function)
    modes = np.linspace(-N / 2, N / 2 - 1, N)
    xn = np.fft.fft(function)
    xn = np.fft.fftshift(xn)
    #xn = GaussianFilter(modes, 10.) * xn
    if(PLOT and PLOTTING):
        plt.subplot(3, 2,5)
        plt.semilogy(modes,np.abs(xn))
        plt.ylim(ymin=1E-16)
        plt.title("FFT(Theta-Beta)")
        plt.subplot(3, 2,6)
        plt.semilogy(modes,np.abs(1j*modes*xn))
        plt.ylim(ymin=1E-16)
        plt.title("i*modes*FFT(Theta-Beta)")
        plt.subplot(3, 2, 2)
        plt.title("Theta")
        plt.plot(alpha, theta)
    xn = np.fft.ifftshift(xn)
    modes = np.fft.ifftshift(modes)
    return np.fft.ifft(1j*modes*xn)

def Calc_Positions(xs,ys,L,theta,vn,vt,dt):
    xP = L/(2*np.pi)*np.real(Spectral_Integration(np.cos(theta)))
    yP = L/(2*np.pi)*np.real(Spectral_Integration(np.sin(theta)))
    xP_Bar = np.sum(xP)/ N
    yP_Bar = np.sum(yP) / N
    # avgVnThetaPrime = np.sum(thetaPrime * vel_normal) * 2 * np.pi / N
    x_Bar = np.sum(xs)/N
    y_Bar = np.sum(ys)/N
    cX = np.sum((vn*np.sin(theta))+vt*-np.cos(theta))/N
    cY = np.sum((vn*(-np.cos(theta)))+vt*-np.sin(theta))/N
    xShift = x_Bar-xP_Bar+cX*dt
    yShift = y_Bar-yP_Bar+cY*dt
    xs = xP + xShift
    ys = yP + yShift
    return(xs,ys)

def HeatFlux(xs,ys,neighborsX,neighborsY,vel_normal,scaleFactor):
    MAX_WIDTH = 10
    thickness = np.sqrt((xs[0,:]-neighborsX)**2+(ys[0,:]-neighborsY)**2)
    plt.figure(1)
    speed = np.sqrt(np.real(Spectral_Derivative(xs[0,:],False))**2 + np.real(Spectral_Derivative(ys[0,:],False))**2)
    for k in range(0, N):
        dist2 = (xs[0, k] - xs[0, :]) ** 2 + (ys[0, k] - ys[0, :]) ** 2
        thickness=1/np.sqrt((xs[0,:]-neighborsX)**2+(ys[0,:]-neighborsY)**2)
        vel_normal[k] = scaleFactor  * np.sum(thickness / np.sqrt(dist2 + eps2) * speed) * 2 * np.pi / N
    #SINK EFFECT
    # vel_normal = vel_normal*(1-thickness/MAX_WIDTH)
    return vel_normal

def RedistributePoints(x,y):
    #Take fft of initial shape
    z = x + 1j*y
    modes = np.linspace(-N / 2, N / 2 - 1, N)
    modes = np.fft.ifftshift(modes)
    zh = np.fft.fft(z)
    dzh = 1j*modes*zh
    dz = np.fft.ifft(dzh)
    dx = np.real(dz)
    dy = np.imag(dz)
    speed = np.sqrt(dx**2+dy**2)
    #Calculate Total Arclength
    L = np.pi * 2 / N * np.sum(speed)
    arclength = np.real(Spectral_Integration(speed))
    #Pad left and right side since periodic (Help with Spline interp)
    arcL = arclength[-7:-1]-L+arclength[1]
    arcR = arclength[0:6]+L
    arclength = np.concatenate((arcL,arclength))
    arclength = np.concatenate((arclength,arcR))
    alphaR = alpha[0:6]+2*np.pi
    extAlpha = np.concatenate(((alpha[-7:-1] - 2 * np.pi+alpha[1]), alpha))
    extAlpha = np.concatenate((extAlpha,alphaR))
    #CUBIC SPLINE INTERPOLATION
    coeff = interp.splrep(arclength,extAlpha,s=0)
    newAlpha = interp.splev(np.linspace(0,N-1,N)*L/N,coeff,der=0)
    # plt.semilogy(range(len(newAlpha)),abs(alpha-newAlpha),'o')
    znew = np.complex(0)
    #EVALUATE AT NEW ALPHA
    for i in range(N):
        znew+= zh[i]/N*np.exp(1j*modes[i]*newAlpha)

    xold = x.copy()
    yold = y.copy()
    x = np.real(znew)
    y = np.imag(znew)
    # plt.semilogy(range(len(x)),abs(x-xold),'r')
    # plt.semilogy(range(len(x)),abs(y-yold))
    plt.show()
    return x,y

def Convection(xs,ys,neighborsX,neighborsY):
    eps = 0.25
    midpointsX = (xs+neighborsX)/2
    midpointsY = (ys+neighborsY)/2
    sinkx = np.zeros(len(xs))
    sinky = np.zeros(len(ys))
    for i in range(len(xs)):
        r2= (xs[i]-midpointsX)**2 + (ys[i]-midpointsY)**2 +eps
        sinkx[i] = np.sum(-(xs[i]-midpointsX)/r2)
        sinky[i] = np.sum(-(ys[i]-midpointsY)/r2)
    print( sinkx, sinky)
    return sinkx, sinky

def UpdateLog(xs,ys,xpoint,ypoint,index,updating):
    if(updating):
        #Create Buffer
        file = open(LOG_FILENAME+".dat",'r')
        sets = int(file.readline())
        fullBuffer = []
        fullPB = []
        for i in range(sets):
            #Read N
            n = int(file.readline())
            #Read # of Trouble Points
            nIssues = int(file.readline())
            bufferx = []
            buffery = []
            pointsBuffer = []
            for j in range(n):
                bufferx.append(float(file.readline()))
                buffery.append(float(file.readline()))

def main():
    global N
    global L
    global eps2
    global theta
    global VERBOSE, PLOTTING
    global LOG_FILENAME, UPDATING
    LOG_FILENAME = "test1"
    UPDATING = False
    VERBOSE = True
    PLOTTING = True
    eps2 = 1e-1
    N = 32  # number of points on front
    T = 2    # time horizon
    m = 100  # number of time steps
    global dt
    dt = T/m # time step size
    global alpha
    alpha = np.linspace(0, 2 * np.pi, N, endpoint=False)
    global lag
    lag = 3
    """NEEDS TO BE A FUNCTION OF LAG"""
    scaleFactor = lag/dt/600
    e_x = 1
    e_y = 1
    print( "Number of points on front is %d" % N)
    print("Time Horizon:",T)
    print("Timesteps:", m)
    global modes
    modes = np.linspace(-N / 2, N / 2 - 1, N)
    modes = np.fft.ifftshift(modes)

    ###############################################################################
    #INITIAL SHAPE
    x = 5*np.cos(alpha)/e_x
    y = 5*np.sin(alpha)/e_y
    #x,y = star(1,9,0.2)
    x,y=RedistributePoints(x,y)

    xs=np.zeros((lag+1,len(x)))
    ys = np.zeros((lag+1,len(y)))
    #Set initial history as starting values
    xs[0,:]=x
    ys[0,:]=y
    #Inner Fireline
    x0 = xs[0, :].copy()
    y0 = ys[0, :].copy()
    x0=x0/1.1
    y0=y0/1.1#-0.02
    # x0 = 0.2*np.cos(alpha)+.02
    # y0 = 0.2*np.sin(alpha)-0.01
    ###############################################################################

    x0,y0=RedistributePoints(x0,y0)
    plt.plot(x,y)
    plt.plot(x0,y0)
    plt.axis('equal')
    plt.show()
    for i in range(1,lag+1):
        xs[i,:]=x0
        ys[i,:]=y0
    time = 0. # current time
    neighborsX, neighborsY = findNeighbors(xs, ys)
    # neighborsX = xs[lag,:]
    # neighborsY = ys[lag,:]
    midpointsX = (xs[0, :] + neighborsX) / 2
    midpointsY = (ys[0, :] + neighborsY) / 2

    """GET THETA VALUES IMMEDIATELY"""
    z = xs[0, :] + 1j * ys[0, :]
    zn = np.fft.fft(z)
    Dz = Spectral_Derivative(z,False)
    global sa
    sa = abs(Dz)
    Dx = np.real(Dz)
    # x-coordinate of the derivative of the shape
    Dy = np.imag(Dz)
    # y-coordinate of the derivative of the shape
    speed = np.sqrt(Dx * Dx + Dy * Dy)
    # Jacobian or arclength of the curve
    nx = Dy / speed
    # x-coordinate of outward normal
    ny = -Dx / speed
    # y-coordinate of outward normal
    L = np.pi * 2 / N * np.sum(speed)
    print("L",L)
    vel_normal = np.zeros(N)
    theta = np.arctan2(Dy, Dx)
    ThetaAdjust(theta)
    # curvature = -np.real(Spectral_Derivative(theta - alpha, False))-1
    """
    #PAPER BURNING
    conX = xs[0, int(N / 8):int(3 * N / 8 + 1)]
    conY = ys[0, int(N / 8):int(3 * N / 8 + 1)]
    #switch order because of counterclockwise alpha
    conX = np.asarray(list(reversed(conX)))
    conY = np.asarray(list(reversed(conY)))
    mean = np.mean(conY)
    max = np.max(abs(conY-mean))
    ###########
    """
    for i in range(0,m):
        plt.clf()
        time = time + dt
        if(VERBOSE):
            print("-------------------------------------")
            print("Time:",time)
            print((theta[N - 1] - theta[0]) / (2 * np.pi*(N-1)/N))
        plt.figure(1)
        velx, vely = bgFlow(xs[0, :], ys[0, :])
        #SHIFT ALL DATA FOR LAG
        for j in range(0,lag):
            xs[lag-j,:]=xs[lag-(j+1),:]
            ys[lag-j,:]=ys[lag-(j+1),:]

        neighborsX, neighborsY = findNeighbors(xs, ys)
        
        ########################################################################
        #NORMAL-VELOCITY DRIVERS

        #HEATFLUX
        #vel_normal=HeatFlux(xs,ys,neighborsX,neighborsY,vel_normal,scaleFactor)/1.1
        vel_normal = np.zeros(len(xs[0,:]))

        #BACKGROUND WIND
        #vel_normal += np.sin(theta) * velx - np.cos(theta) * vely 

        #SINKS
        sinkx, sinky = Convection(xs[0,:],ys[0,:],neighborsX,neighborsY)
        vel_normal += np.sin(theta) * sinkx - np.cos(theta) * sinky 
        ########################################################################

#   NEW PDE
        """Algorithm  (U_t=(U_x)^2 - U_xx - U_xxxx
        Calculate Ux^2
        Calculate 2nd and 4th derivative as (modes^2 - modes^4)*weights
        U = xs[0,:] + 1j*ys[0,:]
        Ux = Spectral_Derivative(U,False)
        uh = np.fft.fft(U)
        vel_normal = abs(Ux**2 + np.fft.ifft((modes**2 - modes**4)*uh))
        plt.plot(np.linspace(0,N,N),vel_normal)
        plt.show()
        """
        theta, vel_tan, L = Calc_Vel_Tan(theta,vel_normal,L)
        xs[0, :], ys[0, :]=Calc_Positions(xs[0,:],ys[0,:],L,theta,vel_normal,vel_tan,dt)

        """
        #PAPER BURNING
        hx = np.gradient(cony)
        hxx = np.gradient(hx)
        conY += dt*(0.1*hxx + 5*hx**2)
        normConY = (conY-np.mean(conY))/np.max(abs(conY-np.mean(conY)))
        sliceY = ys[0, int(N / 8):int(3 * N / 8 + 1)]
        normY = (sliceY-np.mean(sliceY))/np.max(abs(sliceY-np.mean(sliceY)))
        """

        plt.clf()
        # plt.plot(conX,normConY,'r.')
        # plt.plot(conX,normY,'b')

        # plt.figure(5)
        # plt.semilogy(conX,abs(normY-normConY))
        # plt.plot(conX,hx)
        # plt.plot(conX,hxx)
        plt.figure(1)



        if(VERBOSE):
            print("Theta Closed",(theta[-1]-theta[0]))
            print("L=", L)
            print(L/(2*np.pi*np.exp(1.)))
            print("-------------------------------------")

        """PLOTTING"""
        if(PLOTTING):
            plt.subplot(3, 2, 1)
            plt.title("Shape")
            # plt.plot(alpha,vel_normal)
            plt.plot(xs[0,:],ys[0,:],'r--o')
            plt.axis('equal')
            plt.xlim([-1.1,1.1])
            plt.ylim([-1.1,1.1])
            plt.subplot(3, 2,3)
            plt.title("FFT(Vn)")
            # plt.semilogy(np.linspace(-N/2,N/2-1,N,endpoint=False),np.abs(np.fft.fftshift(np.fft.fft(vel_normal))))
            plt.plot(alpha,vel_normal)
            plt.subplot(3, 2, 4)
            plt.title("Vs")
            plt.plot(alpha,vel_tan)

            plt.figure(2)
            plt.plot(xs[0,:],ys[0,:],'r--o')
            plt.plot(xs[lag,:],ys[lag,:],'b')
            plt.plot(xs[lag-1,:],ys[lag-1,:],'g')

            plt.scatter(neighborsX,neighborsY,marker = '*')
            # sliceY = ys[0, (N / 8):(3 * N / 8 + 1)]
            # sliceX = xs[0, (N / 8):(3 * N / 8 + 1)]

            # adjX = sliceX/(sliceX[0])*0.6
            # print("adjx",adjX)
            # adjY = sliceY-np.mean(sliceY)
            # plt.plot(adjX,adjY*0.6)
            # plt.plot(sliceX,sliceY)
            # plt.plot(xs[lag,:],ys[lag,:],'b-o')

            # ax = plt.gca()
            # ax.axes.get_xaxis().set_visible(False)
            # ax.axes.get_yaxis().set_visible(False)
            plt.axis('equal')
            # fname = '_tmp%05d.pdf' % i
            # plt.savefig(fname)

            plt.draw()
            plt.pause(1E-10)
            # plt.draw()
            # plt.pause(1E-16)
            plt.waitforbuttonpress()
    plt.show()

main()



