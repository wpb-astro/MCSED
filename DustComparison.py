import numpy as np 
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt 
from scipy import odr

def f(B,x):
    '''Linear function y = m*x + b'''
    # B is a vector of the parameters.
    # x is an array of the current x values.
    # x is in the same format as the x passed to Data or RealData.
    #
    # Return an array in the same format as y passed to Data or RealData.
    return B[0]*x + B[1]        

def plotting(xnd,xndl,xndu,xd,xdl,xdu,z,name,folder):
    xnderrl = xnd-xndl
    xnderrh = xndu-xnd
    xderrl = xd-xdl
    xderrh = xdu-xd
    #zmin,zmax = min(z),max(z)
    #colors = 2.0*(z-zmin)/(zmax-zmin) - 1.0 #Colors between -1 and 1 depending on redshift range

    linear = odr.Model(f)
    mydata = odr.RealData(xnd,xd,sx=(xnderrh+xnderrl)/2.0,sy=(xderrh+xderrl)/2.0)
    myodr = odr.ODR(mydata,linear,beta0=[1.0,0.0])
    myoutput = myodr.run()
    myoutput.pprint()

    plt.figure()
    sc = plt.scatter(xnd,xd,c=z,label='')
    #plt.errorbar(xnd,xd,yerr=np.array((xderrl,xderrh)),xerr=np.array((xnderrl,xnderrh)),fmt='',label='')
    plt.errorbar([min(xnd)],[max(xd)],yerr=[(np.average(xderrl),np.average(xderrh))],xerr=[(np.average(xnderrl),np.average(xnderrh))],fmt='',label='Typical Error')
    plt.plot(xnd,f(myoutput.beta,xnd),'r--',label=r"y=(%.3f$\pm$%.3f)x + (%.3f$\pm$%.3f)"%(myoutput.beta[0],myoutput.sd_beta[0],myoutput.beta[1],myoutput.sd_beta[1]))
    plt.xlabel("%s z=0.0077"%(name))
    plt.ylabel("%s Metallicity Free"%(name))
    plt.colorbar(sc,label="Redshift")
    #cb.set_label("Redshift")
    plt.legend(loc='best')
    plt.savefig("output/%s/%s_comp.png"%(folder,name))
    
def ploterr(xnd,xndl,xndu,xd,xdl,xdu,name,folder,folder2):
    xnderrl=xnd-xndl
    xnderrh=xndu-xnd
    xderrl=xd-xdl
    xderrh=xdu-xd
    
    xnderrmean = np.average((xnderrh+xnderrl)/2.0)
    xderrmean = np.average((xderrh+xderrl)/2.0)
    xnderrerr = np.std((xnderrh+xnderrl)/2.0)
    xderrerr = np.std((xderrh+xderrl)/2.0)
    
    plt.figure()
    plt.hist([(xnderrh+xnderrl)/2.0,(xderrh+xderrl)/2.0],bins=20,histtype='step',stacked=False,color=['blue','red'],label=[r'%s: $\sigma = %.3f \pm %.3f$'%(folder,xnderrmean,xnderrerr), r'%s: $\sigma = %.3f \pm %.3f$'%(folder2,xderrmean,xderrerr)])
    plt.xlabel("%s Error"%(name))
    plt.legend(loc='best')
    plt.savefig("output/%s/%s_error_hist"%(folder,name))

def plotting2(xd,xdl,xdu,xconst,yd,ydl,ydu,yconst,z,namex,namey,folder):
    xderrl = xd-xdl
    xderrh = xdu-xd
    yderrl = yd-ydl
    yderrh = ydu-yd

    plt.figure()
    sc = plt.scatter(xd,yd,c=z,label='')
    #plt.errorbar(xd,yd,yerr=np.array((yderrl,yderrh)),xerr=np.array((xderrl,xderrh)),fmt='',label='')
    plt.errorbar([min(xd)],[max(yd)],yerr=[(np.average(yderrl),np.average(yderrh))],xerr=[(np.average(xderrl),np.average(xderrh))],fmt='',label='Typical Error')
    plt.vlines(xconst,min(yd),max(yd),colors='r',label=r'$%s=%.3f$'%(namex,xconst))
    plt.hlines(yconst,min(xd),max(xd),colors='b',label=r'$%s=%.3f$'%(namey,yconst))
    plt.xlabel("%s"%(namex))
    plt.ylabel("%s"%(namey))
    plt.colorbar(sc,label="Redshift")
    plt.legend(loc='best')
    plt.savefig("output/%s/%svs%s_comp.png"%(folder,namey,namex))

foldernd = raw_input("Which folder/filename for x axis: ")
folderd = raw_input("Which folder/filename for y axis: ")
fieldnd,photidnd = np.loadtxt("output/%s/%s.dat"%(foldernd,foldernd),dtype=np.dtype("S6,int"),skiprows=2,usecols=(0,1),unpack=True)
#znd,sfr1nd,sfr1ndl,sfr1ndu,sfr2nd,sfr2ndl,sfr2ndu,sfr3nd,sfr3ndl,sfr3ndu,sfr4nd,sfr4ndl,sfr4ndu,sfr5nd,sfr5ndl,sfr5ndu,ebvnd,ebvndl,ebvndu,deltand,deltandl,deltandu,Ebnd,Ebndl,Ebndu,logMnd,logMndl,logMndu = np.loadtxt("output/%s/%s.dat"%(foldernd,foldernd),dtype=float,skiprows=2,usecols=(2,5,4,6,10,9,11,15,14,16,20,19,21,25,24,26,30,29,31,35,34,36,40,39,41,45,44,46),unpack=True)  ## Empirical_direct SFH, no dust, Noll
znd,sfr1nd,sfr1ndl,sfr1ndu,sfr2nd,sfr2ndl,sfr2ndu,sfr3nd,sfr3ndl,sfr3ndu,sfr4nd,sfr4ndl,sfr4ndu,ebvnd,ebvndl,ebvndu,logMnd,logMndl,logMndu = np.loadtxt("output/%s/%s.dat"%(foldernd,foldernd),dtype=float,skiprows=2,usecols=(2,5,4,6,10,9,11,15,14,16,20,19,21,25,24,26,30,29,31),unpack=True)  ## Empirical_direct SFH, no dust, Calzetti
#znd,sfrnd,sfrndl,sfrndu,agend,agendl,agendu,ebvnd,ebvndl,ebvndu,logMnd,logMndl,logMndu = np.loadtxt("output/%s/%s.dat"%(foldernd,foldernd),dtype=float,skiprows=2,usecols=(2,5,4,6,10,9,11,15,14,16,20,19,21),unpack=True) ## Constant SFR, metallicity fixed

fieldd,photidd = np.loadtxt("output/%s/%s.dat"%(folderd,folderd),dtype=np.dtype("S6,int"),skiprows=2,usecols=(0,1),unpack=True)
#zd,sfr1d,sfr1dl,sfr1du,sfr2d,sfr2dl,sfr2du,sfr3d,sfr3dl,sfr3du,sfr4d,sfr4dl,sfr4du,sfr5d,sfr5dl,sfr5du,ebvd,ebvdl,ebvdu,deltad,deltadl,deltadu,Ebd,Ebdl,Ebdu,umin,uminl,uminu,gamma,gammal,gammau,qpah,qpahl,qpahu,logMd,logMdl,logMdu = np.loadtxt("output/%s/%s.dat"%(folderd,folderd),dtype=float,skiprows=2,usecols=(2,5,4,6,10,9,11,15,14,16,20,19,21,25,24,26,30,29,31,35,34,36,40,39,41,45,44,46,50,49,51,55,54,56,60,61,62),unpack=True) ## Empirical direct SFH, with dust, Noll
zd,sfr1d,sfr1dl,sfr1du,sfr2d,sfr2dl,sfr2du,sfr3d,sfr3dl,sfr3du,sfr4d,sfr4dl,sfr4du,ebvd,ebvdl,ebvdu,logMd,logMdl,logMdu = np.loadtxt("output/%s/%s.dat"%(folderd,folderd),dtype=float,skiprows=2,usecols=(2,5,4,6,10,9,11,15,14,16,20,19,21,25,24,26,30,29,31),unpack=True) ## Empirical direct SFH, no dust, Calzetti
#zd,sfrd,sfrdl,sfrdu,aged,agedl,agedu,ebvd,ebvdl,ebvdu,logz,logzl,logzu,logMd,logMdl,logMdu = np.loadtxt("output/%s/%s.dat"%(folderd,folderd),dtype=float,skiprows=2,usecols=(2,5,4,6,10,9,11,15,14,16,20,19,21,25,24,26),unpack=True) ## Constant SFR, metallicity free, Calzetti

ind = []
j=0
for i in range(len(fieldnd)):
    if j==len(fieldd): break
    if fieldnd[i]==fieldd[j] and photidnd[i]==photidd[j]:
        ind.append(i)
        #print fieldnd[i],fieldd[j],photidnd[i],photidd[j]
        j+=1
assert len(ind)==len(fieldd)
assert (fieldnd[ind]==fieldd).all()
assert (photidnd[ind]==photidd).all()
assert (znd[ind]==zd).all()

#plotting(sfrnd[ind],sfrndl[ind],sfrndu[ind],sfrd,sfrdl,sfrdu,zd,"SFR",foldernd)
plotting(sfr1nd[ind],sfr1ndl[ind],sfr1ndu[ind],sfr1d,sfr1dl,sfr1du,zd,"SFR1",foldernd)
plotting(sfr2nd[ind],sfr2ndl[ind],sfr2ndu[ind],sfr2d,sfr2dl,sfr2du,zd,"SFR2",foldernd)
plotting(sfr3nd[ind],sfr3ndl[ind],sfr3ndu[ind],sfr3d,sfr3dl,sfr3du,zd,"SFR3",foldernd)
plotting(sfr4nd[ind],sfr4ndl[ind],sfr4ndu[ind],sfr4d,sfr4dl,sfr4du,zd,"SFR4",foldernd)
#plotting(sfr5nd[ind],sfr5ndl[ind],sfr5ndu[ind],sfr5d,sfr5dl,sfr5du,zd,"SFR5",foldernd)
#plotting(agend[ind],agendl[ind],agendu[ind],aged,agedl,agedu,zd,"Age",foldernd)
plotting(ebvnd[ind],ebvndl[ind],ebvndu[ind],ebvd,ebvdl,ebvdu,zd,"E(B-V)",foldernd)
#plotting(deltand[ind],deltandl[ind],deltandu[ind],deltad,deltadl,deltadu,zd,"delta",foldernd)
#plotting(Ebnd[ind],Ebndl[ind],Ebndu[ind],Ebd,Ebdl,Ebdu,zd,"E_b",foldernd)
plotting(logMnd[ind],logMndl[ind],logMndu[ind],logMd,logMdl,logMdu,zd,"log(M)",foldernd)

#plotting2(umin,uminl,uminu,2.0,gamma,gammal,gammau,0.05,zd,"Umin","Gamma",foldernd)
#plotting2(umin,uminl,uminu,2.0,qpah,qpahl,qpahu,2.5,zd,"Umin","Q_PAH",foldernd)
#plotting2(gamma,gammal,gammau,0.05,qpah,qpahl,qpahu,2.5,zd,"Gamma","Q_PAH",foldernd)

#ploterr(sfrnd[ind],sfrndl[ind],sfrndu[ind],sfrd,sfrdl,sfrdu,"SFR",foldernd,folderd)
#ploterr(agend[ind],agendl[ind],agendu[ind],aged,agedl,agedu,"Age",foldernd)
ploterr(sfr1nd[ind],sfr1ndl[ind],sfr1ndu[ind],sfr1d,sfr1dl,sfr1du,"SFR1",foldernd,folderd)
ploterr(sfr2nd[ind],sfr2ndl[ind],sfr2ndu[ind],sfr2d,sfr2dl,sfr2du,"SFR2",foldernd,folderd)
ploterr(sfr3nd[ind],sfr3ndl[ind],sfr3ndu[ind],sfr3d,sfr3dl,sfr3du,"SFR3",foldernd,folderd)
ploterr(sfr4nd[ind],sfr4ndl[ind],sfr4ndu[ind],sfr4d,sfr4dl,sfr4du,"SFR4",foldernd,folderd)
#ploterr(sfr5nd[ind],sfr5ndl[ind],sfr5ndu[ind],sfr5d,sfr5dl,sfr5du,"SFR5",foldernd)
ploterr(ebvnd[ind],ebvndl[ind],ebvndu[ind],ebvd,ebvdl,ebvdu,"E(B-V)",foldernd,folderd)
#ploterr(deltand[ind],deltandl[ind],deltandu[ind],deltad,deltadl,deltadu,"delta",foldernd,folderd)
#ploterr(Ebnd[ind],Ebndl[ind],Ebndu[ind],Ebd,Ebdl,Ebdu,"E_b",foldernd,folderd)
ploterr(logMnd[ind],logMndl[ind],logMndu[ind],logMd,logMdl,logMdu,"log(M)",foldernd,folderd)

# plt.figure()
# plt.hist(logz,bins=20)
# plt.xlabel("Log(Z)")
# plt.savefig("logZHist.png")

# plt.figure()
# #plt.plot(zd,logz,'bo',label='')
# plt.errorbar(zd,logz,yerr=[logz-logzl,logzu-logz],fmt='bo',label='')
# plt.hlines(-2.1135,min(zd),max(zd),colors='r',label=r'$%s=%.3f$'%('Z',0.0077))
# #plt.errorbar([min(zd)],[max(logz)],yerr=[(np.average(logz-logzl),np.average(logzu-logz))],fmt='',label='Typical Error')
# plt.xlabel("Redshift")
# plt.ylabel("log(Z)")
# plt.savefig("logZvsz.png")