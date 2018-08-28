import numpy as np
import random as rnd
import matplotlib
matplotlib.use("qt4agg")
import matplotlib.pyplot as plt
def pulseresp(Delta,SaltCon,DidBout):
    SaltChange = np.diff(SaltCon)
    On = SaltChange == Delta
    Off = SaltChange ==-1*Delta
    OnTimes = np.where(On)[0]
    OffTimes = np.where(Off)[0]
    Onno = np.shape(On)
    Offno = np.shape(Off)

    if Onno > Offno:
        np.delete(OnTimes, np.array(Onno.shape) - 1)
    count = 0
    PulseBout = np.zeros_like(OnTimes)
    for s in range(0,np.array(OnTimes.shape)-1):
        pulse = np.arange(OnTimes[s],OffTimes[s])
        bouttimes = np.where(DidBout)[0]
        overlap = np.in1d(pulse,bouttimes-1)
        if np.sum(overlap) > 0:
            PulseBout[count] = 1
        else:
            PulseBout[count] = 0
        count += 1

    X = np.sum(PulseBout)/np.double(PulseBout.shape[0])
    return(X)


alist = np.arange(0.01,0.1,0.005)#habituation rate
along = int(alist.shape[0])
blist = np.arange(0.000001,0.00003,0.000001) # dehabituation
blong = int(blist.shape[0])
clist = np.arange(0.005,0.007,0.0001) # response rate
clong = int(clist.shape[0])
ErrorMatrix = np.zeros((along,blong,clong),dtype='float')
for ai in range(0,along):
    a = alist[ai]
    for bi in range(0,blong):
        b = blist[bi]
        for ci in range(0,clong):
            c = clist[ci]
            TestCons = np.array([0,0])
            TestPunish = np.array([0,10])
            MainMean = np.double(np.zeros_like(TestCons))
            SecondMean = np.double(np.zeros_like(TestCons))
            puncount = 0
            for Pun in TestPunish:
                FishNo = 20
                PulseLength = 10
                ClosedLoop = 1
                MainPunish = Pun
                TurnPunish = MainPunish
                ISImean = 20
                MainSalt = 100
                SecProb = 0
                SecondSalt = 0

                ResponseProbability = np.zeros_like(np.arange(FishNo,dtype='double'))
                ResponseProbability2 = np.zeros_like(np.arange(FishNo,dtype='double'))
                T=1800
                for Fish in range(0,FishNo):
                    SaltCon = np.zeros_like(np.arange(T,dtype='double'))
                    BoutRate = np.zeros_like(SaltCon)
                    Integrator = np.zeros_like(SaltCon)
                    SaltFelt = np.zeros_like(SaltCon)
                    DidBout = np.zeros_like(SaltCon)

                    BoutRate[0] = 0.01
                    SaltOff = 60
                    SaltOn = 0
                    ISI = ISImean

                    for t in range(1,T):
                        SaltFelt[t] = np.maximum(SaltCon[t-1] - Integrator[t-1],-2)
                        Error = np.maximum(Integrator[t-1]-SaltCon[t-1],0)
                        Integrator[t] = Integrator[t-1] + a*SaltFelt[t] - b*Error*Integrator[t-1]
                        BoutRate[t] = 0.01 + c*SaltFelt[t]
                        if SaltCon[t-1] == 0:
                            if (t-SaltOff) > ISI:
                                if rnd.random() < SecProb:
                                    SaltCon[t] = SecondSalt
                                    FillSalt = SecondSalt
                                    TurnPunish = 0
                                else:
                                    SaltCon[t] = MainSalt
                                    FillSalt = MainSalt
                                    TurnPunish = MainPunish
                                SaltOn = t-1
                        else:
                            if (t-SaltOn) > PulseLength:
                                SaltCon[t] = 0
                                SaltOff = t-1
                                ISI =np.random.normal(ISImean,5,1)
                                PulseLength = 10
                            else:
                                SaltCon[t] = FillSalt
                        if rnd.random() < BoutRate[t]:
                            DidBout[t] = 1
                            if SaltCon[t] != 0:
                                PulseLength = ClosedLoop*(t-SaltOn + TurnPunish)
                    ResponseProbability[Fish] = pulseresp(MainSalt,SaltCon,DidBout)
                    if SecondSalt > 0:
                        ResponseProbability2[Fish] = pulseresp(SecondSalt, SaltCon, DidBout)
                    else:
                        ResponseProbability2[Fish] = np.array([0])
                MainMean[puncount] = np.mean(ResponseProbability)
                SecondMean[puncount] = np.mean(ResponseProbability2)
                puncount+=1
            e = (np.float(0.3) - MainMean[1])**2 + (np.float(0.6) - MainMean[0])**2
            ErrorMatrix[int(ai),int(bi),int(ci)] = e
    print[a]
plt.plot(TestPunish,MainMean)
plt.plot(TestPunish,SecondMean)
plt.show()




