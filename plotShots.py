def plot_shot_gather_byRange(obsnum,obsfile,shot_beg,shot_end,tlims,throw,ktpfile,linefile,trscale,vred,savepdf,fsize):
    '''
    VJS 9/2014
    Modified 3/2016
    plot_shot_gather.py
    Plot shot gather for one OBS, including picks, by along-track range (not by shot)
    Input:
           obsnum:         	    OBS number for shot gather
           obsfile:                String with full path for obs file
           shot_beg, shot_end:     First and last shot to plot
           tlims:                  Array: [t0 tend]
           throw:                  if ==1, throw out landshots. if == 0, keep land shots.
           ktpfile:                String with full path for ktp file
           linefile:               String with full path for shot file (with shot, x, y, z: along-track.)
           trscale:                Number by which to scale/normalize amplitudes
           vred:                   Reduction velocity, m/s
           savepdf:                If ==1, save pdf form.  If ==0, do not save pdf.
           fsize:                  Figure size, like: fsize=(10,5)
    '''
    
    import sys
    import numpy as np
    import obspy as obs
    import matplotlib.pyplot as plt
    from matplotlib import rc
    from struct import unpack
    sys.path.append('/Users/vjsahakian/Software/PY')
    import tred
    import read_ktp
    
    #Define fonts for plot (default for axes, etc):
    rc('font',**{'family':'sans-serif','sans-serif':['Helvetica'],'size':11})
    
    #Read OBS file
    obsn=obs.read(obsfile)
    
    #Enter time limits of plot:
    t0=tlims[0]
    tend=tlims[1]
    
    #Read in shot info
    shotinfo=np.loadtxt(linefile)
    shotnums=shotinfo[:,0]
    
    #Find shots within range of shot_beg to shot_end
    shotind=np.where((shotnums>=shot_beg) & (shotnums<= shot_end))[0]
    #Pull out x position of shots (shotx)
    shotx=shotinfo[shotind,1]
    
    #shots from file:
    fshots=shotnums[shotind]
    
    #These shots aren't necessarily the same index as those in the 
    #segy file, since these shots include ALL shots for line 7 (including land shots)
    segyshots=np.zeros((len(obsn),1))
    for i in range(len(obsn)):
        header=obsn[i].stats.segy.trace_header.unpacked_header
        segyshots[i]=unpack(">i",header[8:12])[0]
    
    #What indices in the segy file correspond to the shots in the requested range?    
    #this index may also be used for plotting, to refer to the amplitudes
    segshotind=np.where(fshots==segyshots)[0]

    #Get source to receiver range
    rng=np.zeros(np.size(segshotind))
    for k in range(len(segshotind)):
        current_shot=segshotind[k]
        header=obsn[current_shot].stats.segy.trace_header.unpacked_header
        rng[k]=unpack(">i",header[36:40])[0]
    shot_rang_match=[fshots,rng]
    
    #Convert ktp file to lists/arrays with the various pick segments for each branch
    [pg, pn, pb]=read_ktp.read_ktp(ktpfile,throw,shot_end)
    pg,pn,pb=np.array(pg),np.array(pn),np.array(pb)
                    
    #Find source-receiver ranges for each in order to get time with reduction velocity
    #Range arrays:
    pgin,pnin,pbin=[],[],[]
    pgout,pnout,pbout=[],[],[]
    pgr,pnr,pbr=[],[],[]
    
    #For the times and shots:
    ktpin=[]
    sin=[]
    xin=[]
    
    #Pick times arrays:
    ktpgout,ktpnout,ktpbout=[],[],[]
    ktpg,ktpn,ktpb=[],[],[]
    
    #along-line distance arrays:
    xgout,xnout,xbout=[],[],[]
    xg,xn,xb=[],[],[]
    
    #Shot number arrays:
    sgout,snout,sbout=[],[],[]
    sg,sn,sb=[],[],[]
    
    #For Pg
    #for each set of picks in the ktp file:
    for i in range(pg.shape[1]):
        #zero the tmp arrays:
        pgin,sin,ktpin,xin=[],[],[],[]
        #for each shot in each set of the pick file:
        for j in range(len(pg[0][i])):
            #Get the shot number from the pick data:
            snum=pg[0][i][j]
            #Find the index in the shot range match array that corresponds to this shot:
            ind=np.where(shot_rang_match[0]==snum)[0]
            #If such a correlation exists...
            if ind:
                #append the source-receiver range to a tmp:
                pgin.append(shot_rang_match[1][ind][0])
                #also the shot number:
                sin.append(pg[0][i][j])
                #The pick time:
                ktpin.append(pg[1][i][j])
                #And the along-track distance:
                xin.append(shotx[ind][0])
        #Now cat these tmps into an list of lists (sets):
        pgout.append(pgin)
        sgout.append(sin)
        ktpgout.append(ktpin)
        sgout.append(xin)            
    #And turn them into arrays:
    pgr=np.array(pgout)
    sg=np.array(sgout)
    ktpg=np.array(ktpgout)
    xg=np.array(xgout)
     
    #Repeat...       
    #For Pn
    for i in range(pn.shape[1]):
        pnin,sin,ktpin,xin=[],[],[],[]
        for j in range(len(pn[0][i])):
            snum=pn[0][i][j]
            ind=np.where(shot_rang_match[0]==snum)[0]
            if ind:
                pnin.append(shot_rang_match[1][ind][0])
                sin.append(pn[0][i][j])
                ktpin.append(pn[1][i][j])
                xin.append(shotx[ind][0])
        pnout.append(pnin)
        snout.append(sin)
        ktpnout.append(ktpin)
        xnout.append(xin)
    pnr=np.array(pnout)  
    sn=np.array(snout)
    ktpn=np.array(ktpnout)    
    xn=np.array(xnout)
            
    #For Pb
    for i in range(pb.shape[1]):
        pbin,sin,ktpin,xin=[],[],[],[]
        for j in range(len(pb[0][i])):
            snum=pb[0][i][j]
            ind=np.where(shot_rang_match[0]==snum)[0]
            if ind:
                pbin.append(shot_rang_match[1][ind][0])
                sin.append(pb[0][i][j])
                ktpin.append(pb[1][i][j])
                xin.append(shotx[ind][0])
        pbout.append(pbin)
        sbout.append(sin)
        ktpbout.append(ktpin)
        xbout.append(xin)
    pbr=np.array(pbout)  
    sb=np.array(sbout)
    ktpb=np.array(ktpbout)
    xb=np.array(xbout)  
    
    
    #Plot!
    print 'Plotting OBS '+str(obsnum)+' gather'
    plt.close('all')
    f=plt.figure(num=1, figsize=fsize)
    
    for i in range(len(fshots)):
        #Get the times from the segy:
        t=obsn[segshotind[i]].times()
        #Reduce them:
        tnew=tred.tred(t,rng[i],vred)
        #Get the amplitudes:
        amp=obsn[segshotind[i]].data
        #Scale by 0.18 - arbitrary, based on visual.  This is to keep it normalized between the average shot range difference, ~0.08 for obs72
        amp=trscale*(amp/max(abs(amp)))
        #Get positive part for shading
        ineg=amp<=0
        amp_pos=amp.copy()
        amp_pos[ineg]=0
        #Get negative part for shading
        ipos=amp>=0
        amp_neg=amp.copy()
        amp_neg[ipos]=0
        tind=np.where((tnew>=t0) & (tnew<=tend))
        plt.fill_betweenx(tnew[tind],shotx[i],amp_pos[tind]+shotx[i],facecolor='gray') #Shade positive
        #plt.fill_betweenx(tnew,amp_neg+i,i,facecolor='red') #Shade negative
        plt.plot(amp[tind]+shotx[i],tnew[tind],color='0.3') #Plot unshaded wiggle
        
        
    #Plot travel-time picks - number of segments is in the 1 index of *.shape
    #Now, plot picks...
    
    #Pg:
    #Loop over the number of sets in the pick file for this pick type:
    if len(xg)>0:
        for picki in range(pgr.shape[0]):
            #Sort these based on their along-track distance:
            dsort=np.argsort(xg[picki])
            #Make a new array to plot for each set, zero it out now before assigning
            ttmp=np.array([])
            #Loop over each element in the set to get the reduced time for the pick
            for disti in range(len(dsort)):
                ttmp=np.r_[ttmp,tred.tred(ktpg[picki][disti],pgr[picki][disti],vred)]
            #Now plot each set of picks in the file:
            plt.plot(xg[picki],ttmp,color='k',alpha=0.7,linewidth=2)            
    
    #Repeat for other pick types..
    
    #Pn:
    if len(xn)>0:
        for picki in range(pnr.shape[0]):
            dsort=np.argsort(xn[picki])
            ttmp=np.array([])
            for disti in range(len(dsort)):
                ttmp=np.r_[ttmp,tred.tred(ktpn[picki][disti],pnr[picki][disti],vred)]
            #Plot
            plt.plot(xn[picki],ttmp,color='m',alpha=0.7,linewidth=2)
    
    #Pb:
    if len(xb>0):
        for picki in range(pbr.shape[0]):
            dsort=np.argsort(xb[picki])
            ttmp=np.array([])
            for disti in range(len(dsort)):
                ttmp=np.r_[ttmp,tred.tred(ktpb[picki][disti],pbr[picki][disti],vred)]
            #Plot
            plt.plot(xb[picki],ttmp,color='m',alpha=0.7,linewidth=2)
    
    #plt.ylim([0,1.7])
    #plt.show()
    #plt.plot(amp+i,tnew,color='k')
    plt.ylim([t0,tend])
    #plt.xlim([min(obsn[0].data),len(obsn)+max(obsn[len(obsn)-1].data)])
    plt.xlim([min(shotx)-2*trscale, max(shotx)+2*trscale])
    plt.xlabel('Along-track Distance (km)',fontsize=12, fontname='Helvetica')
    plt.ylabel('Time (s) - X/'+str(vred/1000)+'km/s',fontsize=12, fontname='Helvetica')
    plt.title('OBS '+str(obsnum)+' Shot Gather',fontsize=12, fontname='Helvetica')
    
    #Save figure
    plt.savefig('obs'+str(obsnum)+'.png')
    
    if savepdf==1:
        plt.savefig('obs'+str(obsnum)+'.svg',format='svg')
        plt.show()
    else:
        plt.show()
    
    