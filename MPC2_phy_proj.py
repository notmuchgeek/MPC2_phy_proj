'''
List of all the fucntions:



'''

import numpy as np #import the numpy package since this will be useful
import matplotlib.pyplot as plt #plotting package


N=5
e=1.602e-19 #in Coulombs
eps0=8.85e-12

K=1/(4*np.pi*eps0)




# define base functions: dist and unit_v
def dist(r0,r1):
    if len(r0)==3:
        return np.sqrt(abs((r0[0]-r1[0])**2+(r0[1]-r1[1])**2+(r0[2]-r1[2])**2))
    elif len(r0)==2:
        return np.sqrt(abs((r0[0]-r1[0])**2+(r0[1]-r1[1])**2))
    else:
        print('Dimension error')

def unit_v(r0,r1):
    return (r1-r0) / np.linalg.norm((r1-r0))




# define potential, electric field and the summation version 
def potl(r0,r1,q1):
    V=K * (q1*e/dist(r0,r1))
    return V

def Efield(r0,r1,q1):
    #r0 and r1 should be numpy arrays with 3 elements
    E=K* (q1*e/dist(r0,r1)**2)*unit_v(r0,r1)
    #add calculation of E here
    return E

def potl_sum(r0,ri,qi):
    V_sum = 0
    for i in range(len(ri)):
        V = potl(r0,ri[i],qi[i])
        V_sum += V
    return V_sum #return the total potential

def Efield_sum(r0,ri,qi):
    E_sum = 0
    for i in range(len(ri)):
        E = Efield(r0,ri[i],qi[i])
        E_sum += E
    return E_sum #return the total potential




# define total potential energy
def potl_energy_sum(ri,qi):
    U=0 #initialise the total potential energy
    chargeadded=[] #list of indices of charges already added
    for loop in range(len(qi)): #loop over each charge in turn
        #....incomplete code below which you can use as a starting point
        for j in chargeadded: #loop over charges already added (bringing charge loop towards charge j)            
            Uij=e*qi[j] * potl(ri[j],ri[loop],qi[loop])
            #print('Adding PE %.2e [units] of bringing charge %d towards charge %d' % (Uij,loop,j))
            U=U+Uij
        chargeadded.append(loop) #add the index of the added charge to the list
    #print('Total PE is %.2e [units]' % U)
    return U #return the total potential energy



Natomslist=[3,3,1] #list numbers of atoms in each type of molecule#.  type 0 is CO2, type 1 is water
totalatomtypes=sum(Natomslist)
atomindices=[[0,1,2],[3,4,5],[6,]] #give each distinct atom in each molecule a different index. The three atoms in type 0 are numbered 0,1,2 etc
#this will be useful so we can plot the different atoms with different colours
cols='krrryyb'




# define atom postion and transformation operations
def atom_base_positions(type):
    #return the coordinates of each atom of the molecule of this type
    #in spherical polar coordinates, relative to an 'anchor' at the origin
    if type==0: #definitions for CO2
        d=160e-3 #C-O length in nm
        p=0.1*3.3356e-30 #C-O dipole moment magnitude in SI units Cm (using value for CO; C-O bond dipole in CO2 hard to measure as the two dipoles cancel out and the molecule has no net dipole)
        qeff=(p/(d*1e-9))/1.602e-19 #effective charge in units of e
        Na=3
        r,phi,theta,q=np.empty(Na),np.empty(Na),np.empty(Na),np.empty(Na)
        r[0],phi[0],theta[0],q[0]=0.0,0.0  ,0.0   ,+2*qeff #place the C atom at the origin
        r[1],phi[1],theta[1],q[1]=d  ,0.0  ,0.0   ,-qeff #O atom 1, a distance d from the origin at angle theta=0, phi=0 (so along the +z axis)
        r[2],phi[2],theta[2],q[2]=d  ,0.0  ,np.pi ,-qeff #O atom 2, a distance d from the origin at angle theta=pi, phi=0 (so along the -z axis)            
        
    elif type==1: #definitions for water (incomplete)
        d=0.09578 #O-H bond length in nm
        #d=0.09745401671213969 #O-H bond length in nm
        p=1.85*3.3356e-30 #O-H dipole moment magnitude in SI units Cm 
        qeff=(p/(d*1e-9))/1.602e-19 #effective charge in units of e
        Na=3
        r,phi,theta,q=np.empty(Na),np.empty(Na),np.empty(Na),np.empty(Na)
        r[0],phi[0],theta[0],q[0]=0.0,0.0  ,0.0   ,-2*qeff #place the O atom at the origin
        r[1],phi[1],theta[1],q[1]=d  ,0.0  ,151/720*np.pi   ,+qeff 
        r[2],phi[2],theta[2],q[2]=d  ,0.0  ,((1-151/720)*np.pi) ,+qeff

    elif type==2: #definitions for Na+
        Na=1
        r,phi,theta,q=np.empty(Na),np.empty(Na),np.empty(Na),np.empty(Na)
        r[0],phi[0],theta[0],q[0]=0.0,0.0 ,0.0   ,+1 #place the Na atom at the origin
    
    else:
        r,phi,theta,q=np.array([]),np.array([]),np.array([]),np.array([]) #if type is not set then return empty lists

    return r,phi,theta,q

def rotmol_atomposns(moltype,rmtranslation,phimrot,thetamrot):
    #get base coordinates of atoms in molecule (relative to anchor, prior to rotation around anchor point)
    rb,phib,thetab,qb=atom_base_positions(moltype) #type 1 for water

    outr=[] #list of all vector positions of atoms in molecule
    outq=[] #list of atom charges
    for batom in range(len(rb)): #cycle through atoms in molecule and generate rotated coordinates
        rcarts=rb[batom]*np.array([np.sin(thetab[batom]-thetamrot)*np.cos(phib[batom]+phimrot),
                         np.sin(thetab[batom]-thetamrot)*np.sin(phib[batom]+phimrot),
                         np.cos(thetab[batom]-thetamrot)])
        rvector=rmtranslation+rcarts
        outr.append(rvector)
        outq.append(qb[batom])
    return outr,qb

# deifne plot
def mol_plot(moltypes,xa,ya,za,phia,thetaa): #plot out molecule positions and create associated plots. returns total potential
    #create the figure
    f2=plt.figure()
    ax2=plt.subplot(1,1,1,projection='3d')

#    ploth=[] #create list to store plot handles
    
    #initialise one x array, one y array and one z array for each type of atom
    xplot=[]
    yplot=[]
    zplot=[]
    for loop in range(totalatomtypes):
        xplot.append([])
        yplot.append([])
        zplot.append([]) 
    
    #the Python syntax in the next line using zip() allows you to cycle through several arrays (of the same size) at once
    for x,y,z,moltype,phi,theta in zip(xa,ya,za,moltypes,phia,thetaa): #loop through molecules
        rpositions,qvalues=rotmol_atomposns(moltype,[x,y,z],phi,theta) #generate arrays of all atom positions
        
        #the Python syntax in the next line using enumerate() gives you a counter starting at zero as well as cycling through the elements in an array like a normal Python loop does
        for counter, atomindex in enumerate(atomindices[int(moltype)]): #loop through the atoms in this molecule
            rthisatom=rpositions[counter]
            #add the atom coordinates to the relevant dataseries for the plot
            xplot[atomindex].append(rthisatom[0])
            yplot[atomindex].append(rthisatom[1])
            zplot[atomindex].append(rthisatom[2])
   
    #now draw all the atom dataseries onto the plot
    for atomindex in range(totalatomtypes): #loop over all the types of atom
        #ploth.append(ax2.scatter3D(xplot[atomindex],yplot[atomindex],zplot[atomindex],c=cols[atomindex])) #plot atom and store plot handle
        #print(atomindex)
        #print(xplot[atomindex])
        ax2.scatter3D(xplot[atomindex],yplot[atomindex],zplot[atomindex],c=cols[atomindex]) #plot atom and store plot handle
    
    #add lines to represent the bonds between atoms in the molecule
    #plot these as just one data series
    #the code below interleaves NaN values to ensure no line between different molecules
    padar=np.empty(len(xplot[0]))
    padar[:]=np.NaN
    #interleave values into a 1D array
    listx = [xplot[0], xplot[1], padar, xplot[0], xplot[2], padar]
    xi=[val for tup in zip(*listx) for val in tup] #look at the results of this line if you want to work out what it is
    listy = [yplot[0], yplot[1], padar, yplot[0], yplot[2], padar]
    yi=[val for tup in zip(*listy) for val in tup]
    listz = [zplot[0], zplot[1], padar, zplot[0], zplot[2], padar]
    zi=[val for tup in zip(*listz) for val in tup]
    ax2.plot3D(xi,yi,zi,'k') #add the lines joining the atoms of type 0

    #repeat for type 1
    #this assume that there are bonds between the 1st and 2nd atoms in the coordinates list
    #and between the 1st and 3rd
    #(edit the indices in the below code if you set up your atom_base_posns functions differently)

    padar=np.empty(len(xplot[3]))
    padar[:]=np.NaN
    #interleave values into a 1D array
    listx = [xplot[3], xplot[4], padar, xplot[3], xplot[5], padar]
    xi=[val for tup in zip(*listx) for val in tup]
    listy = [yplot[3], yplot[4], padar, yplot[3], yplot[5], padar]
    yi=[val for tup in zip(*listy) for val in tup]
    listz = [zplot[3], zplot[4], padar, zplot[3], zplot[5], padar]
    zi=[val for tup in zip(*listz) for val in tup]
    ax2.plot3D(xi,yi,zi,'r') #add the lines joining the atoms of type 1

    ##add labels
    #ax2.set_xlabel('x (add unit)')
    #ax2.set_ylabel('y (add unit)')
    #ax2.set_zlabel('z (add unit)')
    #ax2.set_title('molecule positions')
    
    #ax2.set_aspect('auto')
        
    ##ax2.view_init(elev=10., azim=30.) #adjust 'camera angle' with this command if desired - angles are in degrees
    #plt.show()
    return f2,ax2 #returns the figure axis handle. This could be useful if you want to edit the figure outside of the function


def allatomposns(moltypes,xa,ya,za,phia,thetaa): 
    rlist=[]
    qlist=[]
    
    #the Python syntax in the next line using zip() allows you to cycle through several arrays (of the same size) at once
    for x,y,z,moltype,phi,theta in zip(xa,ya,za,moltypes,phia,thetaa): #loop through molecules
        rpositions,qvalues=rotmol_atomposns(moltype,[x,y,z],phi,theta) #generate arrays of all atom positions
        for rp,qv in zip(rpositions,qvalues):
            rlist.append(rp)
            qlist.append(qv)
    return rlist,qlist


def yieldPotlEfield(ri,qi,num_samp,start_coord,end_coord):
    samp_x= np.linspace(start_coord,end_coord,num_samp) #sampling length along z axis
    samp_y= np.linspace(start_coord,end_coord,num_samp)
    samp_z= np.linspace(start_coord,end_coord,num_samp)
    samp_yy, samp_zz, samp_xx = np.meshgrid(samp_y, samp_z,samp_x)

    #calculate potl, Efield for each test points          
    potl_molecules = []
    Efield_molecule2D = []

    for i in range(len(samp_xx)):
        pinlistV1=[]
        pinlistE1_2D=[]
        for j in range(len(samp_xx[i])):
            pinlistV2=[]
            pinlistE2_2D=[]

            for k in range(len(samp_xx[i][j])):
                EfieldVectorL2D=[]
                r0 = np.array([samp_xx[i][j][k],samp_yy[i][j][k],samp_zz[i][j][k]])
                #Calcluate 1D potl for each in 3D
                pinlistV2.append(potl_sum(r0,ri,qi))
                
                #Calcluate 2D Efield (x,z) for each in 3D
                EfieldVector = Efield_sum(r0,ri,qi)
                EfieldVectorL2D.append(EfieldVector[0])
                EfieldVectorL2D.append(EfieldVector[1])
                pinlistE2_2D.append(EfieldVectorL2D)
            pinlistV1.append(pinlistV2)
            pinlistE1_2D.append(pinlistE2_2D)
        potl_molecules.append(pinlistV1)
        Efield_molecule2D.append(pinlistE1_2D)
    return samp_x, samp_y, samp_z, samp_yy, samp_zz, samp_xx, np.array(Efield_molecule2D), np.array(potl_molecules)

# samp_x, samp_y, samp_z, samp_yy, samp_zz, samp_xx, Efield_molecule2D,potl_molecules = yieldPotlEfield(13,-0.5,1.5)
