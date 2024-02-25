# -*- coding: utf-8 -*-
"""
Created on Sun Feb 24 2024

@author: Sara Linse
"""
#This Python script reads a pdb file and generates a contact plot
#The user can choose the distance cutoff and atom types investigated

#On the following line you can set the longest cutoff distance in Å. 
cutoff=12.0
#Three more cutoffs, 2, 4 and 6 Å shorter, repectively are set
cutoff2=cutoff-2
cutoff4=cutoff-4
cutoff6=cutoff-6
#On the following line you can set the atom type to be investigated
atomtype='CA'
#On the first line where you see f.open you can change the input pdb filename
#Download this script and place the pdb file in the same folder as the script
#On the second line where you see f.open you can change the output text filename
#If you right cick on the plot after it is made you can save it as a png
import matplotlib.pyplot as plt
import math

#Initialize matrices for data
#If you have more than 10000 atoms or 8 chains, change the values inside range()

#initialize matrix for residuenr, atomnr and chainnr
matrix = []
matrix = [[column for column in range(3)] for row in range(10000)]
#asign some integer values
matrix[0][0]=1 #residue nr
matrix[0][1]=1 #atom nr
matrix[0][2]=1 #chain nr

#initialize resmatrix for storing per chain the number of residues 
resmatrix = []
resmatrix = [[column for column in range(2)] for row in range(8)]
#asign some integer values
resmatrix[1][0]=1 

#initialize cordmatrix for x,y,z coordinates
cordmatrix = []
cordmatrix = [[column for column in range(3)] for row in range(10000)]
#asign some float values
cordmatrix[1][0]=0.1
cordmatrix[1][1]=0.1
cordmatrix[1][2]=0.1

#initialize textmatrix for chainid, residuetype and atomtype
textmatrix = []
textmatrix = [[column for column in range(3)] for row in range(10000)]
#asign some string values
textmatrix[1][0]='A'
textmatrix[1][1]='A'
textmatrix[1][2]='A'

#open pdb file of choice and read all lines
f = open("5zh6.pdb","r")
data = f.readlines()
f.close

#initiate some variables to define their data type
x1=0.1 #x-coordinate of atom1
y1=0.1 #y-coordinate of atom1
z1=0.1 #z-coordinate of atom1
x2=0.1 #x-coordinate of atom2
y2=0.1 #y-coordinate of atom2
z2=0.1 #z-coordinate of atom2
residuenr=1
atomnr=1
I=0 #number of atoms
rc=0 #number of chains
rr=0 #number of residues total
chainid='A'
chain='A'
chainnr=0
atom='A'
attype='CX'
dist=15.0
k=0
m=0
l=0
n=0
myline = []

#Read data, ensure correct data type and sort into the matrices
#if this fails you have a varaint pdb file and need to modify the line ending with split() to match the line format

i=0
#initiate for residue counting
p=1
while p<8:
    resmatrix[p][0]=0
    p=p+1
for line in data:
    myline=line.strip()
    char1 = myline[0]
    char2 = myline[1]
    char3 = myline[2]
    #only sort data from lines that start with ATOM
    if ((char1 == 'A') and (char2 == 'T')and (char3 == 'O')):
        atom,atnr,attype,restype,chainid,resnr,x,y,z,nr1,nr2,at=line.split()
        #print("line=",x,y,z)
        x1=float(x)
        y1=float(y)
        z1=float(z)
        residuenr=int(resnr)
        atomnr=int(atnr)
        i=i+1
        #fill textmatrix with chainid, residuetype and atomtype
        textmatrix[i][0]=chainid
        textmatrix[i][1]=restype
        textmatrix[i][2]=attype
        #fill matrix with residuenr, atomnr and chainnr
        #if more than 8 chains add more if statements below
        matrix[i][0]=residuenr
        matrix[i][1]=atomnr
        #check which chain we are in, convert letters to numbers
        #and store number of residues per chain and last atom        
        if textmatrix[i][0]=='A':
            matrix[i][2]=1 #chain1
            if matrix[i][0]!=matrix[i-1][0]:
                resmatrix[1][0]=resmatrix[1][0]+1 #count residues in chain 1       
            resmatrix[1][1]=i #find last atom in chain 1
            rc=1
        if textmatrix[i][0]=='B':
            matrix[i][2]=2 #chain2    
            if matrix[i][0]!=matrix[i-1][0]:
                resmatrix[2][0]=resmatrix[2][0]+1 #count residues in chain 2       
            resmatrix[2][1]=i #find last atom in chain 2
            rc=2       
        if textmatrix[i][0]=='C':
            matrix[i][2]=3 #chain3
            if matrix[i][0]!=matrix[i-1][0]:
                resmatrix[3][0]=resmatrix[3][0]+1 #count residues in chain 3      
            resmatrix[3][1]=i #find last atom in chain 3
            rc=3
        if textmatrix[i][0]=='D':
            matrix[i][2]=4 #chain4
            if matrix[i][0]!=matrix[i-1][0]:
                resmatrix[4][0]=resmatrix[4][0]+1 #count residues in chain 4       
            resmatrix[4][1]=i #find last atom in chain 4
            rc=4
        if textmatrix[i][0]=='E':
            matrix[i][2]=5 #chain5
            if matrix[i][0]!=matrix[i-1][0]:
                resmatrix[5][0]=resmatrix[5][0]+1 #count residues in chain 5        
            resmatrix[5][1]=i #find last atom in chain 5
            rc=5
        if textmatrix[i][0]=='F':
            matrix[i][2]=6 #chain6
            if matrix[i][0]!=matrix[i-1][0]:
                resmatrix[6][0]=resmatrix[6][0]+1 #count residues in chain 6       
            resmatrix[6][1]=i #find last atom in chain 6
            rc=6
        if textmatrix[i][0]=='G':
            matrix[i][2]=7 #chain7
            if matrix[i][0]!=matrix[i-1][0]:
                resmatrix[7][0]=resmatrix[7][0]+1 #count residues in chain 7      
            resmatrix[7][1]=i #find last atom in chain 7
            rc=7
        if textmatrix[i][0]=='H':
            matrix[i][2]=8 #chain8  
            if matrix[i][0]!=matrix[i-1][0]:
                resmatrix[8][0]=resmatrix[8][0]+1 #count residues in chain 8       
            resmatrix[8][1]=i #find last atom in chain 8
            rc=8
        #fill cordmatrix with x,y,z coordinates
        cordmatrix[i][0]=x1
        cordmatrix[i][1]=y1
        cordmatrix[i][2]=z1
        I=atomnr   
#count total number of residues
n=0
while n<rc:
    n=n+1
    rr=rr+resmatrix[n][0]
print('total number of chains =',rc)
print('total number of residues =',rr)
#produce plot and outputfile with contacs below the user defined cutoff distance          
i=0
j=0
# depict illustration
size=(20+rr)*5/100
plt.rcParams["figure.figsize"] = (size, size)
fig = plt.figure()
ax = fig.add_subplot()
#open a file for saving found contacts
f=open("outputfilename.txt","a+")
f.write('res1')
f.write("\t")
f.write('res2')
f.write("\t")
f.write('distance')
f.write("\t")
f.write("\n")
while (i < I):
    i=i+1
    chainnr=matrix[i][2]
    if textmatrix[i][2]==atomtype:
        x1=cordmatrix[i][0]
        y1=cordmatrix[i][1]
        z1=cordmatrix[i][2]
        resnr=matrix[i][0] 
        print('now processing residue nr',resnr)
        k=matrix[i][2] #chain number        
        j=0
        while (j < I):
            j=j+1
            #print("i,j=",i,j)
            if textmatrix[j][2]==atomtype:
                x2=cordmatrix[j][0]
                y2=cordmatrix[j][1]
                z2=cordmatrix[j][2]
                l=matrix[j][2] #chain number 2
                dist=math.sqrt((x1-x2)**2+(y1-y2)**2+(z1-z2)**2)
                #Check if distance is below the cutoff
                if dist<cutoff:
                    if k==1:
                        ai=matrix[i][0] #residue1 nr
                    if k!=1:
                        m=0
                        while (m<k):
                            m=m+1
                            ai=matrix[i][0]+resmatrix[m][0] #residue1 nr
                    if l==1:
                        di=matrix[j][0] #residue2 nr
                    if l!=1:
                        m=0
                        while (m<l):
                            m=m+1
                            di=matrix[j][0]+resmatrix[m][0] #residue2 nr                       
                    d=str(di) #residue2nr as string
                    a=str(ai) #residue1nr as string
                    dstr=str(dist)
                    #if not the same residue write to file
                    if matrix[i][0]!=matrix[j][0]: 
                        f.write(a)
                        f.write("\t")
                        f.write(d)
                        f.write("\t")
                        f.write(dstr)
                        f.write("\t")
                        f.write("\n")
                    ai2=di
                    di2=-ai
                    di1=-di
                    #plot contacts within four cutoffs
                    #you can change the colors below as you wish
                    if dist<cutoff:
                        plt.plot(ai2, di2, 's', markersize=5, color='peachpuff')
                        plt.plot(ai, di1, 's', markersize=5, color='peachpuff')                
                    if dist<cutoff2:
                        plt.plot(ai2, di2, 's',  markersize=5, color='salmon')
                        plt.plot(ai, di1, 's', markersize=5, color='salmon')
                    if dist<cutoff4:
                        plt.plot(ai2, di2, 's', markersize=5, color='red')
                        plt.plot(ai, di1, 's', markersize=5, color='red')
                    if dist<cutoff6:
                        plt.plot(ai2, di2, 's', markersize=5, color='darkred')
                        plt.plot(ai, di1, 's', markersize=5, color='darkred')
                    plt.xlabel('residue 1')
                    plt.ylabel('-residue 2')
                    plt.title('Contact plot')
# make plot in square format
ax.set_aspect(1.0/ax.get_data_ratio(), adjustable='box')
plt.show()
f.close 
   