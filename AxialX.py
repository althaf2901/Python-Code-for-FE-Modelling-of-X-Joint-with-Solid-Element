import pandas as pd
import numpy as np
import os
print('Enter the drive')
drive = input()
print('Enter the folder name:')
folder= input()
os.mkdir(''''''+str(drive.upper()) +''':\\'''+str(folder)+'\\')
df = pd.read_excel('X_axial_testing_90.xlsx')
column=['filename','alpha','beta','gamma','tau','C_D','e_f_c','e_f_b','num_B_T','num_C_T']
df=df.reindex(column, axis='columns')
df.reset_index(inplace = True)
#Materials Input with Dimensions
import os
import math
E = 210000 #Modulus of elasticity
v = 0.3 #Poisson's ratio
for count in range(len(df)):
    inp_dir =count
    path = ''''''+str(drive.upper()) +''':\\'''+str(folder)+'\\'+str(inp_dir)+'/'
    directory = ''''''+str(drive.upper()) +''':\\'''+str(folder)+'\\'
    os.mkdir(path)
    alpha = df['alpha'][count]
    beta =df['beta'][count]
    gamma = df['gamma'][count]
    tau = df['tau'][count]
    C_D = df['C_D'][count]
    e_f_c = df['e_f_c'][count]
    e_f_b = df['e_f_b'][count]
    num_B_T = df['num_B_T'][count]
    num_C_T = df['num_C_T'][count]
    C_L = 0.5*alpha*C_D
    B_D = C_D*beta
    C_T = 0.5*C_D/gamma
    B_T = tau*C_T
    a_r_c = (1.64-(beta-0.2)*(1.64-0.75)*2)
    a_r_b = (0.8+(beta-0.2)*(1.2-0.8)*2)
    a_r_cc = 1/beta
    a_r_cs = (beta+(1.1-(beta-0.2)*0.3)*0.5*(1+beta))/(0.5*beta*3.14) 
    if alpha<= 12:
        ratio_br = 1
        ratio_ch = 1
    else:
        ratio_br = alpha/8
        ratio_ch = alpha/8 
    C_M_L =C_L/2
    C_R = C_D/2
    C_I_R = C_R - C_T
    B_R = B_D/2
    B_I_R = B_R - B_T
    B_L = round(C_M_L) 
    z = (1.1-(beta-0.2)*0.3)*0.5*(B_R+C_R)
    bcl = min(e_f_b*C_R,0.95*(B_L-C_R))
    z_1 = (e_f_c-z/C_R)*C_R
    theta_cut = 90-math.atan(C_R/(B_R+z-0.5*B_R))*180/math.pi
    ao = 0.2*math.sqrt(B_R*B_T)
    lrmin = max(4,0.4*C_T )
    if ao<=lrmin:
        a=lrmin
    else:
        a=ao
    lrmin = max(4,0.4*B_T )
    if ao<=lrmin:
        a_b=lrmin
    else:
        a_b=ao
    b_b = 0.65*math.sqrt(B_R*B_T)
    b_cc = 0.4*math.sqrt(math.sqrt(B_R*C_R*B_T*C_T))
    b_cs = math.pi*C_R/36
    e_s_c = C_T/num_C_T
    e_s_b = B_T/num_B_T
    num_cir_crown = 1+round(0.5*math.pi*C_R/(e_s_c*a_r_cc))
    num_cir_saddle = 1+round((B_R+z)/(e_s_c*a_r_cs))
    c= 1+round(a_r_c*(b_cc-a)/e_s_c)
    d= 1+round(a_r_b*(b_b-a_b)/e_s_b)
    num_weld_c = 1+round(a_r_c*a/(e_s_c))
    num_weld_b = 1+round(a_r_b*a_b/e_s_b)
    num_br_cut = round((2*d+num_weld_b)+a_r_b*(bcl+a_b-2*b_b)/e_s_b)
    num_ch_cut = round((2*c+num_weld_c)+a_r_c*((z+a-2*b_cs)/e_s_c))
    num_ch_cut_1 = round((a_r_c*((z_1)/e_s_c)))
    e_s_b_cut = (bcl+a-2*b_cs)/(num_br_cut-(num_weld_b+2*d)) 
    e_s_c_cut =  (z+a-2*b_cs)/(num_ch_cut-(num_weld_c+2*c))
    num_brace = round((B_L-C_R-bcl)/(ratio_br*e_s_b))
    num_chord = round((C_M_L-B_R-z-z_1)/(ratio_ch*e_s_c))
    e=round(2*B_I_R/e_s_c)
    num_ch_bt = round((1-(beta-0.2))*num_cir_crown)
    crown_r = round(B_R+(0.5*B_T),2)  
    weld_saddle = round((1-(math.asin(B_R/C_R)*180/math.pi)/45)*0.5*B_T ,1)
    weld_sad1 = max(0.05*B_T,weld_saddle)
    saddle_r = round(B_R + weld_sad1,2)
    dihedral_ang = 90+math.asin(B_R/C_R)*180/math.pi  
    jt_ang_cr= 46.5
    jt_ang_sd= min(54,37.5+(0.225*(dihedral_ang-50))) 
    jt_ang_sd= 37.5+(0.225*(dihedral_ang-50)) 
    with open(path +'model.txt','a+') as f:
      f.truncate(0)
      f.seek(0)
      f.write('/PREP7\n')
      f.write('BOPTN,NWARN,1\n') 
      f.write('ET,1,SOLID186\n')
      f.write('MPTEMP,,,,,,,,\n')
      f.write('MPTEMP,1,0 \n')
      f.write('MPDATA,EX,1,,'+str(E)+'\n')
      f.write('MPDATA,PRXY,1,,'+str(v)+'\n')
      f.write('CYL4,0,0,'+ str(C_I_R)+',,' +str(C_R)+',,'+str(C_M_L)+'\n')
      f.write('wpro,,-90.000000,\n')
      f.write('CYL4,0,0,'+ str(B_I_R)+',,' +str(B_R)+',,'+str(B_L)+'\n')
      f.write('wprota,,90 \n')
      f.write('VSBW,  2 \n')
      f.write('VDELE,  3,,,1\n')
      f.write('wprota,,,90 \n')
      f.write('vsbw,all \n')
      f.write('vdele,2,,,1 \n')
      f.write('vdele,6,,,1 \n')
      f.write('vptn,all \n')
      f.write('vdele,2,,,1 \n')
      #Step 2
      p=1001
      teta = [0,10,20,30,40,50,60,70,80,90]
      q=(crown_r-saddle_r)/9
      r=[]
      for i in range(10):
        r.append(saddle_r + (q*i))
      for i in range(10):
          f.write('k,'+str(p+i)+','+str(r[i]*math.cos(teta[i]*math.pi/180))+','+str(0)+','+str(r[i]*math.sin(teta[i]*math.pi/180))+'\n')
      f.write('ksel,s,,,'+str(p)+','+str(p+9)+'\n')
      f.write('bsplin,all'+'\n')
      f.write('kdele,'+str(p+1)+','+str(p+8)+'\n')
      f.write('LGEN,2,4,,,,'+str(C_D)+',,,'+str(0)+'\n')
      f.write('Askin,4,5'+'\n')
      f.write('btol,0.01'+'\n')  
      f.write('vsba,7,6'+'\n')
      f.write('btol,0.00001'+'\n')  
      #Step 3
      q=2001
      r,y,theta_1,theta_2,incl_ang,y_w,y_w1=[],[],[],[],[],[],[]
      grad = (jt_ang_cr-jt_ang_sd)/9
      for i in range(10):
        y.append(C_R*math.cos(math.asin((B_R*math.cos(teta[i]*math.pi/180))/(C_R))))
        theta_1.append(math.asin(B_I_R*math.cos(teta[i]* math.pi/180)/C_R)*180/math.pi)
        incl_ang.append(jt_ang_sd+grad*i)
        theta_2.append(incl_ang[i]-theta_1[i])
        y_w.append(B_T*(math.tan(theta_1[i]*math.pi/180)+math.tan(theta_2[i]*math.pi/180))+y[i])
      sad_h=y_w[0]
      crown_r1=y_w[9]
      saddle_r1=math.sqrt(sad_h*sad_h+B_R*B_R)
      grad1=(crown_r1-saddle_r1)/9
      for i in range(10):
        r.append(saddle_r1+i*grad1)
        y_w1.append(math.sqrt(r[i]*r[i]-(B_R*math.cos(teta[i]*math.pi/180))*(B_R*math.cos(teta[i]*math.pi/180))))
      for i in range(10):
        f.write('k,'+str(q+i)+','+str(B_R*math.cos(teta[i]*math.pi/180))+','+str(y_w[i])+','+str(-10)+'\n')
      f.write('ksel,s,,,'+str(q)+','+str(q+9)+'\n')
      f.write('bsplin,all'+'\n')
      f.write('kdele,'+str(q+1)+','+str(q+8)+'\n')
      f.write('LGEN,2,4,,,,,'+str(B_D)+',,'+str(0)+'\n')
      f.write('Askin,4,5'+'\n')
      f.write('btol,0.01'+'\n')   
      f.write('vsba,4,6'+'\n')
      f.write('btol,0.00001'+'\n')  
      #Step 4 weld simulation
      f.write('allsel,all'+'\n')
      f.write('askin,30,62'+'\n')
      f.write('al,4,28,60'+'\n')
      f.write('al,5,25,58'+'\n')
      f.write('VA,6,9,11,13,33'+'\n')
      f.write('allsel,all'+'\n')
      #Cutting
      f.write('WPCSYS,-1,0 \n')
      f.write('wprota,,-90 \n')
      f.write('vsbw,all'+'\n')
      f.write('wpoff,,,'+str(C_R+bcl)+'\n') 
      f.write('vsbw,all'+'\n')
      f.write('wprota,,90 \n') 
      f.write('wpoff,,,'+str(B_R+z)+'\n')  
      f.write('vsbw,9'+'\n') 
      f.write('allsel,all'+'\n')
      #cutting brace crown 
      f.write('WPCSYS,-1,0 \n')
      f.write('wprota,,-90 \n')  
      f.write('wpoff,,,'+str(y_w[9]+a_b)+'\n')
      f.write('lsbw,71'+'\n')
      f.write('wpoff,,,'+str(b_b-a_b)+'\n')
      f.write('lsbw,18'+'\n')
      f.write('wpoff,,,'+str(b_b-a_b)+'\n')
      f.write('lsbw,37'+'\n')
      #cutting brace saddle line cut
      f.write('WPCSYS,-1,0 \n')
      f.write('wprota,,-90 \n')
      f.write('wpoff,,,'+str(y_w[0]+a_b)+'\n')
      f.write('lsbw,73'+'\n') 
      f.write('wpoff,,,'+str(b_b-a_b)+'\n')
      f.write('lsbw,83'+'\n') 
      f.write('wpoff,,,'+str(b_b-a_b)+'\n')
      f.write('lsbw,73'+'\n') 
      #cutting chord crown
      f.write('WPCSYS,-1,0 \n')
      f.write('wpoff,,,'+str(crown_r+a)+'\n')
      f.write('lsbw,82'+'\n')
      f.write('wpoff,,,'+str(b_cc-a)+'\n')
      f.write('lsbw,86'+'\n')
      f.write('wpoff,,,'+str(b_cc-a)+'\n')
      f.write('lsbw,82'+'\n')
      #cutting chord saddle
      f.write('WPCSYS,-1,0 \n')
      f.write('wprota,,,90'+'\n')
      f.write('wprota,,'+str((math.asin(saddle_r/C_R) + (a/C_R))*180/math.pi)+'\n')
      f.write('lsbw,40'+'\n')
      f.write('wprota,,'+str(((b_cs-a)/C_R)*180/math.pi)+'\n')   
      f.write('lsbw,82'+'\n')   
      f.write('wprota,,'+str(((b_cs-a)/C_R)*180/math.pi)+'\n') 
      f.write('lsbw,40'+'\n')
      #cutting chord volume
      f.write('WPCSYS,-1,0 \n') 
      f.write('wpoff,'+str(C_R)+',,'+str(B_R+(z))+'\n') 
      f.write('wpro,,-'+str(theta_cut)+'\n')
      f.write('vsbw,11'+'\n') 
      f.write('vsbw,2'+'\n') 
      f.write('WPCSYS,-1,0 \n')
      f.write('allsel,all'+'\n')
      f.write('wpoff,,,'+str(min(B_R+z+z_1,0.95*C_M_L))+'\n') 
      f.write('vsbw,7'+'\n')
      f.write('vsbw,8'+'\n') 
      f.write('WPCSYS,-1,0 \n')
      f.write('wpoff,,,'+str(B_R+z)+'\n') 
      f.write('vsbw,15'+'\n')
      f.write('vdele,7,,,1'+'\n')
      f.write('vdele,8,,,1'+'\n')
      f.write('vdele,16,,,1'+'\n')  
      f.write('vatt,1,,1,0'+'\n')
      # Meshing chord near weld
      f.write('lsel,s,,,17'+'\n')
      f.write('lsel,a,,,77'+'\n') 
      f.write('lsel,a,,,78'+'\n') 
      f.write('lsel,a,,,35'+'\n') 
      f.write('lsel,a,,,93'+'\n') 
      f.write('lsel,a,,,34'+'\n')
      f.write('lsel,a,,,48'+'\n') 
      f.write('lsel,a,,,19'+'\n') 
      f.write('lsel,a,,,46'+'\n') 
      f.write('lsel,a,,,45'+'\n') 
      f.write('lsel,a,,,47'+'\n')
      f.write('lsel,a,,,14'+'\n') 
      f.write('lsel,a,,,53'+'\n') 
      f.write('lsel,a,,,15'+'\n') 
      f.write('lsel,a,,,54'+'\n')  
      f.write('lesize,all,,,'+ str(num_C_T)+',,,,,0'+'\n')    
      f.write('lsel,s,,,73'+'\n')
      f.write('lsel,a,,,89'+'\n')
      f.write('lesize,all,,,'+str(num_weld_c)+',,,,,0'+'\n')
      f.write('lsel,s,,,87'+'\n')
      f.write('lsel,a,,,88'+'\n')
      f.write('lsel,a,,,90'+'\n')
      f.write('lsel,a,,,91'+'\n')
      f.write('lesize,all,,,'+str(c)+',,,,,0'+'\n')
      f.write('lsel,s,,,86'+'\n')
      f.write('lsel,a,,,82'+'\n')
      f.write('lesize,all,,,'+ str(num_ch_cut-(num_weld_c+2*c))+',,,,,0'+'\n')
      f.write('lsel,s,,,94'+'\n')
      f.write('lsel,a,,,95'+'\n')
      f.write('lsel,a,,,43'+'\n')
      f.write('lsel,a,,,81'+'\n')
      f.write('lesize,all,,,'+ str(num_ch_cut)+',,,,,0'+'\n')
      f.write('lsel,s,,,76'+'\n')
      f.write('lsel,a,,,65'+'\n')
      f.write('lesize,all,,,'+str(num_cir_crown)+',,,,,0'+'\n')
      f.write('lsel,s,,,80'+'\n')
      f.write('lsel,a,,,79'+'\n')
      f.write('lesize,all,,,'+str(num_cir_saddle)+',,,,,0'+'\n')
      f.write('vsweep,9'+'\n')
      f.write('vsweep,12'+'\n')
      f.write('vsweep,11'+'\n')
      f.write('vsweep,13'+'\n')  
      f.write('lsel,s,,,49'+'\n')
      f.write('lsel,a,,,27'+'\n')
      f.write('lsel,a,,,22'+'\n')
      f.write('lsel,a,,,42'+'\n')  
      f.write('lesize,all,,,'+ str(num_B_T)+',,,,,0'+'\n')
      f.write('allsel,all'+'\n')  
      f.write('vsweep,4'+'\n') 
      f.write('vsweep,5'+'\n') 
      f.write('vsweep,6'+'\n') 
      # Meshing brace adjacent to weld
      f.write('lsel,s,,,70'+'\n')
      f.write('lsel,a,,,51'+'\n')
      f.write('lsel,a,,,44'+'\n')
      f.write('lsel,a,,,26'+'\n')
      f.write('lesize,all,,,'+ str(num_B_T)+',,,,,0'+'\n')
      f.write('lsel,s,,,84'+'\n')
      f.write('lsel,a,,,85'+'\n')
      f.write('lsel,a,,,38'+'\n')
      f.write('lsel,a,,,71'+'\n')
      f.write('lesize,all,,,'+str(d)+',,,,,0'+'\n')
      f.write('lsel,s,,,20'+'\n')
      f.write('lsel,a,,,37'+'\n')
      f.write('lesize,all,,,'+str(num_weld_b)+',,,,,0'+'\n')
      f.write('lsel,s,,,18'+'\n')
      f.write('lsel,a,,,83'+'\n')
      f.write('lesize,all,,,'+ str(num_br_cut-(num_weld_b+2*d))+',,,,,0'+'\n')
      f.write('lsel,s,,,72'+'\n')
      f.write('lsel,a,,,74'+'\n')
      f.write('lesize,all,,,'+ str(num_br_cut)+',,,,,0'+'\n')
      f.write('vsweep,10'+'\n')
      # Meshing plate
      f.write('lsel,s,,,3'+'\n')
      f.write('lsel,a,,,7'+'\n')
      f.write('lsel,a,,,2'+'\n')
      f.write('lsel,a,,,10'+'\n')
      f.write('lesize,all,,,'+ str(e)+',,,,,0'+'\n')
      f.write('vsel,s,,,1'+'\n')
      f.write('vsweep,all'+'\n')
      # Meshing chord
      f.write('lsel,s,,,110'+'\n')
      f.write('lsel,a,,,112'+'\n')
      f.write('lsel,a,,,111'+'\n')
      f.write('lsel,a,,,109'+'\n')
      f.write('lesize,all,,,'+ str(num_ch_cut_1)+',,,,,0'+'\n')
      f.write('vsel,s,,,14'+'\n')  
      f.write('vsweep,14'+'\n')     
      f.write('EXTOPT,ESIZE,'+str(num_chord)+','+str(ratio_ch)+'\n')
      f.write('vsel,s,,,2'+'\n')
      f.write('vsweep,2,62,40'+'\n')
      # Meshing brace 
      f.write('EXTOPT,VSWE,AUTO,0'+'\n')
      f.write('EXTOPT,ESIZE,'+str(num_brace)+','+str(ratio_br)+'\n')
      f.write('vsel,s,,,3'+'\n')
      f.write('vsweep,3,44,21'+'\n')
      f.write('allsel,all'+'\n')
      #to get node numbering file 
      f.write('ksel,s,,,41,52,1'+'\n')
      f.write('nslk'+'\n')  
      f.write('''nwrite,nodeindex,'txt',,'''+'\n')
      #Axial load
      f.write('nsel,s,loc,x,0'+'\n')
      f.write('D,all, , , , , ,UX, , , , ,'+'\n')
      f.write('nsel,s,loc,z,0'+'\n')
      f.write('D,all, , , , , ,, ,UZ , , ,'+'\n')
      f.write('nsel,s,loc,y,0'+'\n')
      f.write('D,all, , , , , ,,UY , , , ,'+'\n')
      f.write('nsel,s,loc,z,'+str(C_M_L)+'\n')
      f.write('D,all, , , , , ,UX,UY,UZ , , ,'+'\n')
      f.write("nsel,s,loc,y,"+str(B_L-0.01)+","+str(B_L+0.01)+'\n')
      f.write('SF,all,PRES,1'+'\n')
      #Inplane Bending
      f.write('nsel,s,loc,x,0'+'\n')
      f.write('DSYM,SYMM,X, , '+'\n')
      f.write('nsel,s,loc,z,0'+'\n')
      f.write('DSYM,ASYM,Z,,'+'\n') 
      f.write('nsel,s,loc,y,0'+'\n')
      f.write('DSYM,SYMM,Y,,'+'\n') 
      f.write('nsel,s,loc,z,'+str(C_M_L)+'\n')
      f.write('D,all, , , , , ,UX,UY,UZ , , ,'+'\n')
      div = 4
      for i in range(div):
          f.write("nsel,s,loc,y,"+str(B_L-0.01)+","+str(B_L+0.01)+'\n')
          f.write("nsel,r,loc,z,"+ str(i*B_R/div)+","+str( 0.01+(i+1)*B_R/div)+'\n')
          f.write("SF,all,pres,-"+ str((i+0.5)/div)+'\n')
      #Outplane Bending
      f.write('nsel,s,loc,x,0'+'\n')
      f.write('DSYM,ASYM,X, , '+'\n')
      f.write('nsel,s,loc,z,0'+'\n')
      f.write('DSYM,SYMM,Z,,'+'\n') 
      f.write('nsel,s,loc,y,0'+'\n')
      f.write('DSYM,SYMM,Y,,'+'\n') 
      f.write('nsel,s,loc,z,'+str(C_M_L)+'\n')
      f.write('D,all, , , , , ,UX,UY,UZ , , ,'+'\n')
      div = 4
      for i in range(div):
          f.write("nsel,s,loc,y,"+str(B_L-0.01)+","+str(B_L+0.01)+'\n')
          f.write("nsel,r,loc,x,"+ str(i*B_R/div)+","+str( 0.01+(i+1)*B_R/div)+'\n')
          f.write("SF,all,pres,-"+ str((i+0.5)/div)+'\n')
      #Analysis
      f.write('allsel,all'+'\n')
      f.write('/sol'+'\n')
      f.write('antype,0'+'\n')
      f.write('solve'+'\n') 
      f.write('save\n')
      f.write('finish\n')
      if count<len(df)-1:
          f.write('BOPTN,NWARN,1\n')
          f.write('/UIS,MSGPOP,3 \n') 
          f.write('KEYW,PR_SGVOF,1  \n')
          f.write('/NERR,5,10000, ,1,5, \n')
          f.write('/clear \n')
          f.write('''/CWD,'''+"'"+directory+ str(count+1)+"'"+ '\n')
          f.write('''/filnam, '''+ str(count+1)+'\n')
          f.write('''/title, '''+ str(count+1)+'\n')
          f.write('''/INPUT,'model','txt','''+"'"+directory+ str(count+1)+'''\\',, 0'''+'\n')      
# saving the dataframe
str_folder = str(folder)+'_without_result'+'.csv'
df.to_csv(str_folder)

