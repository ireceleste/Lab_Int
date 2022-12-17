import numpy as np




def funzione_dati(t_input, ch):
    """
    

    Parameters
    ----------
    t_input : array di float
        tempo in input dal file .dat
    ch : array di 1 e 2
        canale in cui si è registrato l'evento
    
    Nota: t_input è preso direttamente dal contatore: dato che si riazzera in modo random
    abbiamo fatto un ciclo for per rendere sequenziali i tempi. A questo scopo abbiamo
    creato un array T che contiene i tempi sistemati
    
    Returns
    -------
    t1 : array di tempi del canale 1
    t2 : array di tempi del canale 2
    tc : array di tempi in cui si è registrata una coincidenza tra canale 1 e 2
        registra i tempi per cui il canale 1 e 2 sono scattati a un delta t Delta_t_coinc
    Delta_t12: array di delta t tra telescopio 1 e 2
    """
    t1=np.array([]) # 1+2 del 1
    t2=np.array([]) # 1+2+3 del 2
    t3=np.array([]) # 1+2+3 del 1
    t4=np.array([]) # 1+3 del 1
    t5=np.array([]) # 2+3 del 1
    t6=np.array([]) # 1+2 del 2
    t7=np.array([]) # 1+3 del 2
    t8=np.array([]) # 2+3 del 2
    
    tc= np.array([])
    Delta_t12=np.array([])
    t=t_input
    Delta_t_contatore =0.01 #s : tempo approssimativo tra una misura e l'altra
    # quando il contatore si riazzera -  è circa il tau dell'istogramma dei delta t
    Delta_t_coinc = 100e-9 #s
    
    # ciclo for per sistemare i dati per i tempi
    T = np.array([t[0]])
    
    for i in range(1, len(t)):

        if( t[i-1]<t[i]):
            T = np.append(T, t[i])
        elif( t[i]<t[i-1]):
            T = np.append(T, t[i-1]+ Delta_t_contatore)
            t=t+t[i-1]-t[i]
        else:
            T=np.append(T,t[i])
            #print("Aggiunto 0", i, t[i], t[i-1], ch[i], ch[i-1])

        if (i%100000==0):
            print(f"Fixing times {i/len(t)*100:.0f}%")
   
    print(f"Durata acquisizione = {max(T)} s")

    #ciclo for per spacchettare i tempi tra canale 1 e 2 e per registrare le coincidenze
    print(len(T), len(t))
    
    for i in range(len(T)-1):
        if (ch[i]==1):
             t2=np.append(t2, T[i])
        elif(ch[i]==2):
            t1=np.append(t1, T[i])
        elif (ch[i]==3):
            t3=np.append(t3, T[i])
        elif (ch[i]==4):
            t4=np.append(t4, T[i])
        elif (ch[i]==5):
            t5=np.append(t5, T[i])
        elif (ch[i]==6):
            t6=np.append(t6, T[i])
        elif (ch[i]==7):
            t7=np.append(t7, T[i])
        elif (ch[i]==8):
            t8=np.append(t8, T[i])
        else:
            print('Errore: canale non esistente')
        
        if (ch[i]==2 and i>8 and i<len(T)-9):
            Delta_tmin = 1000
            for j in range(i-8, i+9):
                if (ch[j] == 1):
                    if (abs(T[j]-T[i])<Delta_tmin):
                        Delta_tmin = abs(T[j]-T[i])
            if(Delta_t12<900):
                Delta_t12 = np.append(Delta_t12, Delta_tmin)
            
            if (Delta_tmin < Delta_t_coinc):
                print(Delta_tmin, i)
                tc=np.append(tc, T[i])
                
        if (i%100000==0):
            print(f"Unpacking {i/len(t)*100:.0f}%")
        

    return t1, t2, t3, t4, t5, t6, t7, t8, tc, Delta_t12

       


nome_file= 'efficienze_131222_prova2'
path_file_in='Dati/'+nome_file+'.dat'


ch, t= np.loadtxt(path_file_in, unpack='True')
print(len(t))
#t=t[:500000]

t1, t2, t3, t4, t5, t6, t7, t8, tc, Delta_t12= funzione_dati(t, ch)

np.save('Dati_fixed/'+nome_file+'_t1_fixed.npy', t1)
np.save('Dati_fixed/'+nome_file+'_t2_fixed.npy', t2)
np.save('Dati_fixed/'+nome_file+'_t3_fixed.npy', t3)
np.save('Dati_fixed/'+nome_file+'_t4_fixed.npy', t4)
np.save('Dati_fixed/'+nome_file+'_t5_fixed.npy', t5)
np.save('Dati_fixed/'+nome_file+'_t6_fixed.npy', t6)
np.save('Dati_fixed/'+nome_file+'_t7_fixed.npy', t7)
np.save('Dati_fixed/'+nome_file+'_t8_fixed.npy', t8)


np.save('Dati_fixed/'+nome_file+'_tc_fixed.npy', tc)
np.save('Dati_fixed/'+nome_file+'_Delta_t12_fixed.npy', Delta_t12)


