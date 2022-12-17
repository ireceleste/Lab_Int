import numpy as np
from matplotlib import pyplot as plt
import pylab
from time import time

def istogramma_tempi(t, nbins):
    """
    

    Parameters
    ----------
    t : array dei tempi
    nbins : numero di bins per l'istogramma dei tempi
    numero_figura : numero della figura in cui andrà l'istogramma
    ref_istogramma : stringa da mettere nel titolo su che istogramma è
    numero_acquisizione : numero di ore di acquisizione totale

    Returns
    -------
    n_t : numero di conteggi in ogni bin
    rate_t : rate dei conteggi per ogni bin di t
        array che contiene il numero di eventi registrati in ogni bin temporale
        diviso per la larghezza dell'intervallo temporale'
    """
    
   
    n_t, bins = np.histogram(t, nbins)
    
    Delta_t=bins[10]-bins[9]  # in secondi
    n_tot = sum(n_t)
    
    mean_bins = np.array([(bins[i]+bins[i+1])/2 for i in range(len(bins)-1)])
    return n_t, mean_bins


def istogramma_eff(nd, nt,bins, n_fig, n_tel, n_pm):
    time= bins/3600 #in ore
    epsilon = nt/nd
    depsilon=np.sqrt(nt*(1-nt/nd))/nd
    plt.figure(n_fig)
    plt.title(f"Efficienza PMT {n_pm} telescopio {n_tel}")
    plt.xlabel("Tempo [h]")
    plt.ylabel("Efficienza")
    plt.errorbar(time, epsilon, depsilon, fmt='.', ls='', capsize=2)
    return






# =============================================================================
# Prima acquisizione di prova
# =============================================================================

t0 = time()
n_acquis=0
nbins = 50

nome_file= 'efficienze_131222_prova2'

t1_12= np.load('Dati_fixed/'+nome_file+'_t1_fixed.npy')
t2_123= np.load('Dati_fixed/'+nome_file+'_t2_fixed.npy')
t1_123= np.load('Dati_fixed/'+nome_file+'_t3_fixed.npy')
t1_13= np.load('Dati_fixed/'+nome_file+'_t4_fixed.npy')
t1_23= np.load('Dati_fixed/'+nome_file+'_t5_fixed.npy')
t2_12= np.load('Dati_fixed/'+nome_file+'_t6_fixed.npy')
t2_13= np.load('Dati_fixed/'+nome_file+'_t7_fixed.npy')
t2_23= np.load('Dati_fixed/'+nome_file+'_t8_fixed.npy')


tc= np.load('Dati_fixed/'+nome_file+'_tc_fixed.npy')
Delta_t12= np.load('Dati_fixed/'+nome_file+'_Delta_t12_fixed.npy')

durata_acquisizione = max([max(t1_12), max(t2_123)])/3600 # ore

print(f"Acquisizione di {durata_acquisizione} ore")
c1_12, bins1_12 = istogramma_tempi(t1_12, nbins)
c2_123,bins2_123 = istogramma_tempi(t2_123, nbins)
c1_123, bins1_23 = istogramma_tempi(t1_123, nbins)
c1_13, bins1_13 = istogramma_tempi(t1_13, nbins)
c1_23, bins1_23 = istogramma_tempi(t1_23, nbins)
c2_12, bins2_12 = istogramma_tempi(t2_12, nbins)
c2_13, bins2_13 = istogramma_tempi(t2_13, nbins)
c2_23, bins2_23 = istogramma_tempi(t2_23, nbins)

e1_1 = np.mean(c1_23/c1_123)
e1_2 = np.mean(c1_13/c1_123)
e1_3 = np.mean(c1_12/c1_123)

e2_1 = np.mean(c2_23/c2_123)
e2_2 = np.mean(c2_13/c2_123)
e2_3 = np.mean(c2_12/c2_123)

print("------ Telescopio 1 efficienze------")
print(f"PMT1 = {1/e1_1}")
print(f"PMT2 = {1/e1_2}")
print(f"PMT3 = {1/e1_3}")

print("------ Telescopio 2 efficienze------")
print(f"PMT1 = {1/e2_1}")
print(f"PMT2 = {1/e2_2}")
print(f"PMT3 = {1/e2_3}")

n_fig=n_acquis+1

istogramma_eff(c2_12, c2_123, bins2_12, n_fig = n_fig, n_tel=2, n_pm=3) #epsilon 3 telescopio 2
n_fig+=1
istogramma_eff(c2_13, c2_123, bins2_12,n_fig = n_fig, n_tel=2, n_pm=2) #epsilon 2 telescopio 2
n_fig+=1

istogramma_eff(c2_23, c2_123, bins2_12,n_fig = n_fig, n_tel=2, n_pm=1) #epsilon 1 telescopio 2
n_fig+=1

istogramma_eff(c1_12, c1_123, bins1_12,n_fig = n_fig, n_tel=1, n_pm=3) #epsilon 3 telescopio 1
n_fig+=1

istogramma_eff(c1_13, c1_123, bins1_12,n_fig = n_fig, n_tel=1, n_pm=2) #epsilon 2 telescopio 1
n_fig+=1

istogramma_eff(c1_23, c1_123, bins1_12,n_fig = n_fig, n_tel=1, n_pm=1) #epsilon 1 telescopio 1
n_fig+=1


plt.show()
