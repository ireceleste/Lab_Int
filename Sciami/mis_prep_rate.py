import numpy as npfrom matplotlib import pyplot as pltimport pylabfrom scipy.optimize import curve_fit# =============================================================================# Punto di lavoro del telescopio 1# =============================================================================Al1=1475Al2=1440Al3=1800cps1=np.array([4709, 4602, 3311])cps2=np.array([33891, 33871, 35209])cps3=np.array([1477550, 1458583, 1453225])dcps1=np.sqrt(cps1)dcps2=np.sqrt(cps2)dcps3=np.sqrt(cps3)#epsilon 2cps13=474cps123=452epsilon2=cps123/cps13dcps13=np.sqrt(cps13)dcps123=np.sqrt(cps123*(1-epsilon2)) depsilon2=dcps123/cps13#epsilon 1cps23=742cps123=446epsilon1=cps123/cps23dcps23=np.sqrt(cps23)dcps123=np.sqrt(cps123*(1-epsilon1)) depsilon1=dcps123/cps23#epsilon 3cps12=722cps123=380epsilon3=cps123/cps12dcps12=np.sqrt(cps12)dcps123=np.sqrt(cps123*(1-epsilon3)) depsilon3=dcps123/cps12epsilon=epsilon1*epsilon2depsilon=np.sqrt((epsilon1*depsilon2)**2+(epsilon2*depsilon1)**2)print(f"------Efficienze Telescopio 1------ \nPMT1 = {epsilon1:.3f} pm {depsilon1:.3f} \nPMT2 = {epsilon2:.3f} pm {depsilon2:.3f} \nPMT3 = {epsilon3:.3f} pm {depsilon3:.3f}")print(f"\nEfficienza epsilon 1 * 2 = {epsilon:.3f} pm {depsilon:.3f} ")print("\n-----------Telescopio 1 cose attese-------------")x=41 #cmy=20 #cmS= x*y #cm2N = 1*S/(60) # attesi 1 raggio cosmico in 1 minuto al cm^2print("Numero di raggi cosmici al secondo attesi nel telescopio 1: ", N)print("Numero di raggi cosmici misurati al secondo attesi nel telescopio 1:", N*epsilon)#stimiamo le coincidenze accidentali, mettendo stima brutale dei ratew=50e-9wmin=2e-9Dt=100R1=np.mean(cps1)/DtR2=np.mean(cps2)/DtR3 = np.mean(cps3)/DtdR1=np.sqrt(sum(dcps1**2))/DtdR2=np.sqrt(sum(dcps2**2))/DtdR3=np.sqrt(sum(dcps3**2))/Dtcps_doppie_acc=R1*R2*(2*w-2*wmin)dcps_doppie_acc=(2*w-2*wmin) * np.sqrt((R1*dR2)**2+(R2*dR1)**2)print(f"CPS PMT1 = {R1:.0f} \pm {dR1:.0f}")print(f"CPS PMT2 = {R2:.0f} \pm {dR2:.0f}")print(f"CPS PMT3 = {R3:.0f} \pm {dR3:.0f}")print(f'Stima delle doppie accidentali al secondo: {cps_doppie_acc:.6f} pm {dcps_doppie_acc:.6f}')R1_tel1 = R1R2_tel1 = R2dR1_tel1 = dR1dR2_tel1 = dR2# =============================================================================# Punto di lavoro del telescopio 2# =============================================================================Al1=1740Al2=1750Al3=1711cps1=np.array([7001, 7061, 7032])cps2=np.array([18192, 18068, 18235])cps3=np.array([15646, 12937, 34845])dcps1=np.sqrt(cps1)dcps2=np.sqrt(cps2)dcps3=np.sqrt(cps3)#epsilon 2cps13=871cps123=807dcps13=np.sqrt(cps13)dcps123=np.sqrt(cps123*(1-epsilon2)) depsilon2=dcps123/cps13epsilon2=cps123/cps13#epsilon 1cps23=1771cps123=825epsilon1=cps123/cps23dcps23=np.sqrt(cps23)dcps123=np.sqrt(cps123*(1-epsilon1)) depsilon1=dcps123/cps23epsilon2=cps123/cps13#epsilon 3cps12=1491cps123=832epsilon3=cps123/cps12dcps12=np.sqrt(cps12)dcps123=np.sqrt(cps123*(1-epsilon3)) depsilon3=dcps123/cps12epsilon=epsilon1*epsilon2*epsilon3depsilon=np.sqrt((epsilon1*depsilon2)**2+(epsilon2*depsilon1)**2)print(f"\n------Efficienze Telescopio 2------ \nPMT1 = {epsilon1:.3f} pm {depsilon1:.3f} \nPMT2 = {epsilon2:.3f} pm {depsilon2:.3f} \nPMT3 = {epsilon3:.3f} pm {depsilon3:.3f} \nProdotto = {epsilon}")# =============================================================================# Incertezza sul rate# =============================================================================print("\n-----------Telescopio 2 cose attese-------------")x=48.2 #cmy=40 #cmS= x*y #cm2N = 1*S/(60) # attesi 1 raggio cosmico in 1 minuto al cm^2print("Numero di raggi cosmici al secondo attesi nel telescopio 2: ", N)print("Numero di raggi cosmici misurati al secondo attesi nel telescopio 2:", N*epsilon)#stimiamo le coincidenze accidentali, mettendo stima brutale dei ratew=50e-9wmin=2e-9Dt=100R1=np.mean(cps1)/DtR2=np.mean(cps2)/DtR3=np.mean(cps3)/DtdR1=np.sqrt(sum(dcps1**2))/DtdR2=np.sqrt(sum(dcps2**2))/DtdR3=np.sqrt(sum(dcps3**2))/Dtprint(f"CPS PMT1 = {R1:.0f} \pm {dR1:.0f}")print(f"CPS PMT2 = {R2:.0f} \pm {dR2:.0f}")print(f"CPS PMT3 = {R3:.0f} \pm {dR3:.0f}")cps_triple_acc=R1*R2*R3*(2*w-2*wmin)*wdcps_triple_acc=(2*w-2*wmin) * w * np.sqrt((R1*R3*dR2)**2+(R3*R2*dR1)**2 +(R1*R2*dR3)**2)print("Stima delle triple accidentali al secondo: ", cps_triple_acc)print("\n-----------Stima coincidenza accidentale tra i 2 telescopi-------------")R1_tel2 = R1R2_tel2 = R2R3_tel2 = R3dR3_tel2 = dR3dR1_tel2 = dR1dR2_tel2 = dR2t_coincid_telescopi = 100e-9t_acquis = 40*3600#n_coinc_acc = R1_tel1*R2_tel1*R1_tel2*R2_tel2*R3_tel2 *(2*w-2*wmin)*(2*w-2*wmin)*t_coincid_telescopi *t_acquis *w#n_coinc_acc = cps_triple_acc * cps_doppie_acc * t_coincid_telescopi * t_acquisrate_tel1 = 11.876186339752811rate_tel2 = 16.04636922292368n_coinc_acc = rate_tel1 * rate_tel2 * t_coincid_telescopi * t_acquisrate_coinc_acc = rate_tel1 * rate_tel2 * t_coincid_telescopi print("N coincidenze accidentali telescopi in 40 ore s", n_coinc_acc)print("Rate coincidenze accidentali telescopi al secondo", rate_coinc_acc)