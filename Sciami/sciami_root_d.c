#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <math.h>
#include <TGraph.h>
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <TMath.h>
#define N 42925
#define V 429240
#define NBINS 3000


typedef struct stuttura1{
	double n_t[NBINS];
	double rate_t[NBINS];
	double Delta_t;
}isto_temps;

typedef struct stuttura2{
	double mu_cps;
	double sigma_cps;
}analisi;


double max(double  t[], int len){
    double max=t[0];
    for(int i=0;i<len;i++){
      if(max<=t[i]){
        max=t[i];

        }
    }

    return max;
}

double min(double t[], int len){
    double min=t[0];
    for(int i=0;i<len;i++){
      if(min>=t[i]){
        min=t[i];

        }
    }

    return min;
}

isto_temps istogramma_tempi(double t[N], int nbins, int len, int numero_figura, char ref_istogramma[30], double numero_acquisizione, char nome_figura[30]){
    /*
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
    */
    isto_temps results; //struct per inserire i risultati
    printf("%lf \n", t[len-1]);
    TH1F *histogram= new TH1F("histogram","Istogramma",nbins, t[0],t[len-1]+1);
    for(int i=0; i<len;i++){
        histogram->Fill(t[i]);
    }
    new TCanvas;
    histogram->Draw();
    
 
    double medie_bins[nbins];
    for(int i=0;i<nbins;i++){
        medie_bins[i]=histogram->GetBinCenter(i+1); //è giusto??? //media dei bins
        results.n_t[i]=histogram->GetBinContent(i+1);
 
        
    }
    
    double Deltat=(medie_bins[10]-medie_bins[9]); 
 
    for(int i=0; i<nbins; i++){
        results.rate_t[i]=results.n_t[i]/Deltat;
    }
    results.Delta_t=Deltat;
    
    double n_tot = 0;
    double x[nbins];
    for(int i=0;i<nbins;i++){
        n_tot+=results.n_t[i];
        x[i]=medie_bins[i]/3600;
    }      
    
    new TCanvas;
    TGraph *graphic=new TGraph(nbins, x, results.rate_t);
    graphic->SetMarkerStyle(20);
    graphic->Draw();
    double maxt=max(t,len);
    double dt=maxt/nbins;
    double rate = n_tot/maxt;  
    
    printf("Numero totale conteggi %s è ntot=%lf in Deltat = %lf s ; rate = %lf cps \n", ref_istogramma, n_tot,dt,rate);
    return results;
}

isto_temps istogramma_tc(double t[N], int nbins, int len, int numero_figura, char ref_istogramma[30], double numero_acquisizione, char nome_figura[30]){
    /*
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
    */
    isto_temps results; //struct per inserire i risultati
    printf("%lf \n", t[len-1]);
    TH1F *histogram= new TH1F("histogram","Istogramma",nbins, t[0],t[len-1]+1);
    for(int i=0; i<len;i++){
        histogram->Fill(t[i]);
    }
    new TCanvas;
    histogram->Draw();
    
 //è giusto?? cambia in caso  (nbins)
    double medie_bins[nbins];
    for(int i=0;i<nbins;i++){
        medie_bins[i]=histogram->GetBinCenter(i+1); //è giusto??? //media dei bins
        results.n_t[i]=histogram->GetBinContent(i+1);
 
        
    }
    
    double Deltat=(medie_bins[10]-medie_bins[9]); 
 
    for(int i=0; i<nbins; i++){
        results.rate_t[i]=results.n_t[i]/Deltat;
    }
    results.Delta_t=Deltat;
    
    double n_tot = 0;
    double x[nbins];
    for(int i=0;i<nbins;i++){
        n_tot+=results.n_t[i];
        x[i]=medie_bins[i]/3600;
    }      
    
    new TCanvas;
    TGraph *graphic=new TGraph(nbins, x, results.n_t);
    graphic->SetMarkerStyle(20);
    graphic->Draw();
    double maxt=max(t,len);
    double dt=maxt/nbins;
    double rate = n_tot/maxt;  
    
    printf("Numero totale conteggi %s è ntot=%lf in Deltat = %lf s ; rate = %lf cps \n", ref_istogramma, n_tot,dt,rate);
    return results;
}



analisi analisi_cps(double counts[NBINS], int nbins,double Delta_t, int numero_figura, char ref_istogramma[30], double numero_acquisizione, char nome_figura[30]){
    analisi results;
/*   Parameters
    ----------
    counts : array di conteggi
        contiene i conteggi restituiti dall'istogramma per i tempi
    Delta_t : float
        intervallo temporale corrispondente ai conteggi: rate = counts/delta_t.
    numero_figura : numero della figura che verrà creata
    ref_istogramma : specifiche dell'istogramma nel titolo della figura
     numero_acquisizione : float
         numero di ore di acquisizione totale
    save_fig : bool, optional
        se True, la figura viene salvata. The default is False.
    nome_figura : string, optional
        nome del file in cui verrà salvata la figura. The default is None.
    
        

    Returns
    -------
    mu_cps : float
        media della poissoniana uscita dal fit divisa per delta_t: dà la media del rate
    sigma_mu_cps : float
         errore sulla media della poissoniana uscita dal fit divisa per delta_t: dà la sigma della media del rate.

    */ 
    printf("\n------Stima di pdf %s-------\n",ref_istogramma);
    
    //stima di media per confrontarla con risultati fit
    double tot=0;
    for(int i=0; i<NBINS;i++){
        tot+=counts[i];
    }
    double mean_exp=tot/NBINS;
    double tot1=0;
    for(int i=0; i<NBINS;i++){
        tot1+=pow((counts[i]-mean_exp),2)/NBINS;
    }
    printf("tot1=%lf \n", tot1);
    double std_exp=pow(tot1,0.5);
    
    
    //nota: anche mean e std hanno una loro varianza (sono statistiche)
    double minimo=min(counts,NBINS);
    double massimo=max(counts,NBINS);
    printf("Media dei counts = %lf pm %lf ", mean_exp,std_exp);
    TH1F *histogram= new TH1F("histogram","Istogramma",nbins, minimo ,massimo+1);
    for(int i=0; i<NBINS;i++){ //NBINS O 100????????
        histogram->Fill(counts[i]);
    }
    new TCanvas;
    histogram->Draw();
    

    //preparazione dei dati per il fit: prendo la media dei bins (intera!)
   
    double medie_bins[NBINS];
    double occorrenze_rate[NBINS];
    double x[NBINS];
    double y[NBINS];
    int len_x=0;
    int index[NBINS];
    for(int i=0;i<NBINS;i++){
        
        occorrenze_rate[i]=histogram->GetBinContent(i+1);
        printf("%lf \n",occorrenze_rate[i]);
        medie_bins[i]=histogram->GetBinCenter(i+1);
        if(occorrenze_rate[i]>0){
            x[len_x]=medie_bins[i];
            index[len_x]=i;
            len_x++;
            
        
        }
        
    }
    double sigma_y[len_x];
    double sum_y=0;//occorrenze totali

    for(int j=0; j<len_x;j++){
        y[j]=occorrenze_rate[index[j]];
        sigma_y[j]=pow(y[j], 0.5); //errore su y: assumo conteggi poissoniani
        sum_y+=y[j];
        
    }
    
    //normalizzo y
    /*for(int j=0; j<len_x;j++){
        y[j]=y[j]/sum_y;
        sigma_y[j]=sigma_y[j]/sum_y;
        
    }*/
    //fit della poissoniana 
    
    printf("%d \n y",len_x);
    new TCanvas;
    TGraphErrors *graphic = new TGraphErrors(len_x,x,y,nullptr,sigma_y);
    
    graphic->SetMarkerStyle(20);
    char a[300]="Istogramma";
    char * title;
    title=strcpy(a,ref_istogramma);
    
    graphic->SetTitle(title);
    // Define a function 
    //TF1 *f = new TF1(“myfunc”, "[0]*sin(x)" , parameters);

    graphic->GetHistogram()->GetYaxis()->SetTitle("Occorrenze in ogni bin");
    graphic->GetHistogram()->GetXaxis()->SetTitle("Conteggi in #Delta t s"); //come metto numero qua?????
    graphic->Fit("gaus"); //NON RIESCO A TROVARE GAUS
    TLegend *legend = new TLegend(1.,1.,1.,1.,"Lab. Lesson 1");
    legend->AddEntry(graphic,"Exp. Points","PE");
    //egend->AddEntry(gaus,"Th. Law","L");
    legend->Draw();

    graphic->Draw("APE");
    gStyle->SetOptFit(1111);
    
    /*
    

    popt, pcov = curve_fit(fit_function, x, y, sigma = sigma_y, p0=[np.mean(counts), occorrenze_totali*0.02])
    
    # estraggo i dati del fit
    mu_fit=popt[0]
    sigma_mu_fit = np.sqrt(pcov.diagonal())[0]
    var_fit = np.sqrt(mu_fit)  # la varianza di una poissoniana è uguale a mu
    sigma_var_fit = sigma_mu_fit/np.sqrt(mu_fit)/2


    mu_cps = mu_fit/Delta_t
    sigma_mu_cps = sigma_mu_fit/Delta_t

    print(f"Media della poissoniana fittata:  = {mu_fit:.3f} pm {sigma_mu_fit:.3f}")
    print(f"Varianza della poissoniana fittata:  = {var_fit:.3f} pm {sigma_var_fit:.3f}")

    print(f"Media della poissoniana divisa per Delta t: rate = {mu_cps:.3f} pm {sigma_mu_cps:.3f}")
    

    # Grafico istogramma conteggi al secondo: nel grafico divido le x per delta_t
    
    plt.figure(numero_figura)
    ax=plt.subplot()
    plt.title(f'Istogramma {ref_istogramma} - acquisizione di {numero_acquisizione:.0f} h')
    plt.ylabel('Occorrenze in ogni bin')
    plt.xlabel(f'Conteggi in $\Delta t$ = {Delta_t:.1f} s')
    #plt.xlim([-3, max(counts)+1])
    #plt.ylim([0,max(y)*1.3])
    
    t = plt.text(0.02, 0.66, 
                 f"$\overline{{x}}$ = {mean_exp:.1f} \n"
                 +f"$\sigma_{{x}}$ = {std_exp:.2f} \n"
                 +f"$\mu$ = {mu_fit:.1f} $\pm$ {sigma_mu_fit:.1f} \n"
                 +f"$\sigma$ = {var_fit:.2f} $\pm$ {sigma_var_fit:.2f}\n"
                 +f"cps = {mu_cps:.2f} $\pm$ {sigma_mu_cps:.2f}", transform=ax.transAxes)
    t.set_bbox(dict(facecolor='red', alpha=0.5))
    
    plt.bar( x, y, width = width_bin, align='center', label='Dati' , zorder=1, capsize=2)
    plt.errorbar(x, y, yerr = sigma_y, fmt='.',  zorder=3, capsize=2, color='k')
    x_plot = np.arange(min(x)-1, max(x)+1, 1)
    plt.plot(x_plot,fit_function(x_plot, *popt), ls = '-', label='Fit', color='red', zorder=2)
    #plt.legend(loc='best')
    
    */
    results.mu_cps=0;
    results.sigma_cps=0;
    
    return results;
}
void ciao(){

    double Epsilon_tel1= 0.573;
    double N_tel1 = 13.666; //numero di raggi cosmici al secondo attesi nel tel1
    double N_tel2 = 32.13;//numero di raggi cosmici al secondo attesi nel tel2
    double Epsilon_tel2 = 0.246;
    int ch[N];
    double t[N];
    int i;
    int len_t=0;
    FILE * fin=fopen("/home/sarinaleggy/Documenti/root_mio/sciami_291122_prova4.dat","r");
    for(i=0; i<N; i++) {
      fscanf(fin,"%d %lf",ch+i,t+i);
    }
    fclose(fin);
    
    

    int len_T=0;
    double t1[N], t2[N],tc[N];
    int ch1[N], ch2[N];
    double Delta_t_contatore=0.01;
    double Delta_t_coincidenza=0.0000001;
    double T[N];
    T[0]=t[0];
    for(i=1; i<N;i++){
      if(t[i-1]<=t[i] ){
        T[i]=t[i];
        len_T++;
      }
      else if( t[i-1]>t[i] ){
        T[i]=t[i-1]+Delta_t_contatore;
        len_T++;
        double Delta=t[i-1]-t[i];
        for(int j=i+1;j<N;j++){
          t[j]=t[j]+Delta;

        }

      }
    }


    int len_t2=0; 
    int len_t1=0;
    int len_tc=0;
    int len_t12=0;
    double Delta_t12[N];
 
    for(int k=0;k< len_T-1;k++){
        if(ch[k]==2){
            t2[len_t2]=T[k];
            len_t2++;
        }
        else if(ch[k]==1){
            t1[len_t1]=T[k];
            len_t1++;

        }
        if((ch[k]==1) && (k>8) && (k<len_T-9)){
            double Delta_tmin = 1000;
            for(int j=k-8;j<k+9;j++){
                if (ch[j] == 2){
                    if (fabs(T[j]-T[k])<Delta_tmin){
                        Delta_tmin = fabs(T[j]-T[k]);
                    }
                }
            }
            Delta_t12[len_t12]= Delta_tmin;
            len_t12++;

            if(Delta_tmin < Delta_t_coincidenza){
                tc[len_tc]=T[k]; 
                len_tc++;
                
    
            }
        }
    
    }
    printf("lent12=%d \n", len_t12);
    /*for(int i=1; i<len_t1;i++){
        if(t1[i-1]>t1[i]){
            printf("errore");
        }

    }
    printf("lunghezza t1: %d \n",len_t1);
    printf("lunghezza t2: %d \n",len_t2);
    */
         
    /*
    TH1F *hist1= new TH1F("hist1","Istogramma1",300, t1[0],t1[len_t1]);
    for(int i=0; i<len_t1;i++){
        hist1->Fill(t1[i]);
    }
   
    hist1->Draw();
    
    TH1F *hist2= new TH1F("hist2","Istogramma2",300, t2[0],t2[len_t2]);
    
    for(int i=0; i<len_t2;i++){
        hist2->Fill(t2[i]);
    }
    new TCanvas;
    hist2->Draw();
    double a[3000];
    for(int i=0;i<3000;i++){
        a[i]=hist2->GetBinContent(i+1);
        
    }
    TH1F *hist3= new TH1F("hist3","Istogramma3",300, 0,1500);
    for(int i=0; i<len_t2;i++){
        hist3->Fill(a[i]);
    }
    new TCanvas;
    hist3->Draw();
    
   */
    int n_acquis=0;

    double max_t2,max_t1;
    max_t2=max(t2,len_t2);
    max_t1= max(t1,len_t1);
    double durata_acquisizione;
    if(max_t2>=max_t1){
        durata_acquisizione = max_t2/3600; //ore
    }
    else if(max_t2<max_t1){
        durata_acquisizione = max_t1/3600; //ore
    }
    
    //manca tc
    isto_temps parte1, parte2,partetc;
    analisi rate1, rate2;

    parte1= istogramma_tempi(t1,NBINS,len_t1, n_acquis+1, (char*)"telescopio 1",durata_acquisizione, (char*)"ist_tel1");

    parte2= istogramma_tempi(t2, NBINS,len_t2,n_acquis+2, (char*)"telescopio 2",durata_acquisizione, (char*)"ist_tel2");

    partetc= istogramma_tc(tc, 100,len_tc, n_acquis+3,(char*)"coincidenze", durata_acquisizione,(char*)"ist_coinc");
    
    //Fit dell'istogramma dei conteggi di eventi per unità di tempo per i due canali

    rate1= analisi_cps(parte1.n_t,100, parte1.Delta_t, n_acquis+4, (char*)"conteggi telescopio 1",durata_acquisizione,  (char*)"ist_poiss_tel1");
    
    rate2= analisi_cps(parte2.n_t,100, parte2.Delta_t, n_acquis+5,(char*)"conteggi telescopio 2", durata_acquisizione, (char*)"ist_poiss_tel2");

    //Fit delle differenze di tempo tra un evento e l'altro

    /*analisi_delta_t(t1, n_acquis+6, ref_istogramma=r"$\Delta t$ telescopio 1", save_fig=True, nome_figura=f"ist_dtempi_tel1_{durata_acquisizione:.0f}h", numero_acquisizione = durata_acquisizione)

    analisi_delta_t(t2, n_acquis+7, ref_istogramma=r"$\Delta t$ telescopio 2", save_fig=True, nome_figura=f"ist_dtempi_tel2_{durata_acquisizione:.0f}h", numero_acquisizione = durata_acquisizione)

    analisi_dt_12(Delta_t12, n_acquis+8, ref_istogramma = r"$\Delta t$ telescopi", save_fig=True, nome_figura=f"ist_dtempi_tel12_{durata_acquisizione:.0f}h", numero_acquisizione=durata_acquisizione)
    */
}