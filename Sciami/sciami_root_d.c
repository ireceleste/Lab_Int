#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH2.h>
#include <TF1.h>
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
#include <cstring>
#include <vector>
#include <string>
#include <iostream>
#include <bits/stdc++.h>
#include <algorithm>
#include <cmath>
#include <TAxis.h>
#include <iomanip>
#include <TBox.h>
#include <TPaveText.h>
using namespace std;
#define NBINS 3000

typedef struct stuttura1{
	vector<double> n_t;
	vector<double> rate_t;
	double Delta_t;
}isto_temps;

typedef struct stuttura2{
	double mu_cps;
	double sigma_cps;
}analisi;

double max(vector<double> t){
    double max=t[0];
    for(int i=0; i<size(t); i++){
        if(t[i]>max){
            max=t[i];
        }
    }
    return max;
}

double min(vector<double> t){
    double min=t[0];
    for(int i=0; i<size(t); i++){
        if(t[i]<min){
            min=t[i];
        }
    }
    return min;
}

isto_temps istogramma_tempi(vector<double> t, int nbins, int len, int numero_figura, char ref_istogramma[30], double numero_acquisizione, char nome_figura[30]){
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
    isto_temps results; 
    char v[300]="Istogramma di ";
    char * title;
    title=strcat(v,ref_istogramma);
    
    TH1F *histogram= new TH1F("histogram",title,nbins, t[0]/3600,(t[len-1]+1)/3600);
    for(int i=0; i<size(t);i++){
        histogram->Fill(t[i]/3600);
    
    }
    
    vector<double> medie_bins;
    for(int i=0;i<nbins;i++){
        medie_bins.push_back(histogram->GetBinCenter(i+1)); //è giusto??? //media dei bins
        results.n_t.push_back(histogram->GetBinContent(i+1));
   
        
    }
    double Deltat = (histogram->GetBinWidth(1))*3600; //in secondi
    histogram->Scale(1/(Deltat));
    
    new TCanvas;
    
    histogram->GetXaxis()->SetTitle("tempo [h]");
    histogram->GetYaxis()->SetTitle("Rate [Hz]");
    histogram->Draw();
    ///////////////////////
    for(int i=0; i<nbins; i++){
        results.rate_t.push_back(results.n_t[i]/Deltat);
    }
    results.Delta_t=Deltat;

    double n_tot = 0;
    for(int i=0;i<nbins;i++){
        n_tot+=results.n_t[i];
    }     

    /////////////////////////
    int maxt= max(t);
    double dt=maxt/nbins;
    double rate = n_tot/maxt;  
    printf("Numero totale conteggi %s è ntot=%lf in Deltat = %lf s ; rate = %lf cps \n", ref_istogramma, n_tot,dt,rate);
    return results;
}

isto_temps istogramma_tc(vector<double> t, int nbins, int len, int numero_figura, char ref_istogramma[30], double numero_acquisizione, char nome_figura[30]){
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
    char v[300]="Istogramma di ";
    char * title;
    title=strcat(v,ref_istogramma);
    
    TH1F *histogram= new TH1F("histogram",title,nbins, t[0]/3600,(t[len-1]+1)/3600);
    for(int i=0; i<size(t);i++){
        histogram->Fill(t[i]/3600);
    
    }
    
    vector<double> medie_bins;
    for(int i=0;i<nbins;i++){
        medie_bins.push_back(histogram->GetBinCenter(i+1)); //è giusto??? //media dei bins
        results.n_t.push_back(histogram->GetBinContent(i+1));
   
        
    }
    double Deltat = (histogram->GetBinWidth(1))*3600; //in secondi
    //histogram->Scale(1/(Deltat));
    
    new TCanvas;
    string dou=to_string((int) Deltat);
    string str1="Occorrenze in bin di ";
    string str2=" s";
    string titlehist=str1+dou+str2;

    histogram->GetXaxis()->SetTitle("tempo [h]");
    histogram->GetYaxis()->SetTitle(titlehist.c_str());
    histogram->Draw();
    ///////////////////////
    for(int i=0; i<nbins; i++){
        results.rate_t.push_back(results.n_t[i]/Deltat);
    }
    results.Delta_t=Deltat;

    double n_tot = 0;
    for(int i=0;i<nbins;i++){
        n_tot+=results.n_t[i];
    }     

    /////////////////////////
    int maxt= max(t);
    double dt=maxt/nbins;
    double rate = n_tot/maxt;  
    printf("Numero totale conteggi %s è ntot=%lf in Deltat = %lf s ; rate = %lf cps \n", ref_istogramma, n_tot,dt,rate);
    return results;
}

analisi analisi_cps(vector<double>  counts, int nbins,double Delta_t, int numero_figura, char ref_istogramma[30], double numero_acquisizione, char nome_figura[30]){
    analisi results;

    printf("\n------Stima di pdf %s-------\n",ref_istogramma);
    
    //stima di media per confrontarla con risultati fit
    double tot=0;
    for(int i=0; i<size(counts);i++){
        tot+=counts[i];
    }
    double mean_exp=tot/size(counts);
    double tot1=0;
    for(int i=0; i<size(counts);i++){
        tot1+=pow((counts[i]-mean_exp),2)/size(counts);
    }

    double std_exp=pow(tot1,0.5);
    
    
    //nota: anche mean e std hanno una loro varianza (sono statistiche)
    double minimo= min(counts);
    double massimo=max(counts);
    printf("Media dei counts = %lf pm %lf ", mean_exp,std_exp);
    string strt1="Istogramma ";
    string strt2=ref_istogramma;
    string strt3=" - acquisizione di ";
    string strt4=to_string((int)numero_acquisizione);
    string strt5=" h";
    string titlevero=strt1+strt2+strt3+strt4+strt5;
    TH1F *histogram= new TH1F("histogram",titlevero.c_str(),nbins,minimo,massimo+1);
    for(int i=0; i<size(counts);i++){ //NBINS O 100????????
        histogram->Fill(counts[i]);
    }
    new TCanvas;
    string str1="Conteggi in #Delta t = ";
    string str2=to_string((int) Delta_t);
    string str3=" s";
    string assex=str1+str2+str3;
    histogram->GetXaxis()->SetTitle(assex.c_str());
    histogram->GetYaxis()->SetTitle("Occorrenze in ogni bin");
    histogram->Draw();

    //preparazione dei dati per il fit: prendo la media dei bins (intera!)
   
    vector<double> medie_bins, occorrenze_rate, x, y, sigma_y;
    int len_x=0;
    vector<int> index;
    for(int i=0;i<nbins;i++){
        
        occorrenze_rate.push_back(histogram->GetBinContent(i+1));
        medie_bins.push_back(histogram->GetBinCenter(i+1));
        if(occorrenze_rate[i]>0){
            x.push_back(medie_bins[i]);
            index.push_back(i);
            len_x++;
            
        
        }
        
    }

    double sum_y=0;//occorrenze totali

    for(int j=0; j<size(x);j++){
        y.push_back(occorrenze_rate[index[j]]);
        sigma_y.push_back(pow(y[j], 0.5)); //errore su y: assumo conteggi poissoniani
        sum_y+=y[j];
        
    }
 
    //fit della poissoniana 
    new TCanvas;

    TGraphErrors *graphic = new TGraphErrors(size(x),&x[0],&y[0],nullptr,&sigma_y[0]);
    
    graphic->SetMarkerStyle(20);
    
    graphic->SetTitle(titlevero.c_str());
    // Define a function 
    //TF1 *f = new TF1("testhi"," TMath::Poisson" , min(x,len_x),max(x,len_x));

    graphic->GetHistogram()->GetXaxis()->SetTitle(assex.c_str());
    graphic->GetHistogram()->GetYaxis()->SetTitle("Occorrenze in ogni bin"); //come metto numero qua?????
     //NON RIESCO A TROVARE poisson
    graphic->Fit("gaus","L");
    TLegend *legend = new TLegend(1.,1.,1.,1.,"Lab. Lesson 1");
    legend->AddEntry(graphic,"Exp. Points","PE");
    //egend->AddEntry(gaus,"Th. Law","L");
    legend->Draw();

    graphic->Draw("APE");
    gStyle->SetOptFit(1111);
    results.mu_cps=0;
    results.sigma_cps=0;
   
    return results;
}

//====================================================================================
// Istogramma delle differenze di tempo tra un evento e l'altro dello stesso telescopio
//====================================================================================
void analisi_delta_t(vector<double> t,int len, int nbins, int numero_figura, char ref_istogramma[30], char nome_figura[30],double numero_acquisizione){
    vector<double> Delta_t;
    Delta_t.push_back(t[1]-t[0]);
    for(int i=1;i<size(t)-1;i++ ){
        Delta_t.push_back(t[i+1]-t[i]);
    }
   


    printf("\n------Stima di %s-------",ref_istogramma);
    double tot=0;
    for(int i=0; i<size(t);i++){
        tot+=Delta_t[i];
    }
    double mean_exp=tot/size(t);
    double tot1=0;
    for(int i=0; i<size(t);i++){
        tot1+=pow((Delta_t[i]-mean_exp),2)/(size(t));
    }
    double std_exp=pow(tot1,0.5);
    double minimo=min(Delta_t);
    double massimo=max(Delta_t);
    printf("Media dei Delta_t = %lf pm %lf ", mean_exp,std_exp);
    string strt1="Istogramma ";
    string strt2=ref_istogramma;
    string strt3=" - acquisizione di ";
    string strt4=to_string((int)numero_acquisizione);
    string strt5=" h";
    string titlevero=strt1+strt2+strt3+strt4+strt5;
    TH1F *histogram= new TH1F("histogram",titlevero.c_str(),nbins, minimo ,massimo+1);
    
    for(int i=0; i<size(Delta_t);i++){ 
        histogram->Fill(Delta_t[i]);
    }
    new TCanvas;
    histogram->Draw();

    
    //stima di media per confrontarla con risultato fit
    vector<double> medie_bins, occorrenze_rate, x, y, sigma_y;
    int len_x=0;
    vector<int> index;
    for(int i=0;i<nbins;i++){
        occorrenze_rate.push_back(histogram->GetBinContent(i+1));
        medie_bins.push_back(histogram->GetBinCenter(i+1));
        if(occorrenze_rate[i]>0){
            x.push_back(medie_bins[i]);
            index.push_back(i);
            len_x++;
            
        
        }
        
    }

    double sum_y=0;//occorrenze totali

    for(int j=0; j<size(x);j++){
        y.push_back(occorrenze_rate[index[j]]);
        sigma_y.push_back(pow(y[j], 0.5)); //errore su y: assumo conteggi poissoniani
        sum_y+=y[j];
        
    }
    

    //fit della exp 
    

    TCanvas* c=new TCanvas(titlevero.c_str());
    c->SetLogy();
    
    TGraphErrors *graphic = new TGraphErrors(size(x),&x[0],&y[0],nullptr,&sigma_y[0]);
    
    graphic->SetMarkerStyle(20);
    
    graphic->SetTitle(titlevero.c_str());
    // Define a function 
    TF1 *f = new TF1("myfunc", "[0]*exp(-x*[1])" , min(x),max(x));
    /*
    plt.ylabel(f'{(bins[4]-bins[3]):.2f} s')
   */
    
    string str1="Occorrenze in bin di ";
    string str2=to_string((int) ((medie_bins[4]-medie_bins[3])*1000));
    string str3=" ms";
    string titleee=str1+str2+str3;
    graphic->GetHistogram()->GetYaxis()->SetTitle(titleee.c_str());
    graphic->GetHistogram()->GetXaxis()->SetTitle("Delta tempo tra eventi [s]"); 
    graphic->Fit("myfunc"); 
    TLegend *legend = new TLegend(1.,1.,1.,1.,"Lab. Lesson 1");
    legend->AddEntry(graphic,"Exp. Points","PE");
    //egend->AddEntry(gaus,"Th. Law","L");
    legend->Draw();

    graphic->Draw("APE");
    gStyle->SetOptFit(1111);
       
}

// =============================================================================
// Istogramma delle differenze di tempo tra i due telescopi
// =============================================================================
void analisi_delta_t12( vector<double> Delta_t, int nbins, int numero_figura, char ref_istogramma[30], char nome_figura[30],double numero_acquisizione){
    
    printf("\n------Stima di %s-------",ref_istogramma);
    double tot=0;
    for(int i=0; i<size(Delta_t);i++){
        tot+=Delta_t[i];
    }
    double mean_exp=tot/size(Delta_t);
    double tot1=0;
    for(int i=0; i<size(Delta_t);i++){
        tot1+=pow((Delta_t[i]-mean_exp),2)/(size(Delta_t));
    }
    double std_exp=pow(tot1,0.5);
    double minimo=min(Delta_t);
    double massimo=max(Delta_t);
    printf("Media dei Delta_t = %lf pm %lf ", mean_exp,std_exp);
    string strt1="Istogramma ";
    string strt2=ref_istogramma;
    string strt3=" - acquisizione di ";
    string strt4=to_string((int)numero_acquisizione);
    string strt5=" h";
    string titlevero=strt1+strt2+strt3+strt4+strt5;
    TH1F *histogram= new TH1F("histogram",titlevero.c_str(),nbins, minimo ,massimo+1);
    for(int i=2; i<size(Delta_t);i++){ 
        histogram->Fill(Delta_t[i]);
    }
    new TCanvas;
    histogram->Draw();

    
    //stima di media per confrontarla con risultato fit
    vector<double> medie_bins, occorrenze_rate, x, y, sigma_y;
    int len_x=0;
    vector<int> index;
    for(int i=0;i<nbins;i++){
        occorrenze_rate.push_back(histogram->GetBinContent(i+1));
        medie_bins.push_back(histogram->GetBinCenter(i+1));
        if(occorrenze_rate[i]>0){
            x.push_back(medie_bins[i]);
            index.push_back(i);
            len_x++;
            
        
        }
        
    }

    double sum_y=0;//occorrenze totali

    for(int j=0; j<size(x);j++){
        if(occorrenze_rate[index[j]]>0){
            y.push_back(occorrenze_rate[index[j]]);
            sigma_y.push_back(pow(y[j], 0.5)); //errore su y: assumo conteggi poissoniani
            sum_y+=y[j];
        }
        
        
    }
    printf("minimo di y= %lf", min(y));

    //fit della exp 
    

    TCanvas* c=new TCanvas(titlevero.c_str());
    c->SetLogy();
    TGraphErrors *graphic = new TGraphErrors(size(x),&x[0],&y[0],nullptr,&sigma_y[0]);
    printf("maxsigmay %lf",max(sigma_y));
    graphic->SetMarkerStyle(20);
    
    
    graphic->SetTitle(titlevero.c_str());
    TF1 *f = new TF1("myfuncion", "[p0]*exp(-x*[p1])" , min(x),max(x));
    f->SetParameters(y[1],1/mean_exp);
    f->SetParNames("N","#lambda");
    gStyle->SetOptFit();
     string str1="Occorrenze in bin di ";
    string str2=to_string((int) ((medie_bins[4]-medie_bins[3])*1000));
    string str3=" ms";
    string titleee=str1+str2+str3;
    graphic->GetHistogram()->GetYaxis()->SetTitle(titleee.c_str());
    graphic->GetHistogram()->GetXaxis()->SetTitle("Delta tempo tra eventi [s]"); 
    
    graphic->Fit("myfuncion"); 

    TLegend *legend = new TLegend(1.,1.,1.,1.,"Lab. Lesson 1");
    legend->AddEntry(graphic,"Exp. Points","PE");
    //egend->AddEntry(gaus,"Th. Law","L");
    legend->Draw();

    graphic->Draw("APE");
    gStyle->SetOptFit(1111);
       
}

void ciao(){

    double Epsilon_tel1= 0.573;
    double N_tel1 = 13.666; //numero di raggi cosmici al secondo attesi nel tel1
    double N_tel2 = 32.13;//numero di raggi cosmici al secondo attesi nel tel2
    double Epsilon_tel2 = 0.246;
 
    int a;
    double b;
    int i;
    int len_t=0;

    FILE * fin=fopen("/home/sarinaleggy/Documenti/lab/Lab_Int/Sciami/Dati/sciami_291122_prova5.dat","r");
    if( fin == NULL ){
        cout << "non riesco ad aprire il primo file" << endl;
        return;
    }

    cout << fin << endl;
    vector<int> ch;
    vector<double> t;

    double time(0.);
    int channel(0);
    int nVar;
    i=0;
    while( (nVar = fscanf(fin,"%d %lf",&channel,&time))!=EOF){
      if (nVar== 2){
        t.push_back(time);
        ch.push_back(channel);
      }

    }
  
    fclose(fin);


    int len_T=0;
    vector<double> t1, t2,tc,T;
    vector<int> ch1, ch2;
    double Delta_t_contatore=0.01;
    double Delta_t_coincidenza=0.0000001;
    T.push_back(t[0]);
    for(i=1; i<size(t);i++){
      if(t[i-1]<=t[i] ){
        T.push_back(t[i]);
        len_T++;
      }
      else if( t[i-1]>t[i] ){
        T.push_back(t[i-1]+Delta_t_contatore);
        len_T++;
        double Delta=t[i-1]-t[i];
        for(int j=i+1;j<size(t);j++){
          t[j]=t[j]+Delta;

        }
      }

    }


    int len_t2=0; 
    int len_t1=0;
    int len_tc=0;
    int len_t12=0;
    vector<double> Delta_t12;
 
    for(int k=0;k< size(T)-1;k++){
        if(ch[k]==1){
            t2.push_back(T[k]);
            len_t2++;
        }
        else if(ch[k]==2){
            t1.push_back(T[k]);
            len_t1++;

        }
        if((ch[k]==1) && (k>8) && (k<size(T)-9)){
            double Delta_tmin =900;
            for(int j=k-8;j<k+9;j++){
                if (ch[j] == 2){
                    if (fabs(T[j]-T[k])<Delta_tmin){
                        Delta_tmin = fabs(T[j]-T[k]);
                        
                    }
                }
            }
            if(Delta_tmin>800){
                printf("\n errore all'indice k= %d \n",k);
            }
            Delta_t12.push_back(Delta_tmin);
            len_t12++;
           

            if(Delta_tmin < Delta_t_coincidenza){
                tc.push_back(T[k]); 
                len_tc++;
                
    
            }
        }
    
    }
 
    
    
    int n_acquis=0;

    int max_t2=* max_element(t2.begin(),t2.end()),max_t1= * max_element(t1.begin(),t1.end());
    
    double durata_acquisizione;
    if(t2[max_t2]>=t1[max_t1]){
        durata_acquisizione = max_t2/3600; //ore
    }
    else if(t2[max_t2]<t1[max_t1]){
        durata_acquisizione = max_t1/3600; //ore
    }
    printf("%lf \n",durata_acquisizione);
    
    //manca tc
    isto_temps parte1, parte2,partetc;
    analisi rate1, rate2;
    printf("size=%lu",size(t1));
    parte1= istogramma_tempi(t1,NBINS,size(t1), n_acquis+1, (char*)"telescopio 1",durata_acquisizione, (char*)"ist_tel1");

    printf("maxt12 %f \n", max(Delta_t12));
    printf("sizeT %lu \n", size(T));

    parte2= istogramma_tempi(t2, NBINS,size(t2),n_acquis+2, (char*)"telescopio 2",durata_acquisizione, (char*)"ist_tel2");
 
    partetc= istogramma_tc(tc, 100,size(tc), n_acquis+3,(char*)"coincidenze", durata_acquisizione,(char*)"ist_coinc");
    
    //Fit dell'istogramma dei conteggi di eventi per unità di tempo per i due canali

    rate1= analisi_cps(parte1.n_t,177, parte1.Delta_t, n_acquis+4, (char*)"conteggi telescopio 1",durata_acquisizione,  (char*)"ist_poiss_tel1");
  
    rate2= analisi_cps(parte2.n_t,177, parte2.Delta_t, n_acquis+5,(char*)"conteggi telescopio 2", durata_acquisizione, (char*)"ist_poiss_tel2");

    //Fit delle differenze di tempo tra un evento e l'altro

    analisi_delta_t(t1,len_t1, 100,n_acquis+6, (char*)"#Delta t telescopio 1", (char*)"ist_dtempi_tel1",durata_acquisizione);
    //per il file prova 5 eliminiamo i primi 5 punti perch e privi di senso

    analisi_delta_t(t2,len_t2,100, n_acquis+7, (char*)"#Delta t telescopio 2",  (char*)"ist_dtempi_tel2", durata_acquisizione);

    analisi_delta_t12(Delta_t12, 100, n_acquis+8,(char*) "#Delta t dei 2 telescopi", (char*)"ist_dtempi_tel12",durata_acquisizione);

}