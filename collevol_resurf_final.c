
/* Collisional evolution */

/* This code traces the collisional evolution of both volatiles and solids, as described in Bonsor et al, 2023 */

/* This code requires inparam.in input file to have the exact format specified in the example file */

#include <stdio.h>
#include <math.h>
#include<stdlib.h>

double sum_array(double a[], int num_elements)
{
  int i;
  double sum=0;
   for (i=1; i<num_elements; i++)
   {
	 sum = sum + a[i];
   }
   return(sum);
}

double round_to_digits(double value, int digits)
{
    if (value == 0.0) // otherwise it will return 'nan' due to the log10() of zero
        return 0.0;

    double factor = pow(10.0, digits - ceil(log10(fabs(value))));
    return round(value * factor) / factor;   
}

/* Find the average density of a body - this is a variable which changes with time as planetesimals in the kth bin loose volatiles and become denser! */
double average_density(double rho_s, double rho_v, double fv)
{
    double rho_k;
    rho_k=pow(((1-fv)/rho_s + fv/rho_v), -1);
    return rho_k;
}

// find the efficiency at which volatiles are released
double find_chi_tot (double constant_chi, double rho_s, double rho_v, double fvk, double mass_i_k, double mass_i[], int nbin, int ifrag_k, double alpha)
{
    double chi_tot, mass_j, chi_totp=0, chi_f;
    double scrap=0;
    double rho_k=average_density(rho_s, rho_v, fvk);
    double new_const=M_PI*4.*rho_k/(3.*mass_i_k);
    int j;

    for (j=ifrag_k; j<nbin; j++)
        {
            mass_j=mass_i[j];
            chi_f =M_PI*(4. *rho_k/(mass_j))* ( pow((3*mass_j/(4.*rho_k*M_PI)), 2./3.)*constant_chi -pow((3.*mass_j/(4.*rho_k*M_PI)), 1./3.)*pow(constant_chi, 2.) + pow(constant_chi, 3.)/3.);

            if (chi_f<1)
            {
                chi_totp=chi_totp + pow(mass_j, 2-alpha) *chi_f;
            }
            else
            {
                chi_totp=chi_totp + pow(mass_j, 2-alpha);
            }
            scrap= scrap +pow(mass_j,  2- alpha);
        }
    chi_tot=chi_totp/scrap;
   
    if (chi_tot>1)
    {
        chi_tot=1;
    }

    if (pow(6*mass_i[ifrag_k]/(rho_k*M_PI), 1./3.) < constant_chi)
    {
        chi_tot=1;
    }
    return chi_tot;
    
}






int main(void)
{
    double delta_t; // the timestep
    double r_in; // the inner belt radius
    double r_out; // the outer belt radius
    double imax; // the maximum inclination in the disc
    double ecc; // the maximum eccentricity present in the disc
    double au=1.5e11;
    double vrel; // the relative velocity of collisions - assumed to be a constant value
    double delta; // the logarithmic spacing between bins
    double alpha_p; // n(D) dD \propto D^-alpha dD
    double alpha; //n(M) dM \propto M^-alpha_p dM
    double mtot_0; // the total mass in the belt
    double dbl, dmax; // the minimum and maximum size planetesimal in the collisional cascade
    
    double mearth=5.97360e+24;
    double grav=6.67259e-11;
    double msun=1.98900e+30;
    
    double ntime1;
    int ntime;
    int itime, i, imk, nmax, ifv, id, outputinterval;
    char * word;
    int nbin;
    double fv_0; // the initial volatile fraction of planetesimals
    //start from the beginning or read file?
    int start, chi_calc, solids_only, track_origin, resurfacing;  // solids only?
    int time_s=0, time_v=0, time_gas=0;
    int nf=1e2, diam_scrap;
    int print_chi; // determines whether you print out chi or not ..
    int no_cat; // determines whether you include catastrophic collisions or not ..
    print_chi=0; //default value - do not print
    no_cat=0; // if no_cat=1 then no gas is lost in catastrophic collisions, only resurfacing collisions
    float x;
    
    
    double y2, y1, r, y1a, y1b;
    double chi_totp, new_const, scrap, scrap1;
    double fvarray[nf], n_k;
    
    double rho_s=3e3;// kg/m3
    double rho_v=1000; //
    double grain_a=0.01 ; //grain size
    double total_mdot, total_mgas, total_mgas_0; // total mass loss from solids in each time step, total mass in gas or volatiles (should be conserved)
    double constant_chi, chi_k, mass_j; // noted has h in the manuscript - this is the depth from which volatiles are removed in a collision
    
    //Open File
    FILE *ifp, *ofp,*ofpv, *ofgas, *ins, *inv, *ingas, *ofparam, *ofins, *ofinv, *ofer, *ofchi, *forigin, *ofdiameter, *ofratec, *ofrater;
    
    char *mode = "r";
    /* Collision output */
    char outputFilename[] = "collouts.dat", outputFilenamev[] = "colloutv.dat", outputFilenamegas[] = "colloutgas.dat";

    /* Input files, if restarting - see below for instructions regarding format */
    char infilenames[] = "in_s.dat", infilenamev[] = "in_v.dat", infilenamegas[] = "in_gas.dat";
    /* error file - check for error messages! */
    char inparam[]="inparam.in", erfile[]="errors.dat";
    
    /* Storing the initial parameters of the collisional cascade for reference */
    char initials[]="initials.dat", initialv[]="initialv.dat";
    
    /* Storing the collision rate data */
    ofratec = fopen("ratec.dat", "w");
    ofrater = fopen("rater.dat", "w");
    
    ofp = fopen(outputFilename, "w");
    ofpv = fopen(outputFilenamev, "w");
    ofdiameter = fopen("diam.dat", "w");
    ofparam=fopen(inparam, "r");
    ofer=fopen(erfile, "w");
   
    
    if (ofp == NULL || ofpv == NULL || ofratec == NULL || ofrater == NULL || ofparam == NULL|| ofer == NULL)
    {
        printf("Error! Could not open file\n");
        exit(-1); // must include stdlib.h
    }
    
    //read in the parameters from a file named inparam.in [NB check that the format of this file is correct - this is hard coded in (unhelpfully)...]
    
    fscanf(ofparam, "%s\t %d\t %s\n",&word, &start, &word);
    
    if (start == 1){ printf(" Starting from the beginining\n");} else {printf(" starting from input files\n");}
    
    fscanf(ofparam, "%s\t %d\t %s\n",&word, &chi_calc, &word);
    
    fscanf(ofparam, "%s\t %d\t %s\n",&word, &solids_only, &word);
    
    if (solids_only==1)
    {
        printf("Solids Only\n");
        chi_calc=0;
    }
    fscanf(ofparam, "%s\t %d\t %s\n",&word, &track_origin, &word);
    printf("Track_origin: %d\n", track_origin);
    
    if (track_origin==1)
    {
        printf("Tracking origin\n");
    }
    else
    {
        printf("Not tracking origin\n");
    }
    
    if (chi_calc == 1){ printf(" Calculating Chi\n");} else {printf(" Reading chi from input files\n");}
    printf("Made it to chi_cal %d \n", chi_calc);
    
    fscanf(ofparam, "%s\t %lf\t %s\n",&word, &delta_t,&word);
    
    printf("dt %e \n", delta_t);
    
    fscanf(ofparam, "%s\t %lf\t %s\n ",&word, &r_in,&word);
    fscanf(ofparam, "%s\t %lf\t %s\n",&word, &r_out,&word);
    fscanf(ofparam, "%s\t %lf\t %s\n",&word, &imax,&word);
    printf("r_in %e \t rout %e \t imax %e \n", r_in, r_out, imax);
    fscanf(ofparam, "%s\t %lf\n",&word, &ecc);
    printf("ecc %e \n", ecc);
    fscanf(ofparam, "%s\t %lf\t %s\n",&word, &vrel,&word);
    fscanf(ofparam, "%s\t %lf\t %s\n",&word, &delta,&word);
    
    printf("delta %e\t  \n", delta);
    
    
    fscanf(ofparam, "%s\t %lf\n",&word, &alpha_p);
    printf("Word %e\t %e\t  \n", word, alpha_p);
    fscanf(ofparam, "%s\t %lf\t %s\n",&word, &mtot_0,&word);
    fscanf(ofparam, "%s\t %lf\t %s\n",&word, &rho_s,&word);
    fscanf(ofparam, "%s\t %lf\t %s\n",&word, &rho_v,&word);
    fscanf(ofparam, "%s\t %lf\t %s\n",&word, &dmax,&word);
    fscanf(ofparam, "%s\t %lf\t %s\n",&word, &dbl,&word);
    
    fscanf(ofparam, "%s\t %lf\n",&word, &ntime1);
    printf("Ntime %e\n", ntime1);
    fscanf(ofparam, "%s\t %d\n",&word, &outputinterval);
    printf("Outputinterval %d\n", outputinterval);
    printf(" ntime %e, outputinterval %d\n", ntime1, outputinterval);
    printf(" ntime %e, outputinterval %d\n", ntime1, outputinterval);
    
    fscanf(ofparam, "%s\t %lf\n",&word, &fv_0);
    printf("Fv(0) %e\n", fv_0);
    fscanf(ofparam, "%s\t %d\n ",&word, &resurfacing);
    fscanf(ofparam, "%s\t %lf\n ",&word, &constant_chi);
    fscanf(ofparam, "%s\t %d\n ",&word, &print_chi);
    fscanf(ofparam, "%s\t %d\n ",&word, &no_cat);
    printf("No Catastrophic collisions %d\n", no_cat);
    
    printf(" chi %e, resurfacingl %d\n", constant_chi, resurfacing);
    printf("printing chi to file %d\n", print_chi);
    printf("Delta t %e, Rin %e, Rout %e\n", delta_t, r_in, r_out);
    printf("Imax %e, Ecc %e, vrel %e\n", imax, ecc,vrel);
    printf("delta %e, alpha_p %e, mtot0 %e\n", delta, alpha_p, mtot_0);
    printf("rho  s %e, v %e\n",  rho_s, rho_v);
    if (print_chi) {ntime=2;}
    printf("dmax %e\t dbl %e\t ntime %e\t outputinterval %d\n", dmax, dbl, ntime1, outputinterval);
    
    
    //////////////////////////////////////////////////////////////////////
    
    /// find the number of bins to be used based on size of largest planetesimal (dmax) and blow-out size, dbl (nbin)
    i=0;
    x=dmax;
    while (x>dbl)
    {
        x=pow(1-delta, 1./3.)*x;
        i=i+1;
    }
    nbin=i;
    printf(" \n Nbin: %d\t,  %e\t DBL %e\n", nbin, x, dbl);
    
    //////////////////////////////////////////////////////////////////////
    
    /// declare variable arrays which are a function of nbin
    
    double const1,q, rho_k, rmid,mass_i[nbin], mass_s_k[nbin], diam_i[nbin] ;
    /* The catastrophic collision rate, the initial mass in solids in the bin, ..., the total mass in solids in the kth bin, the volatile fraction in the kth bin */
    double rate_c[nbin],   mtot_initial[nbin],mtot_initial1[nbin], mtot_k[nbin], fv[nbin];
    double pik, vol;
    /* the mass gained by the kth bin from catastrophic collision in solids, in volatiles, the gas lost from the kth bin*/
    double mdots_k,  mdotv_k, mdotgas_k,  m_all,  mgas, mvol_change;
    mgas=0;
    double mgas_dust, msolids_dust; // the mass lost directly to dust (volatiles are added to mgas, whilst solids are tracked to test mass conservation)
    double msolids_check, mvol_check; // checking the mass in solids and volatiles at each timestep
    int max_bin_c = 2*nbin; //the maximum bin size to which everything is normalised (i.e. no dust can be produced below this value) set to 1e4
    printf(" \n Nbin: %d\t,  %e\t DBL %e\n", nbin, x, dbl);
    double diam_i_read[nbin],mtotv_k[nbin], rate_r[nbin];
    
    
    
    /* Intrinsic collision probability*/
    rmid=r_in+(r_out-r_in)/2.; /* the belt mid-radius*/
    vrel=ecc*pow((grav*msun)/(rmid*au), 0.5); /*the average velocity of collisions, assumed constant*/
    printf("ecc: %e\n", ecc);
    printf("vrel: %e\n", vrel);
    //////////////////////////////////////////////////////////////////////
    
    /* assign initial array */
    diam_i[1]=dmax;
    
    // assign the average planetesimal mass for each bin
    if (solids_only != 1)
    {
        mass_i[1]=average_density(rho_s, rho_v, fv_0)*pow(diam_i[1],3.)*M_PI/6.;
        mass_s_k[1]=(1-fv_0)*mass_i[1];
    }
    if (solids_only == 1)
    {
        mass_i[1]=rho_s*pow(diam_i[1],3.)*M_PI/6.;
        mass_s_k[1]=mass_i[1];
    }
    printf("Mk : %e\t %e\t %e\t %e \n", mass_i[1], mass_s_k[1], diam_i[1], fv_0);
    
    for (i=2;i<nbin;i++)
    {
        mass_i[i]=mass_i[i-1]*(1-delta);
        if (solids_only == 1)
        {
            mass_s_k[i]=mass_i[i];
            diam_i[i]=pow(mass_i[i]*6./(M_PI*rho_s),1./3.);
        }
        else
        {
            mass_s_k[i]=(1-fv_0)*mass_i[i];
            diam_i[i]=pow(mass_i[i]*6./(M_PI*average_density(rho_s, rho_v, fv_0)) ,1./3.);
        }
        nmax=i;
    }
    
    printf("Mk : %e\t %e\t %e\t %e \n", mass_i[1], mass_s_k[1], diam_i[1], fv_0);
    //////////////////////////////////////////////////////////////////////
    
    /*Primordial size distribution with a constant power law*/
    
    alpha=(alpha_p+2.0l)/3.0l;
    printf("alpha is 11/6 for alpha_p=7/2, standard collisional cascade. Alpha*6 %e \n", alpha*6);
    fprintf(ofp, "%e \t %e \t ",delta_t, alpha);
    for (i=1;i<nbin;i++)
    {
        mtot_initial[i]=pow(mass_i[i],(2.0l-alpha));
    }
    const1=mtot_0*mearth/sum_array(mtot_initial,nbin);
    
    for (i=1;i<nbin;i++)
    {
        mtot_initial1[i]=mtot_initial[i]*const1;
    }
    m_all=sum_array(mtot_initial1, nbin);
    printf("Testing whether the total mass sums to the initial mass input %e %e \n", m_all/mearth, mtot_0);
    
    //// check that there is not more mass in the first bin than the mass of one object!
    if (mtot_initial1[1] < mass_i[1]) { printf("Mass ERROR %e %e\n", mtot_initial1[1], mass_i[1]); }
    //////////////////////////////////////////////////////////////////////
    
    /* QD star and QS star Dispersal and shattering thresholds */
    
    float Qa = 620 ;/*$J kg$^{−1}$, */
    float a = 0.3;
    float Qb = 5.6e-3;/*$J kg$^{−1}$ and $*/
    float b = 1.5 ;
    double qdstar[nbin], qsstar[nbin], mck[nbin], mlr[nbin], chij, chi_tot, chi_tot_all[nbin];
    int irk[nbin], ick[nbin], j, k, p, ifrag[nbin], igravity;
    
    for (i=1;i<nbin;i++)
    {
        qdstar[i]=Qa *pow(diam_i[i],-a) + Qb *pow(diam_i[i],b);
        qsstar[i]=Qa *pow(diam_i[i],-a);
        if (qdstar[i]>1.01*qsstar[i]) {igravity=i;}
        mck[i] = ( 2*qdstar[i]/pow(vrel,2) )*mass_i[i];
        for (j=1; mass_i[j]>mck[i]; j++)
        {
            ick[i]=j;
        }
        if (mass_i[ick[i]]> mass_i[i]) {printf("Error Mck is greater than Mk %d %d \n", ick[i], i);}
        
        if (mass_i[1] <mck[i]) {ick[i]=1;}
        mlr[i] = ( 2*qsstar[i]/pow(vrel,2) )*mass_i[i]; /* finding the mass of hte largest remnant */
        for (j=1; mass_i[j]>mlr[i]; j++){irk[i]=j;}
        if (irk[i]-ick[i] <3)
        {
            irk[i]=ick[i];
        }
        // find the fragment M/2
        for (j=1; mass_i[j]>mass_i[i]/2.; j++){ifrag[i]=j;}
        
    }
 
    printf("Transition from gravity to strength regime %d %e %d %e\n", igravity, diam_i[igravity]);
    /* find an array of chi tot - assume that it independent of fv and use for the rest of the loop! This is a practical measure to speed up computation */
    for (i=1;i<nbin;i++)
    {
        chi_tot_all[i]=find_chi_tot (constant_chi, rho_s, rho_v, fv_0, mass_i[i], mass_i, nbin, ifrag[i], alpha);
    }
    
    //////////////////////////////////////////////////////////////////////
    vol=8*M_PI*pow(rmid*au,3)*ecc*sin(imax)*(1+pow(ecc,2)/3.);
    pik=M_PI*vrel/vol;
    printf("IIII Pik: %e\n", pik);
    fprintf(ofp, "%e \t %e \t ",pik, rmid);
    
    //////////////////////////////////////////////////////////////////////
    
    /* Redistribution Function */
    int kminusi[nbin];
    double const_norm_f_s, redistribution, test_redistribution;
    double f_s[nbin];
    test_redistribution=0;
    const_norm_f_s=0;
 
    for (i=round(log(2)/delta); i<nbin + max_bin_c; i++)
    {
        redistribution= pow((1.0l-delta),(i*(2.0l-alpha)))*delta*(2.0l-alpha)* pow((1./2.),((alpha-2.0l)));
        const_norm_f_s=const_norm_f_s + redistribution;
    }
    printf("Test Redistribution Initial %e\n", const_norm_f_s);
    //Normalize f_s to test_redistribution, such that it always sums to one...
    for (i=1; i<nbin; i++)
    {
        f_s[i]= (1./const_norm_f_s) *pow((1.0-delta),(i*(2.0-alpha)))*delta*(2.0-alpha)* pow((1./2.),((alpha-2.0)));
    }
    test_redistribution=0;
    for (i=round(log(2)/delta); i<nbin; i++)
    {
        test_redistribution=test_redistribution + (1./const_norm_f_s) *pow((1.0l-delta),(i*(2.0l-alpha)))*delta*(2.0l-alpha)* pow((1./2.),((alpha-2.0l)));
    }
    printf("Test Redistribution Initial normalized %e  This is the proportion that is lost to dust: %e\n", test_redistribution, 1.0l-test_redistribution);
    
    /// Do a quick test
    
    for (j=1; j<nbin; j++) // track the mass lost from the jth bin that ends up as dust or volatiles lost directly to dust then gas..
    {
        // the sum of f_s from ilr to nbin was contributed to the bins
        mgas_dust=0;
        if ((nbin-j)> log(2)/delta)
        {
            
            for (i=j+ round(log(2)/delta); i<nbin; i++)
            {
                p=i-j;
                mgas_dust = mgas_dust + (1./const_norm_f_s) *pow((1.0-delta),(p*(2.0-alpha)))*delta*(2.0-alpha)* pow((1./2.),((alpha-2.0)));
                
            }
            
            for (i=nbin; i<nbin+max_bin_c; i++)
            {
                p=i-j;
                mgas_dust = mgas_dust + (1./const_norm_f_s) *pow((1.0-delta),(p*(2.0-alpha)))*delta*(2.0-alpha)* pow((1./2.),((alpha-2.0)));
                
            }
            if (mgas_dust !=1.0) {fprintf(ofer, "Error gas mass not equal to one %d %e %e\n", j, mgas_dust, 1.0-mgas_dust);}
           
        }
    }
 
    
    
    //volatile fraction
    fprintf(ofp, "%e \n ",fv_0);
    
    //////////////////////////////////////////////////////////////////////
    
    //either use initial size distribution, or read from a file !
    if (start==1) /* Starting from t=0 */
    {
        //initialise solids
        if (solids_only !=1)
        {
            for (i=1; i<nbin; i++){mtot_k[i]=mtot_initial1[i]*(1-fv_0);}
        }
        else
        {
            for (i=1; i<nbin; i++){mtot_k[i]=mtot_initial1[i];}
        }
        //initialise volatiles
        total_mgas_0=0;
        if (solids_only !=1)
        {
            for (i=1; i<nbin; i++)
            {
                mtotv_k[i]=mtot_initial1[i]*fv_0;
                fv[i]=mtotv_k[i]/(mtotv_k[i]+mtot_k[i]);
                total_mgas_0=total_mgas_0+ mtotv_k[i];
                //printf("%e\t", fv[i]);
            }
            printf("Total Mass in Volatiles %e\n",total_mgas_0);
        }
        
        printf("Mtot last bin%e\n",mtot_k[nbin-1]);
        
        // open a file to write the initial conditions too...
        ofins = fopen(initials, "w");
        ofinv = fopen(initialv, "w");
        
        if (ofins == NULL || ofinv == NULL)
        {
            printf("Error! Could not open file\n");
            exit(-1); // must include stdlib.h
        }
        
        
        for (i=1; i<nbin; i++)
        {
            fprintf(ofins, "%e \t ",mtot_k[i]);
            if (solids_only !=1) {	fprintf(ofinv, "%e \t ",mtotv_k[i]);}
            if (i==nbin-1) { printf("Mtoti: %e\n",mtot_k[i]);}
            //printf("%d\n",i);
            //	printf("Nbin: %d\n",nbin);
        }
        
        
        printf("Mtot %d\n", i-1);
        
        fclose(ofins);
        if (solids_only !=1) {  fclose(ofinv);}
    }
    else  /* for reading a size distribution from file (in_s.dat, in_v.dat) - restarting the simulations */
    {
        
        if (infilenames == NULL || infilenamev == NULL || infilenamegas == NULL)
        {
            printf("Error! Could not open file\n");
            exit(-1); // must include stdlib.h
        }
        
        printf("Mtot %d\n", i-1);
        
        /// Use these commands to collate the output and restart the simulation, with the files in_s.dat and in_v.dat
        
        /*;
         grep -A 1 "Timestep: 175000" collouts.dat > in_s.dat
          grep -A 1 "Timestep: 175000" colloutv.dat > in_v.dat
          grep -A 2 Gas: collouts.dat > gas.dat
          grep -B 2 "Timestep: 175000" gas.dat > in_gas.dat
          cp collouts.dat collouts_old.dat
          cp colloutv.dat colloutv_old.dat
         Make sure you copy to old!
         
         */
        
        ins = fopen(infilenames, "r");
        inv = fopen(infilenamev, "r");
        ingas = fopen(infilenamegas, "r");
        
        fscanf(ins, "%s\t %d\n",&word, &time_s);
        printf("%c\t %d \n",word, time_s );
        fscanf(inv, "%s\t %d\n", &word, &time_v);
        printf("Starting at S: %d \t V: %d\n", time_s, time_v);
        
        fscanf(ingas, "%s\t %lf\n", &word, &mgas);
        printf("Starting at S: %d \t V: %d \t \n", time_s, time_v);
        
        fprintf(ofp, "\n %d\n   ",time_s);
        fprintf(ofpv, "\n %d\n   ",time_s);
        
        fprintf(ofp, "\n Gas: %e\n   ", mgas);
        
        fclose(ofp);
        ofp = fopen(outputFilename, "w");
        
        for (i=1; i<nbin; i++)
        {
            fscanf(ins, "%lf\t", &mtot_k[i]);
            fscanf(inv, "%lf\t", &mtotv_k[i]);
            fv[i]=mtotv_k[i]/(mtotv_k[i]+mtot_k[i]);
            
            fprintf(ofp, "%e \t ",mtot_k[i]);
            fprintf(ofpv, "%e \t ",mtotv_k[i]);
            
            
            if (i==1) { printf("Mtoti: %e\n",mtot_k[i]);}
        }
        
        total_mgas_0=0;
        if (solids_only !=1)
        {
            for (i=1; i<nbin; i++)
            {
                total_mgas_0=total_mgas_0+ mtotv_k[i];
            }
        }
        
        
        if (time_s !=time_v) {printf("ERROR %d %d\n", time_s, time_v);}
        
    }
    
    ///////////////////////////////////////////////////////////////////////////
    ///
    ///
    
    // Check that Volatile Mass is conserved
    if (solids_only!=1)
    {
        total_mgas=0;
        for (j=1; j<nbin; j++)
        {
            total_mgas=total_mgas + mtotv_k[j];
        }
        total_mgas=total_mgas + mgas;
        if (total_mgas != total_mgas_0)
        {
            printf("Before we start: Whoops Gas mass is not conserved %e \t %e \n", total_mgas, total_mgas_0);
        }
    }

    
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    /***************This is the main loop where mass lost and gained by each bin is calculated at every timestep ********************************************/
    
    
    /// initially no gas m
    if (start==1){mgas=0;}
    total_mdot=0;
    
    printf("About to start timestep loop");
    
    /* sum over each timestep*/
    for (itime=time_s; itime<ntime1; itime++)
    {
        
        test_redistribution=0;
        if (itime % outputinterval == 0 || itime< 10*outputinterval) /* Write to output files every outputinterval only! */
        {
            fprintf(ofp, "\n Gas: %d %e \n Mdot %d %e \n Timestep: %d\n   ",itime, mgas, itime,  total_mdot, itime);
            printf("Mdot: %d %e\n   ", itime, total_mdot);
            
            printf("\n Printing to file: Gas: %e \n Timestep: %d\n   ",mgas, itime);
            fprintf(ofpv, "\n Timestep: %d\n   ",itime);
            fprintf(ofratec, "\n Timestep: %d\n   ",itime);
            fprintf(ofrater, "\n Timestep: %d\n   ",itime);
            fprintf(ofdiameter, "\n Timestep: %d\n   ",itime);
        }
        total_mdot=0;
        mgas_dust=0;
        
        // reiniatilise the arrays!
        for (j=1; j<nbin; j++)
        {
            rate_c[j]=0;
            rate_r[j]=0;
        }
        
        for (j=1; j<nbin; j++)
        {
            /////////////////////////////////////////////////////////////////////////////////////////////////////
            /* Find the density of planetesimals in the kth bin, rho_k and hte number of colliders, n_k*/
            if (solids_only==1)
            {
                rho_k=rho_s;
                n_k=mtot_k[j]/mass_s_k[j];
                if ((diam_i[j] - pow(6*(mtot_k[j])/(M_PI*rho_k*n_k), 1./3.))/diam_i[j]> 1) { printf("ERROR in DIAM %d %e %e\n", j, diam_i[j], pow(6*(mtot_k[j])/(M_PI*rho_k*n_k), 1./3.));}
            }
            else
            {
                rho_k = average_density(rho_s, rho_v, fv_0); // pow(((1-fv_0)/rho_s + fv_0/rho_v), -1)
                n_k=mtot_k[j]/mass_s_k[j];
                if ((diam_i[j] - pow(6*(mtot_k[j]+mtotv_k[j])/(M_PI*rho_k*n_k), 1./3.))/diam_i[j]> 1) { printf("ERROR in DIAM %d %e %e\n", j, diam_i[j], pow(6*(mtot_k[j]+mtotv_k[j])/(M_PI*rho_k*n_k), 1./3.));}
                
                diam_i[j]=pow(6*(mtot_k[j]+mtotv_k[j])/(M_PI*rho_k*n_k), 1./3.);
                
            }
            /////////////////////////////////////////////////////////////////////////////////////////////////////
            
            /// check for problems with diameters
            if (diam_i[j] <0)
            {
                printf("ERROR in DIAM NEGATIVE %e \n", diam_i[j]);
            }
            if (diam_i[j] != diam_i[j])
            {
                printf("ERROR in DIAM  NAN %d %d %e \t %e \t%e \t %e \t %e  %e %e \n", j, itime, diam_i[j], mtot_k[j], mass_s_k[j], mtotv_k[j], fv[j], rho_k, n_k);
            }
            
            if (itime % outputinterval == 0 || itime< 10*outputinterval)
            {
                fprintf(ofdiameter, "%e\t   ",diam_i[j]);
            }
            /*************************************************/
            
            /* Calculate the catastrophic collision rate in the kth bin, rate_c*/
            
            
            for (i=1;i<ick[j]+1; i++)
            {
                rate_c[j]=rate_c[j]+(mtot_k[i]/mass_s_k[i]) *pow((diam_i[i]+diam_i[j]),2.)*pik/4.;
            }

            /* Calculate the resurfacing collision rate in the kth bin, rate_r*/
            if (solids_only!=1)
            {
                for (i=ick[j];i<irk[j]; i++)
                {
                    rate_r[j]=rate_r[j]+(mtot_k[i]/mass_s_k[i])*pow((diam_i[i]+diam_i[j]),2.)*pik;
                }
            }
            ///////////////////////////////
            /// check for problems with rates
            if (rate_c[j] <0)
            {
                rate_c[j]=0;
                printf("ERROR in rate_c NEGATIVE %e %d  %d %d \n", rate_c[j], ick[j], j, itime);
            }
            if (rate_c[j] != rate_c[j])
            {
                printf("ERROR in rate_c NAN %e \t %e \t%e \t %e \t%e  %d %d %d %e %e \n", rate_c[j], mtot_k[j], mtotv_k[j], rho_k, n_k, ick[j], j,itime, mass_s_k[j], diam_i[j]);
                
                for (i=1;i<ick[j]; i++)
                {
                    rate_c[j]=rate_c[j]+(mtot_k[i]/mass_s_k[i]) *pow((diam_i[i]+diam_i[j]),2.)*pik/4.;
                   
                }
                rate_c[j]=0;
            }
            
            if (solids_only!=1)
            {
                if (rate_r[j] <0)
                {
                    rate_r[j]=0;
                    printf("ERROR in rate_r NEGATIVE %e \n", rate_r[j]);
                }
                if (rate_r[j] != rate_r[j])
                {
                    printf("ERROR in rate_r NAN %e \t %e \t%e \t %e \t%e \n", rate_r[j], mtot_k[j], mtotv_k[j], rho_k, n_k);
                    rate_r[j]=0;
                }
            }
            /*************************************************/
        } // j
        
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        
        /*  sum over all bins */
        for (k=1; k<nbin; k++)
        {
            if (solids_only==1)
            {
                rho_k=rho_s;
            }
            else
            {
                rho_k=average_density(rho_s, rho_v, fv_0);//  pow(((1-fv_0)/rho_s + fv_0/rho_v), -1);
            }
            n_k=mtot_k[k]/mass_s_k[k];
            // diameters can change as volatiles lost
            // adjust all diameters before calculations!
            
            /* calculate loss from kth bin*/
            
            // mass gain - based on loss from other bins! (j)  */
            imk=k-log(2)/delta;
            
            /***********************************************************/
            
            mdots_k=0; /* the gain in solids due to catastrophic collisions in the kth bin */
            mdotgas_k=0;
            mdotv_k=0;
            mvol_change=0;
            
            if (imk>0 )
            {
                for (i=1; i<imk; i++)
                {
                    p=k-i;
                    mdots_k=mdots_k+f_s[p]*mtot_k[i]*rate_c[i];
                    
                    if (solids_only!=1)
                    {
                        // find chi_k
                        
                        // label for the bin of M_k/2
                        new_const=0;
      
                        rho_k=average_density(rho_s, rho_v, fv[i]); //pow(((1-fv_0)/rho_s + fv_0/rho_v), -1); //(1-fv[i])*rho_s + fv[i]*rho_v;
                        mass_j=mass_i[i]; // solid mass of objects in the jth bin
                        /* the efficiency at which volatiles are lost*/
                        chi_k=M_PI*(4. *rho_k/(3.*mass_j))* ( 3.*pow((3*mass_j/(4.*rho_k*M_PI)), 2./3.)*constant_chi -3.*pow((3.*mass_j/(4.*rho_k*M_PI)), 1./3.)*pow(constant_chi, 2.) + pow(constant_chi, 3.));
                        
                        if (pow((3*mass_j/(4.*rho_k*M_PI)), 2./3.)*constant_chi + pow(constant_chi, 3.) < pow((3.*mass_j/(4.*rho_k*M_PI)),1./3.)*pow(constant_chi, 2.))
                        {
                            chi_k=1;
                        }
                        if (chi_k>1)
                        {
                            chi_k=1;
                        }
                        
                          if (no_cat==1)
                        {
                            mdotgas_k=mdotgas_k;
                            mdotv_k=mdotv_k+f_s[p]*mtotv_k[i]*rate_c[i];
                            test_redistribution = test_redistribution + f_s[p];
                        }
                        else
                        {
                            mdotv_k=mdotv_k+(1-chi_k)*f_s[p]*mtotv_k[i]*rate_c[i];
                            mdotgas_k=mdotgas_k+chi_k*f_s[p]*mtotv_k[i]*rate_c[i];
                            test_redistribution = test_redistribution + f_s[p];
                        }
                        mvol_change=mvol_change+ f_s[p]*mtotv_k[i]*rate_c[i];
                    } //solids!=1
                } //i ++
            } //imk>0
            if (print_chi & (resurfacing!=1))
            {

                printf("%d, %d, %d, %d, %e, %e, %e, %e, %e \n",k, ifrag[k], nbin,ifrag[k]-nbin, diam_i[k], mtotv_k[k], chi_totp, new_const, chi_k);
            }
            /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            
            /* ADD/LOOSE MASS in SOLIDS to kth bin*/
            /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            
            if (no_cat)
            {
                mtot_k[k]=mtot_k[k];
            }
            else
            {
                mtot_k[k]=mtot_k[k]+(mdots_k-rate_c[k]*mtot_k[k])*delta_t;
            }
            total_mdot=total_mdot + (mdots_k-rate_c[k]*mtot_k[k])*delta_t;
            
            /// ADD/LOOSE MASS IN VOLATILES TO KTH BIN
            if (solids_only!=1)
            {
                if (resurfacing!=1)
                {
                    mtotv_k[k]=mtotv_k[k]+(mdotv_k-rate_c[k]*mtotv_k[k])*delta_t;
                }
                else
                {
                    // find chitot
                    //chi_tot= find_chi_tot (constant_chi, rho_s, rho_v, fv[k], mass_i[k], mass_i, nbin, ifrag[k], alpha);
                    chi_tot=chi_tot_all[k]; // take chi from the array that was calculated at the start (note that this ignores changes to chi from
                    ////////////////////////////////////////////////////
                    
                    if (no_cat)
                    {
                        mtotv_k[k]=mtotv_k[k]-mtotv_k[k]*rate_r[k]*chi_tot*delta_t; // calculate the change in volatile mass due to catastrophic collisions
                    }
                    else
                    {
                        mtotv_k[k]=mtotv_k[k]+(mdotv_k-rate_c[k]*mtotv_k[k] -mtotv_k[k]*rate_r[k]*chi_tot)*delta_t; // calculate the change in volatile mass due to catastrophic and resurfacing collisions
                    }
                    if (print_chi) // print to file the values of chi, in order to create Fig. 4, if the print_chi option is selected
                    {
                        printf("%d, %d, %d, %d, %e, %e, %e, %e, %e \n",k, ifrag[k], nbin,ifrag[k]-nbin, diam_i[k], mtotv_k[k], chi_totp, new_const, chi_tot);
                    }
                }
                
                if (mtotv_k[k] <0){mtotv_k[k]=0;}
                if (mtotv_k[k] != mtotv_k[k]){mtotv_k[k]=0;}
            }//handle Nan
            if (mtot_k <0)
            {
                
                fprintf(ofer, "Error: Solid mass k: %d mtot_k: %e itime: %d\n", k, mtot_k[k], itime);
                mtot_k[k]=0;
            }
            if (mtot_k[k] != mtot_k[k])
            {
                mtot_k[k]=0;
            } //handle Nan
            
            
            // find new volatile fractions
            if (solids_only!=1)
                {
                    fv[k]=mtotv_k[k]/(mtot_k[k]+mtotv_k[k]);
                    if (itime % outputinterval == 0 || itime< 10*outputinterval) // print the total mass in solids + volaitles and collision rate to file
                            {
                                fprintf(ofp, "%e \t ",mtot_k[k]);
                                fprintf(ofratec,"%e \t ",rate_c[k]);
                                fprintf(ofpv, "%e \t ",mtotv_k[k]);
                            }
                }
            else
                {
                    if (itime % outputinterval == 0 || itime< 10*outputinterval) // print the total mass and collision rate to file
                        {
                            fprintf(ofp, "%e \t ",mtot_k[k]);
                            fprintf(ofratec,"%e \t ",rate_c[k]);
                        }
                }
            if (resurfacing!=1 & no_cat!=1)
                    {
                        mgas= mgas +mdotgas_k*delta_t;
                    }
            else
                {
                    if (resurfacing)
                    {
                        if (mtotv_k[k]>0)
                        {
                            mgas= mgas +(mdotgas_k + chi_tot*mtotv_k[k]*rate_r[k])*delta_t; // Calculate the change in gas mass at each timestep
                        }
                        if (itime % outputinterval == 0 || itime< 10*outputinterval)
                        {
                            fprintf(ofrater,"%e \t ",rate_r[k]); // print the resurfacing collision rate to file
                        }
                    }
                    if (no_cat)
                    {
                        mgas= mgas + chi_tot*mtotv_k[k]*rate_r[k]*delta_t;
                    }
                }
            
        } //k
        
        /// Add to gas the mass in volatiles that arrives directly in bins below dbl in size (i.e. straight to dust)
        mgas_dust=0;
        msolids_dust=0;

        for (j=1; j<nbin; j++) // track the mass lost from the jth bin that ends up as dust or volatiles lost directly to dust then gas..
        {
            // the sum of f_s from ilr to nbin was contributed to the bins
            if ((nbin-j)> log(2)/delta)
            {
                for (i=nbin; i<nbin+max_bin_c; i++)
                {
                    p=i-j;
                    mgas_dust = mgas_dust + (1./const_norm_f_s) *pow((1.0-delta),(p*(2.0-alpha)))*delta*(2.0-alpha)* pow((1./2.),((alpha-2.0)))*mtotv_k[j]*rate_c[j];
                    
                }
            }
            
            
        }
                
/* Calculate the total gas mass lost to date*/
            mgas= mgas + mgas_dust*delta_t;

        
        // Check that Volatile Mass is conserved
        if (solids_only!=1)
        {
            if (itime % outputinterval == 0 || itime< 10*outputinterval)
            {
                total_mgas=0;
                for (j=1; j<nbin; j++)
                {
                    total_mgas=total_mgas + mtotv_k[j];
                }
                total_mgas=total_mgas + mgas;
                // Let's fix this by adding any remainder volatile mass that is not tracked to gas
                if (total_mgas != total_mgas_0 && total_mgas_0>total_mgas)
                {
                    
                    mgas =mgas + (total_mgas_0- total_mgas);
                    total_mgas =total_mgas +(total_mgas_0- total_mgas);
                    
                }
            }
        }

	      m_all=sum_array(mtot_k, nbin);

	    } //itime

   fclose(ofp);
   fclose(ofpv);
   fclose(ofgas);
   fclose(ofer);

      return 0;
}
