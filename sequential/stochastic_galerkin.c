#include "stochastic_galerkin.h"

//Compute factorial
int factorial(int n){
    
    int a = 1;

    if ( n == 1 || n == 0 ){
        
        return 1;

    }

    else{

        for (int i = 1; i <= n; i++){
      
            a = a * i;
      
        }
      
    return a;

    }

}

//Compute inner products (C_klm)
double compute_C_klm(int k, int l, int m){
    
    int klm;
    
    int kl_m_2, lm_k_2, mk_l_2;
    
    int lf, mf;
    int kl_m_2_f, lm_k_2_f, mk_l_2_f;
    
    klm = k + l + m;
    
    kl_m_2 = (k + l - m)/2;
    lm_k_2 = (l + m - k)/2;
    mk_l_2 = (m + k - l)/2;
    
    if (klm % 2){
        
        return 0;
        
    }
    
    else{
        
        lf       = factorial(l);
        mf       = factorial(m);
        kl_m_2_f = factorial(kl_m_2);
        lm_k_2_f = factorial(lm_k_2);
        mk_l_2_f = factorial(mk_l_2);
        
        return (lf * mf) / (kl_m_2_f * lm_k_2_f * mk_l_2_f);
        
    }
    
}


//Compute inner products for m = 1 (C_kl1)
double compute_C_kl1(int k, int l){
    
    int kl1;
    int kl_1_2  , l1_k_2  , onek_l_2;
    int kl_1_2_f, l1_k_2_f, onek_l_2_f;
    
    int lf;
    
    kl1 = k + l + 1;
    
    kl_1_2   = (k + l -1)/2;
    l1_k_2   = (l + 1 -k)/2;
    onek_l_2 = (1 + k -l)/2;
    
    //Compute C_kl1
    if (kl1 % 2){
        
        return 0;
        
    }
    
    else{
        
        lf         = factorial(l);
        kl_1_2_f   = factorial(kl_1_2);
        l1_k_2_f   = factorial(l1_k_2);
        onek_l_2_f = factorial(onek_l_2);
        
        return lf / (kl_1_2_f * l1_k_2_f * onek_l_2_f);
        
    }
    
}

