//-------------------------------------------------------------------------------------------------
//this file includes microscopic parameter definitions and microscopic process function definitions
//-------------------------------------------------------------------------------------------------
#pragma once
#include<iostream>
using namespace std;
#include<cmath>

//Microscopic parameter 
double kbf1 = 0.40667;  //mol/L/s
double kbf2 = 1210.7;  //mol/L/s
double kbf3 = 20361.1;  //mol/L/s
double beta = 0.741640597999647; //number of HHTP units per volume (nm^3) , 1.0/3.897/0.346
double N_nuc_hhtp = 20;//number of monomer for nucleation 
double d =0.346; //the distance between the two layers nm
double S = 1; //disk-shaped crystals  
double kappa = 3.897; //the lateral area divided by the number of vertices in the 2D lattice nm^2
double  d_nuc_init = 7.002431562087647; //nm   (kappa*d*20) = pi*( d_nuc_init/2)^2*h_nuc_init
double  h_nuc_init = 0.7002431562087647; //nm    d_nuc_init = 10*h_nuc_init
double PI = 3.1415926536; //pi
double induction_time(double kbf1, double c_monomer);
double nucleation_rate(double kbf1, double kbf2, double kbf3, double c_monomer);
double nucleation_rate_c_monomer(double kbf1, double kbf2, double kbf3, double c_monomer, double N_nuc_hhtp);
double growth_rate_ver(double kbf1, double kbf2, double kbf3, double c_monomer, double d);
double growth_rate_lat(double kbf1, double kbf2, double kbf3, double c_monomer, double S, double kappa);
double growth_rate_lat_ver_c_monomer(double kbf1, double kbf2, double kbf3, double c_monomer, double c_monomer_nucleus, double dia, double hei, double beta);


double induction_time(double kbf1, double c_monomer)
{
	double induction_time_;
	induction_time_ = 0.068916294722367 / (kbf1 * c_monomer);
	return induction_time_; //s

}


double nucleation_rate(double kbf1, double kbf2, double kbf3, double c_monomer)
{
	double nucleation_rate_;
	double gamma1 = kbf2 / kbf1;
	double gamma2 = kbf3 / kbf1;
    // c_monomer: for HHTP only, in mol/L
    // kbf1, kbf2, kbf3: in L/mol/s
    // gamma1, gamma2: unitless
	nucleation_rate_ = 0.030363169641638444 * c_monomer * c_monomer * kbf1 * pow(gamma1, 0.05) * pow(gamma2, 0.025); // mol/L/s
	return nucleation_rate_;
}

double nucleation_rate_c_monomer(double kbf1, double kbf2, double kbf3, double c_monomer, double N_nuc_hhtp)
{
	double nucleation_rate_c_monomer_; // for HHTP only
    // c_monomer: for HHTP only, in mol/L
    // kbf1, kbf2, kbf3: in L/mol/s
	nucleation_rate_c_monomer_ = nucleation_rate(kbf1, kbf2, kbf3, c_monomer) * N_nuc_hhtp; // mol/L/s
	return nucleation_rate_c_monomer_;
}


double growth_rate_ver(double kbf1, double kbf2, double kbf3, double c_monomer, double d ) //nm/s
{
    double gamma1 = kbf2 / kbf1;
	double gamma2 = kbf3 / kbf1;
	double kbf23 = 2008.0000177867732 - 199.49404467629432 * pow(gamma2,1.0/7)+ pow(pow(gamma1, 3) * pow(gamma2, 1.0/3) / (53248.465539123325 * gamma1 + 13199.729668575232 * gamma2), 0.5);
	double growth_rate_ver_ = d * kbf1 * c_monomer * kbf23;
	return growth_rate_ver_;
}


double growth_rate_lat(double kbf1, double kbf2, double kbf3, double c_monomer, double S, double kappa)//nm/s
{
	double growth_rate_lat_;
	double gamma1 = kbf2 / kbf1;
	double gamma2 = kbf3 / kbf1;
	growth_rate_lat_ = S * kbf1 * c_monomer * pow(pow(gamma1, 2) * (gamma2) * kappa/ (0.8339246759138563 * gamma1 + 1.259398254488272 * gamma2), 0.5);
	return growth_rate_lat_;
}


double growth_rate_lat_ver_c_monomer(double kbf1, double kbf2, double kbf3, double c_monomer, double c_monomer_nucleus, double dia, double hei,double beta, double S, double kappa, double d)
{
	double growth_rate_lat_c_monomer_;
	growth_rate_lat_c_monomer_ = PI /4* beta * (growth_rate_ver(kbf1, kbf2, kbf3, c_monomer, d) * dia * dia + 2 * dia * hei * growth_rate_lat(kbf1, kbf2, kbf3, c_monomer,S,kappa)) * c_monomer_nucleus;
	return growth_rate_lat_c_monomer_;
}




