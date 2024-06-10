#pragma once
#include<iostream>

using namespace std;

//-------------------------------------------
// User-defined parameters
//-------------------------------------------

int addition_number = 5;//number of addition
int addition_interval = 1200; //time interval s
double C_hhtp_add = 0.005; //monomer concentration mol/L
double V_add[] = {0.020000000000000004,0.06,0.1777776,0.1066668,0.0355556}; //volume of addition L



//--------------------------------------------------------
// crystal 
// -------------------------------------------------------

//struct
struct crystal_struct {
	double dia; //nuclei_diameter nm
	double hei; // nuclei_height nm
	double concentration; //concentration mol/L
} crystal[10000000];


//---------------------------------------------------------------------
// The following parameters are code internal calculation parameters
// According to different units
//---------------------------------------------------------------------

//step
//const int max_step = 86400000; //max_step
float i_step = 0; //step

//s
double time_unit = 0.001; //unit_time s
double t;    //time s
double t_in; //induction time  s
int last_add_monomer_time = 0;//The last time of monomer addition s

//L
double V;    //volume L

//mol/L
double c_monomer_hhtp=0; //monomer concentration mol/L
double c_consum_nucleation = 0; //monomer concentration consumption from nucleation mol/L
double c_consum_growth = 0; //monomer concentration consumption from growth mol/L
double c_consum = 0;//monomer concentration consumption mol/L
double c_nucleus = 0;  //total nuclei concentration mol/L
double min_c = 1e-11;//min monomer concentration mol/L

//mol
double monomer_consumption_mol = 0;//total monomer consumption mol
double mol_consum_growth_all = 0; //monomer consumption from nucleation mol
double mol_consum_nucleation_all = 0; //monomer consumption from nucleation mol
double mol_addition = 0; //monomer addition mol

//nm
double ave_d = 0;//average diameter
double max_dia = 0;//largest diameter

//number
int max_step_rale = 0; //final max_step
int i_addition = 0;//number of addition
int point_number = 0; //number of addition = addition_number


//counter(Judgment boundary condition)
int point_end = 0; //=1 end simulation
int point_addition = 1;//=0 monomer addition
int point_it = 1; // =0 induction time
int point_out = 1; //Output data : monomer concentration = 0






