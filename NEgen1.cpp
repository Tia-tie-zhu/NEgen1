#include<iostream>
#include<fstream>
#include<numeric>
#include<algorithm>
#include<iomanip>
#include "cof_parameter.h"
#include "cof_function.h"
#include<queue>
using namespace std;
using std::accumulate;

int main() {


	c_monomer_hhtp = C_hhtp_add ; //monomer concentration mol/L
	V = V_add[i_addition];//volume L
	mol_addition =  C_hhtp_add * V_add[i_addition];//monomer mol
	t = induction_time(kbf1, c_monomer_hhtp);//time s
	i_addition++;//number of additions


	//output
	cout << i_step << "   ";
	cout << t << "   ";
	cout << V << "   ";
	cout << "0   0  "; //average_dia, max_dia
	cout <<  c_nucleus << "   ";
	cout <<  c_monomer_hhtp << "   ";
	cout <<  c_consum_growth << "   ";
	cout <<  c_consum_nucleation << "    ";
	cout <<  mol_consum_growth_all << "   ";
	cout <<  mol_consum_nucleation_all << "    ";
	cout <<  V * c_monomer_hhtp << "     ";
	cout <<  V * C_hhtp_add << endl;
	

	int i = 0;
	while(true)
	{
		
		i_step = i; //step
		t = t + time_unit; //time

		c_consum_growth = 0;//monomer consumption from growth mol/l

		//induction time
		if (point_it == 0) {
			point_it = 1;
			t_in = induction_time(kbf1, c_monomer_hhtp);
			double dia;
			double hei;
			double con;
			//growth
			double it_c_consum_growth = 0; // mol/L
			for (int it_i = 0; it_i < t_in / time_unit + 1; it_i++) {
				for (int j = 0; j < i; j++) {

					dia = crystal[j].dia; //diameter nm

					hei = crystal[j].hei; //height nm

					con = crystal[j].concentration; //nucei cocentration mol/L


					crystal[j].dia = crystal[j].dia + time_unit * growth_rate_lat(kbf1, kbf2, kbf3, c_monomer_hhtp, S, kappa);//diameter nm

					crystal[j].hei = crystal[j].hei + time_unit * growth_rate_ver(kbf1, kbf2, kbf3, c_monomer_hhtp, d); //height nm

					it_c_consum_growth = it_c_consum_growth + PI / 4 * (crystal[j].dia * crystal[j].dia * crystal[j].hei - dia * dia * hei) * beta * con;//monomer cocentration in induction time mol/L
				}

				c_monomer_hhtp = c_monomer_hhtp - it_c_consum_growth;//monomer concentration mol/L
				c_consum_growth = it_c_consum_growth;//monomer consumption from growth mol/l
				mol_consum_growth_all = mol_consum_growth_all + it_c_consum_growth * V;//monomer concentration mol


				//output
				if (it_i % 10000 == 0 || c_monomer_hhtp <  min_c) {

					//diameter
					ave_d = 0.00000;
					max_dia = 0.00000;
					for (int k = 0; k < i + 1; k++)
					{
						ave_d = ave_d + crystal[k].concentration * crystal[k].dia;
						if (crystal[k].dia > max_dia)
							max_dia = crystal[k].dia;
					}


				    //output
					cout << i << "   ";
					cout << t << "   ";
					cout << V << "   ";
					cout <<  ave_d / c_nucleus << "   ";
					cout <<  max_dia << "    ";
					cout <<  c_nucleus << "   ";
					cout <<  c_monomer_hhtp << "   ";
					cout <<  c_consum_growth << "     0      "; //c_consum_nucleation = 0
					cout <<  mol_consum_growth_all << "   ";
					cout <<  mol_consum_nucleation_all << "    ";
					cout <<  V * c_monomer_hhtp << "     ";
					cout << i_addition << endl; //number of additions


				}

				it_c_consum_growth = 0; //nucei cocentration in induction time mol/L

				//c_monomer_hhtp =0
				if (c_monomer_hhtp < min_c) {
					break;
				}
				else {
					i = i + 1;
					t = t + time_unit;

				}
			}

		}


		//growth
		double dia;
		double hei;
		double con;
		for (int j = 0; j < i; j++) {
			
			dia = crystal[j].dia; //diameter nm

			hei = crystal[j].hei; //height nm

			con = crystal[j].concentration; //nuclei concentration mol/L


			crystal[j].dia = crystal[j].dia + time_unit * growth_rate_lat(kbf1, kbf2, kbf3, c_monomer_hhtp, S, kappa);//diameter nm

			crystal[j].hei = crystal[j].hei + time_unit * growth_rate_ver(kbf1, kbf2, kbf3, c_monomer_hhtp, d);//height nm

			c_consum_growth = c_consum_growth + 3.1415926536 / 4 * (crystal[j].dia * crystal[j].dia * crystal[j].hei - dia * dia * hei) * beta * con;//monomer cocentration from growth in induction time mol/L

		}
		mol_consum_growth_all = mol_consum_growth_all + c_consum_growth * V;//monomer cocentration from growth mol


		c_consum_nucleation = nucleation_rate_c_monomer(kbf1, kbf2, kbf3, c_monomer_hhtp, N_nuc_hhtp) * time_unit;//monomer cocentration from nucleation mol/L
		mol_consum_nucleation_all = mol_consum_nucleation_all + c_consum_nucleation * V;//monomer cocentration from growth and nucleation mol

		c_consum = c_consum_growth + c_consum_nucleation;//monomer cocentration from nucleation and growth mol/L


		//Determine whether the monomer concentration is 0
		if (c_monomer_hhtp < min_c) { //monomer concentration = 0

			//interval_time - reaction_time
			t = t + addition_interval - (t - last_add_monomer_time);//s

			//new time_
			double time_ = time_unit * c_monomer_hhtp / (c_consum / time_unit); //s
			c_consum_growth = c_consum_growth * (time_ / time_unit); //monomer cocentration from growth  mol/L
			c_consum_nucleation = c_consum_nucleation * (time_ / time_unit); //monomer cocentration from nucleation mol/L
			monomer_consumption_mol = mol_consum_growth_all + mol_consum_nucleation_all;//monomer cocentration from growth and nucleation mol

			//new nuclei

			crystal[i].dia = d_nuc_init;//diameter nm
			crystal[i].hei = h_nuc_init;//height nm
			crystal[i].concentration = nucleation_rate(kbf1, kbf2, kbf3, c_monomer_hhtp) * time_;//nuclei concentration mol/L
			c_nucleus = c_nucleus + crystal[i].concentration;//nuclei concentration mol/L

			point_addition = 0;
			point_out = 0;

			//c_monomer_hhtp = 0;

			if (point_number == 1) {
				point_end = 1;
			}
				
		}
		else {//monomer concentration != 0
			monomer_consumption_mol = mol_consum_growth_all + mol_consum_nucleation_all;

			//new nuclei
			crystal[i].dia = d_nuc_init;//diameter nm
			crystal[i].hei = h_nuc_init;//height nm
			crystal[i].concentration = nucleation_rate(kbf1, kbf2, kbf3, c_monomer_hhtp) * time_unit;//nuclei concentration mol/L
			c_nucleus = c_nucleus + crystal[i].concentration;//nuclei concentration mol/L

			c_monomer_hhtp = c_monomer_hhtp - c_consum;//monomer concentration mol/L

		}


        //output
		if (i%10000 ==0 || point_out == 0) {
			point_out = 1;

			//diameter
			ave_d = 0.00000;
			max_dia = 0.00000;
			for (int k = 0; k < i + 1; k++)
			{

				ave_d = ave_d + crystal[k].concentration * crystal[k].dia;
				if (crystal[k].dia > max_dia)
					max_dia = crystal[k].dia;
			}
			
			//output
			cout << i << "   ";
			cout << t << "   ";
			cout << V << "   ";
			cout <<  ave_d / c_nucleus << "   ";
			cout <<  max_dia << "    ";
			cout <<  c_nucleus << "   ";
			cout <<  c_monomer_hhtp << "   ";
			cout <<  c_consum_growth << "  "; 
			cout <<  c_consum_nucleation << "  ";
			cout <<  mol_consum_growth_all << "   ";
			cout <<  mol_consum_nucleation_all << "    ";
			cout <<  V * c_monomer_hhtp << "     ";
			cout << i_addition << endl; //number of additions

		}

        //add
		if (i_addition < addition_number) {
			if (point_addition == 0) {
				point_addition = 1;
				//new nuclei concentration
				double new_c_nucleus = 0; //nuclei concentration mol/L
				for (int j = 0; j < i + 1; j++)
				{
					crystal[j].concentration = crystal[j].concentration * V / (V_add[i_addition] + V);//nuclei concentration mol/L
					new_c_nucleus = crystal[j].concentration + new_c_nucleus;//nuclei concentration mol/L
				}
				c_nucleus = new_c_nucleus;//nuclei concentration mol/L

				//new monomer concentration
				c_monomer_hhtp = (c_monomer_hhtp * V + C_hhtp_add * V_add[i_addition]) / (V_add[i_addition] + V);//monomer concentration mol/L
				mol_addition = mol_addition + C_hhtp_add * V_add[i_addition];//monomer  mol
				V = V + V_add[i_addition];//volume L
				last_add_monomer_time = t;//time s
				i_addition = i_addition + 1;//number of addition
			    point_it = 0;
			}

			if (i_addition == addition_number) {
				point_number = 1;
			}
		}


		max_step_rale = i;
		if (point_end == 1)
		{
			break;
		}
			

		i = i + 1;



	}

}





