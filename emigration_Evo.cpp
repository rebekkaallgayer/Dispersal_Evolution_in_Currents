#include <iostream>
#include <fstream>  //stream class to both read and write from/to files
#include <string>
#include <math.h> //for atan2 and pow
#include <random> //for distributions
#include <sstream> //for adding strings together

using namespace std;
//variables

int npatches = 40;
int MaxInd = 200;
float InitDisp = 0.1;
int StartPop = 30;
float Lambda = 1.5;
int K = 50;
int b = 1;
float a = (pow(Lambda, 1 / b) - 1) / K;
float current = .5; //between 0 and 1 ( when 0 or 1 all the inds move same direction; when 0.5 equal in both directions)
float mrate = 0.01; //mutation rate
float mortality = 0.01; //cost of dispersal/ mortality rate
float LDD_prob = 0;

int reps = 100;

//if you want to dump outputs at intervals, instead of having it all outputted at the end, write_intervals=true
bool write_intervals = false; 
int interval = 0;
//if you want to track where genotypes originated
bool track_geno = false;
//what kind of boundary do you want: reflective, absorbing, or wrapped, or mixed!
string boundary = "reflective";

//if you want the distribution of genotypes in each patch in gen 1000
bool geno_dist = false;

string working_dir = "~/Documents/current_emigration";


float* create_matrix(int rows, int cols) {
	float *tmp = new float[rows*cols];
	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < cols; j++) {
			tmp[i*cols + j] = -9;
			
		}
	}

	return tmp;
}



void fill_matrix(float* res_mat, float value, int row, int col, int mat_col) { //value=whatever you want to fill your matrix with, row and col however much you want to fill
	for (int r = 0; r < row; r++) {
		for (int c = 0; c < col; c++) {
			res_mat[r*mat_col+c] = value;
		}
	}
}


std::default_random_engine p_generator;
std::default_random_engine u_generator;

void fill_rand_disp(float* res_mat, int row, int col, int mat_col) { //value=whatever you want to fill your matrix with, row and col however much you want to fill
	uniform_real_distribution<float> disp_dist(0.1, 0.9);
	float value;
	
	for (int r = 0; r < row; r++) {
		for (int c = 0; c < col; c++) {
			value= disp_dist(u_generator);
			res_mat[r*mat_col + c] = value;
		}
	}
}

void reproduction(vector<int>& psize, float* disp_mat, int generation, float current, int rep) { //population sizes and dispersal probabilities
	float* temp_disp = create_matrix(npatches, MaxInd); //nrow=patches, col=individuals, temporary dispersal probs
	vector<int> jsize(npatches, 0); //number of juveniles in each patch
	
	for (int p = 0; p < npatches; p++) {//loop through patches
		vector<float> muts; //this is to store the mutations that occur in this patch, if track_geno==true

		float ExpOff = Lambda * pow((1 + a * psize[p]), -b); //calculate expected number of offspring per ind in patch
															//psize(p) is the size of the population in patch p
		for (int i = 0; i < psize[p]; i++) { //loop through inds in patch
			poisson_distribution<int> pois_dist(ExpOff);
			int NOff = pois_dist(p_generator); //how many offspring per individual
			if (NOff > 0) { //if there was offspring,
				for (int joff = 0; joff < NOff; joff++) { //each juvenile
					jsize[p]++; //increase the number of juveniles in jsize
					float inherited_prob= disp_mat[p*MaxInd+i]; //juvenile inherits the disp prob from the parent
					uniform_real_distribution<float> unif_dist1(0.0, 1.0);
					uniform_real_distribution<float> unif_dist2(-0.2, 0.2);

					float new_prob;
					float ran_mrate = unif_dist1(u_generator); //random mutation rate
					if (ran_mrate < mrate) { //if a mutation occurs, ie if mrate=0.01, 1% will get a mutation
						float ran_disp = unif_dist2(u_generator); //random dispersal probability
						new_prob =inherited_prob + ran_disp; //change the dispersal prob of the offspring
						if(new_prob<0){ //if the change has made it negative, just make it 0
							new_prob = 0;
							temp_disp[p*MaxInd+(jsize[p]-1)] = 0;
						}
						else if (new_prob > 1) {
							new_prob = 1;
							temp_disp[p*MaxInd+(jsize[p]-1)] = 1;
						}
						else {
							temp_disp[p*MaxInd + (jsize[p]-1)] = new_prob;
						}
						if (track_geno == true) {
							muts.push_back(new_prob);
						}

						
					}
					else {
						temp_disp[p*MaxInd + (jsize[p] - 1)] = inherited_prob;
					}
				}
			}
		}
		if (track_geno == true) {
			ofstream mut_file;
			stringstream file_name;
			file_name << working_dir<< "/Results/track_geno/current_" << current << "_mut_patch_" << p<<".txt";
			if (generation == 0 && rep==0) {
				mut_file.open(file_name.str().c_str());
			}
			else {
				mut_file.open(file_name.str().c_str(), ios_base::app);
			}
			for (int line = 0; line < muts.size(); line++) {
					mut_file << muts[line] << " " << generation << " " << rep<< " "<< endl;
			}
			mut_file.close();
		}//end of track_geno
		
	}

	psize = jsize; //over-writing pop size vector because there are no overlapping generations! however many juveniles are produced, that is the size of the population for the next discrete timestep
	for (int p = 0; p < npatches; p++){
		for (int i = 0; i < MaxInd; i++) {
			disp_mat[p*MaxInd + i] = temp_disp[p*MaxInd + i];
		}
	}
	
	//delete tempmat
	delete[] temp_disp;
}

void dispersal(float current, vector<int>& psize, float* disp_mat) {
	float* tempmat = create_matrix(npatches, MaxInd); //nrow=patches, col=individuals, temporary dispersal probs
	vector<int> dsize(npatches, 0);

	for (int p = 0; p < npatches; p++) { //loop through patches
		for (int i = 0; i < psize[p]; i++) { //loop through individuals
			int patch_dest=p; //start with current patch until we know if it disperses or not

			uniform_real_distribution<float> unif_dist1(0.0, 1.0);
			float ran_disp = unif_dist1(u_generator);
			float ran_curr= unif_dist1(u_generator);
			float ran_cost = unif_dist1(u_generator);
			float ran_ldd = unif_dist1(u_generator);

			//for rare long distance dispersal, use all individuals, not just dispersers
			if (LDD_prob !=0 && ran_ldd <= LDD_prob) { //if the random number is less than the long distance dispersal prob
					int ldd_patch = rand() % npatches;
					patch_dest = ldd_patch;
					dsize[patch_dest]++; //increase destination patch density by 1
					tempmat[patch_dest*MaxInd + (dsize[patch_dest] - 1)] = disp_mat[p*MaxInd + i]; //carry over dispersal probability (this simply updates the matrix for the individuals that died)

					//note: there is no added risk of mortality here
			}

			//if they dont do LDD:
			else if (disp_mat[p*MaxInd + i] > ran_disp) { //if the individual disperses, if disp prob is .9, 90% of pop will disperse
					//movement between patches
					if (current >= ran_curr) { patch_dest = p + 1; } //downstream ie if current=0.9, 90% of individuals will go downstream
					else { patch_dest = p - 1; } //upstream, 10% will go upstream

					//for reflective boundaries
					if (boundary=="reflective") {
						if (patch_dest < 0) { patch_dest = 0; }
						if (patch_dest > npatches - 1) { patch_dest = npatches - 1; } //it has to be -1 because C++ counts from 0 and i cant do dsize[40], it is out of bounds
					}
					
					//for absorbing boundaries
					else if(boundary=="absorbing"){
						if (patch_dest < 0) { patch_dest = -9; }
						if (patch_dest > npatches - 1) { patch_dest = -9; }
					}

					//for wrapped boundaries
					else if (boundary == "wrapped") {
						if (patch_dest < 0) { patch_dest = npatches-1; } //if it goes further up than the beginning, it reaches the end
						if (patch_dest > npatches - 1) { patch_dest = 0; } //if it goes past the "end" of the stream, it goes back to the beginning
					}

					//for mixed boundaries
					else if (boundary == "mixed") {
						if (patch_dest < 0) { patch_dest = 0; } //simulates the "head" of the river
						if (patch_dest > npatches - 1) { patch_dest = -9; } //if it goes past the "end" of the stream, it is lost to unfavourable habitat
					}

					//cost of dispersal
					if (ran_cost < mortality) { patch_dest = -9; } //ie if mortality=.01 that means 1% die
								//if this happens, individual is simply not added to dsize (dont decrease psize because it has no effect except to mess up the for loop numbers, so the last indivs will be ignored!)

					if (patch_dest != -9) {
						dsize[patch_dest]++; //increase destination patch density by 1
						tempmat[patch_dest*MaxInd + (dsize[patch_dest] - 1)] = disp_mat[p*MaxInd + i]; //carry over dispersal probability (this simply updates the matrix for the individuals that died)
						//because dsize starts from 0, i dont -1 when the indiv moves
						//i dont modify psize because then it gets complicated assigning disp probs
					}
			}
			//if they stay where they are:
			else if (disp_mat[p*MaxInd + i] < ran_disp) {
				dsize[p]++; //add one here if it didnt emigrate because dsize starts from 0!
				tempmat[p*MaxInd+(dsize[p]-1)] = disp_mat[p*MaxInd+i];
			}
		}
	}
	
	psize = dsize; //theyve all moved around but it is still the same gen that just reproduced, so over-write the old population density vector
	for (int i = 0; i < npatches; i++) {
		for (int j = 0; j < MaxInd; j++) {
			disp_mat[i*MaxInd + j] = tempmat[i*MaxInd + j];
		}
	}


	//delete tempmat
	delete[] tempmat;
}



int main() {

	////////////initialise arrays to hold results/////////////
	float* disp_mat= create_matrix(npatches, MaxInd); 

	//make a vector of length npatches, filled with StartPop
	vector<int> psize(npatches, StartPop); 

	float* results=create_matrix(npatches, reps); 

	float* popmatrix = create_matrix(npatches, reps);

	float* varmatrix = create_matrix(npatches, reps);

	ofstream fout;
	ofstream geno_file;
	ofstream pop_file;
	ofstream var_file;

	//this loop can be used to run all currents for multiple values of mutation rate, ldd, dispersal cost, or any other key parameter to investigate
	/*for (int d = 0; d < 2; d++) {
		LDD_prob = 0.006 + d* 0.014;
		cout << "ldd: " << LDD_prob << endl;*/


		for (int ci = 0; ci < 6; ci++) {
			current = 0.5 + ci * 0.1;
			cout << "current: " << current << endl;


			for (int r = 0; r < reps; r++) {
				cout << "rep: " << r << endl;
				// for each rep, work with an empty matrix/vector then initialise it
				//disp probs 
				fill_matrix(disp_mat, -9, npatches, MaxInd, MaxInd); 
				fill_matrix(disp_mat, InitDisp, npatches, StartPop, MaxInd);

				//to start with a range of dispersal probabilities
				//fill_rand_disp(disp_mat, npatches, StartPop, MaxInd);


				//pop sizes in each patch
				for (int i = 0; i < psize.size(); i++) {
					psize[i] = StartPop;
				}

				//for 1000 gens and producing output only at the end
				if (write_intervals == false && track_geno==false) { 

					for (int gen = 0; gen < 1000; gen++) { //let the model run for 1,000 generations
						reproduction(psize, disp_mat, gen, current, r);
						dispersal(current, psize, disp_mat);
					}

					for (int p = 0; p < npatches; p++) {
						//cout << "patch " << p << "(density " << psize[p] << ") :";
						float var = 0; float av; float sd;
						if (!psize[p] == 0) {
							for (int i = 0; i < psize[p]; i++) {
								if(i==0){ results[p*reps + r] = disp_mat[p*MaxInd + i];} //because matrixes are filled with -9, first indiv overwrites placeholder
								else{results[p*reps + r] = results[p*reps + r] + disp_mat[p*MaxInd + i];}	//all other indivs get added to first
							}
							av= results[p*reps + r] / psize[p];//this takes the mean disp probability for each patch for each rep
							results[p*reps + r] = av;

							for (int i = 0; i < psize[p]; i++) {
								var += pow((disp_mat[p*MaxInd + i] - av), 2);
							}
							var /= psize[p];
							sd = sqrt(var);
						}
						else {
							results[p*reps + r] = -9;
						}
						//cout << results[p*reps + r] << endl;
						popmatrix[p*reps + r] = psize[p];
						varmatrix[p*reps + r] = var;
					}

					stringstream ss;
					ss << working_dir<< "/Results/change_curr/mean_disp_prob_c" << current << ".txt";
					fout.open(ss.str().c_str());

					for (int line = 0; line < npatches; line++) {
						for (int col = 0; col < reps; col++) {
							fout << results[line*reps + col] << " ";
						}
						fout << endl;
					}

					fout.close();

					//save population sizes
					//stringstream pp;
					//pp << "../Results/pop_sizes/popsize_c" << current << ".txt";
					//pop_file.open(pp.str().c_str());

					//for (int line = 0; line < npatches; line++) {
					//	for (int col = 0; col < reps; col++) {
					//		pop_file << popmatrix[line*reps + col] << " ";
					//	}
					//	pop_file << endl;
					//}

					//pop_file.close();

					////save variance
					//stringstream vv;
					//vv << "../Results/mean_vars/var_c" << current << ".txt";
					//var_file.open(vv.str().c_str());

					//for (int line = 0; line < npatches; line++) {
					//	for (int col = 0; col < reps; col++) {
					//		var_file << varmatrix[line*reps + col] << " ";
					//	}
					//	var_file << endl;
					//}

					//var_file.close();
				
				} //end of if(write_intervals==false)

				//if you want genotype distribution of gen1000
				if (geno_dist == true) {
					stringstream gg;
					gg << working_dir<< "/Results/track_geno/c" << current << "_rep" << r << ".txt";
					geno_file.open(gg.str().c_str());

					for (int line = 0; line < npatches; line++) {
						for (int col = 0; col < MaxInd; col++) {
							geno_file << disp_mat[line*MaxInd + col] << " ";
						}
						geno_file << endl;
					}

					geno_file.close();

				}//end of if(geno_dist==true)

				//if you're tracking genos, and reps=1!
				if (write_intervals == false && track_geno == true) {
					for (int gen = 0; gen < 1000; gen++) { //let the model run for 1,000 generations
						reproduction(psize, disp_mat, gen, current, r);
						dispersal(current, psize, disp_mat);
					}

					stringstream ss;
					ss << working_dir<< "/Results/track_geno/disp_probs_c" << current <<"rep"<<r<< ".txt";
					fout.open(ss.str().c_str());

					for (int line = 0; line < npatches; line++) {
						for (int col = 0; col < MaxInd; col++) {
							fout << disp_mat[line*MaxInd + col] << " ";
						}
						fout << endl;
					}

					fout.close();
				}//end of track_geno==true

				//for 20,000 gens and outputs every 1000 years
				if (write_intervals == true) { 

					for (int gen = 0; gen < 20001; gen++) { //i need to have 5001 so that i get an output for 5000 years
						reproduction(psize, disp_mat, gen, current, r);
						dispersal(current, psize, disp_mat);

						int gen_counter = 0;
						if (gen == 0 || gen % 1000 == 0) {
							cout << "starting to write files" << endl;
							for (int p = 0; p < npatches; p++) {

								//because i am re-using the same results matrix, i need to clear it everytime i want to fill it
								if (!psize[p] == 0) {
									if (gen != 0) { //if it's not the first time step, reset the results matrix
										results[p*reps + r] = 0;
									}
									for (int i = 0; i < psize[p]; i++) {
										results[p*reps + r] = results[p*reps + r] + disp_mat[p*MaxInd + i]; //fill it with the new values
									}
									results[p*reps + r] = results[p*reps + r] / psize[p];//this takes the mean disp probability for each patch for each rep
								}
								else {
									results[p*reps + r] = 0; //if there are no individuals, there are no values
								}


								if (r == 0) { //if it's the first rep, create the file with the first results
									ofstream fout;
									stringstream ss;
									ss << working_dir<< "/Results/change_mut/mut_" << mrate << "/mean_disp_prob_c" << current << "_gen_" << gen << "rep" << r << ".txt";
									fout.open(ss.str().c_str());
									for (int line = 0; line < npatches; line++) {
										//for (int col = 0; col < reps; col++) {
										fout << results[line*reps + r] << endl;
									}
									fout.close();
									//this should now have one column (for the first rep) in the file
								}

								else {
									ifstream fin; ofstream fout;
									stringstream ii; stringstream oo;
									ii << working_dir<< "/Results/change_mut/mut_" << mrate << "/mean_disp_prob_c" << current << "_gen_" << gen << "rep" << r - 1 << ".txt";
									fin.open(ii.str().c_str());
									oo << working_dir<< "/Results/change_mut/mut_" << mrate << "/mean_disp_prob_c" << current << "_gen_" << gen << "rep" << r << ".txt";
									fout.open(oo.str().c_str());

									string str;
									int line = 0;
									while (getline(fin, str)) {
										fout << str << " " << results[line*reps + r] << "\n";
										line++;
									}

									fin.close(); fout.close();
								}

							} //end of patches
						} //end of writing files
					} //end of gen

				} //end of if(write_intervals==true)
			} //end of reps

			//get rid of all the files except the last rep, so you only end up with one file per generational output. i have to put this extra because i want all the reps to finish before i clear the files
			if (write_intervals == true) { 

				for (int k = 0; k < reps - 1; k++) {
					for (int g = 0; g < 20001; g++) {
						if (g == 0 || g % 1000 == 0) {
							stringstream ii;
							ii << working_dir<<"/Results/change_mut/mut_" << mrate << "/mean_disp_prob_c" << current << "_gen_" << g << "rep" << k << ".txt";
							remove(ii.str().c_str());
						}
					}

				}
				cout << "deleted files for current " << current << endl;
			}
	

		} //end of currents
	//} //end of mrate
		cout << "model is done" << endl;
	cin.get();
	return 0;
}
