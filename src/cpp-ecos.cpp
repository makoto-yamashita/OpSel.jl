/*--------------------------------------*/
/*--------------------------------------*/
/*  Tree Computation by SOCP            */
/*--------------------------------------*/
/*--------------------------------------*/

#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <vector>
#include <algorithm>
#include <sys/time.h>

// for SOCP solver, ECOS
#include <ecos.h>

using namespace std;

#define message(msg) {cout << msg << " at line: " << __LINE__ << " file: " << __FILE__ << endl;}
#define EXIT_ERROR (-1)
#define MAX_CHARACTER 1024

void setTimeVal(struct timeval& targetVal)
{
  static struct timezone tz;
  gettimeofday(&targetVal,&tz);
}

double getRealTime(const struct timeval& start,
		   const struct timeval& end)
{
  const long int second = end.tv_sec - start.tv_sec;
  const long int usecond = end.tv_usec - start.tv_usec;
  return ((double)second) + ((double)usecond)*(1.0e-6);
}

#define TimeStart(START__) \
  static struct timeval START__; setTimeVal(START__)
#define TimeEnd(END__) \
  static struct timeval END__; setTimeVal(END__)
#define TimeCal(START__,END__) getRealTime(START__,END__)


void usage(char* argv0)
{
  cout << "Usage   : " << argv0 << " flagLowerbound N_s N inputfile outputfile ecos-wrapper ecos-input ecos-info ecos-x ecos-y ecos-z" << endl;
  cout << "Example : " << argv0 << " 0 14 2800 sP_limited_F1s.csv sP_limited_F1s-result.csv ecos-wrapper.exe ecos-input.txt ecos-info.txt ecos-x.txt ecos-y.txt ecos-z.txt" << endl;
  cout << " +++++++  flagLowerbound -> 0 (lowerbound is 0) or  1 (lowerbound is the last column of inputfile)" << endl;
}


class oneLine
{
public:
  idxint id;
  idxint mother;
  idxint father;
  idxint original_id;
  idxint original_place;
  idxint original_mother;
  idxint original_father;
  double coefficient;
  double lowerbound;
  double upperbound;
  double contribution;

  idxint print() {
    printf("id = %ld, mother = %ld, father = %ld, coeffcient = %e, lb = %e, ub = %e\n",
	   id, mother, father, coefficient, lowerbound, upperbound);
    printf("orig_id = %ld, orig_place = %ld, orig_mother = %ld, orig_father = %ld\n",
	   original_id, original_place, original_mother, original_father);
    printf("contribution = %e\n", contribution);

  }

  // for sorting by original_id
  static bool compare(oneLine* left, oneLine* right) {
    return (left->original_id) < (right->original_id);
  }

  // for sorting by original_place (to output the final solution)
  static bool compare_back(oneLine* left, oneLine* right) {
    return (left->original_place) < (right->original_place);
  }

};

// To construct Ainv in the sparse style
class j_element
{
public:
  idxint i;
  pfloat value;
  // constructor
  j_element(idxint i, pfloat value) {
    this->i     = i;
    this->value = value;
  }
  void print() {
    printf("i=%ld, value=%e\n",i,value);
  }
  static bool compare(j_element* left, j_element* right) {
    return (left->i) < (right->i);
  }
};

// y = Ainv*x (Ainv must be [n times n], x must be longer than n [first n elements are used])
void multi_Ainv_x(idxint n, idxint* Ainv_jc, idxint* Ainv_ir, pfloat* Ainv_pr, 
		  pfloat* x, pfloat* y) {
  for (idxint index1=0; index1<n; ++index1) {
    y[index1] = 0.0;
  }

  for (idxint j=0; j<n; ++j) {
    pfloat x_j = x[j];
    idxint Ainv_j_start = Ainv_jc[j];
    idxint Ainv_j_end   = Ainv_jc[j+1];
    // printf("Ainv[%ld] -> [%ld:%d] :: x[%ld] = %e\n", j, Ainv_j_start, Ainv_j_end, j, x_j);
    for (idxint Ainv_index = Ainv_j_start; Ainv_index < Ainv_j_end; ++Ainv_index) {
      idxint i = Ainv_ir[Ainv_index];
      y[i] += Ainv_pr[Ainv_index] * x_j;
    }
  }
}


int main(int argc, char** argv)
{
  TimeStart(MAIN_START);
  if (argc != 12) {
    cout << "---------- Error!! Usage is not correct -------------" << endl;
    cout << "---  argc is now " << argc << ", but it should be 12 --" << endl;
    usage(argv[0]);
    exit(EXIT_ERROR);
  }

  int argIndex = 0;
  
  argIndex++;
  int flagLowerbound = -1; // dummy initialize
  if (sscanf(argv[argIndex], "%d", &flagLowerbound) == 0) {
    cout << "Error!! cannot read integer flagLowerbound from " << argv[argIndex] << endl;
    usage(argv[0]);
    exit(EXIT_ERROR);
  }
  if (flagLowerbound != 0 && flagLowerbound != 1) {
    cout << "Error!! flagLowerbound should be 0 or 1" << endl;
    usage(argv[0]);
    exit(EXIT_ERROR);
  }

  argIndex++;
  double N_s = 0; // dummy initialize
  if (sscanf(argv[argIndex], "%lf", &N_s) == 0) {
    cout << "Error!! cannot read floating-number N_s from " << argv[argIndex] << endl;
    usage(argv[0]);
    exit(EXIT_ERROR);
  }
  double theta = 0.5/N_s;

  argIndex++;
  int Nint = 0; // dummy initialize
  if (sscanf(argv[argIndex], "%d", &Nint) == 0) {
    cout << "Error!! cannot read integer N from " << argv[argIndex] << endl;
    usage(argv[0]);
    exit(EXIT_ERROR);
  }
  idxint N = Nint;
  printf("flagLowerbound = %d, N_s = %lf, theta = %e, N = %d\n", flagLowerbound, N_s, theta, N);

  argIndex++;
  char* filename_input = argv[argIndex];
  FILE* fp_input = fopen(filename_input, "r");
  if (fp_input == NULL) {
    cout << "Error!! cannot open " << filename_input << endl;
    exit(EXIT_ERROR);
  }
  argIndex++;
  char* filename_output = argv[argIndex];
  FILE* fp_output = fopen(filename_output, "w");
  if (fp_output == NULL) {
    cout << "Error!! cannot open " << filename_output << endl;
    exit(EXIT_ERROR);
  }
  argIndex++; char* ecos_wrapper = argv[argIndex];
  argIndex++; char* ecos_input   = argv[argIndex];
  argIndex++; char* ecos_info    = argv[argIndex];
  argIndex++; char* ecos_x       = argv[argIndex];
  argIndex++; char* ecos_y       = argv[argIndex];
  argIndex++; char* ecos_z       = argv[argIndex];

  printf("All arguments are analyzed (argIndex = %d)\n", argIndex);

  printf("Reading %s ...\n", filename_input);
  vector<oneLine*> readLines;
  idxint readLineCount = 0;
  idxint readItemCount = 0;
  char readline[MAX_CHARACTER];
  while (fp_input != NULL) {
    if (fgets(readline, MAX_CHARACTER, fp_input) == NULL) {
      printf("Reached End of File, %ld lines  \n", readLineCount);
      break;
    }
    readLineCount++;
    if (readline[0] < '0' || readline[0] > '9') {
      #if 0
	printf("Line %ld is comment, skip\n", readLineCount);
      #endif
      continue;
    }
    #if 0
      printf("Read [%ld] => %s\n", readLineCount, readline);
    #endif
    oneLine* oneline = new oneLine;
    oneline->original_place = readItemCount;
    oneline->id     = -1; // dummy initialize;
    oneline->mother = -1; //dummy initialize;
    oneline->father = -1; //dummy initialize;
    oneline->lowerbound   = 0.0;
    oneline->contribution = 0.0;
    idxint tokenCount = 0;
    char* token = NULL;
    token = strtok(readline,",");
    if (token != NULL) {
      oneline->original_id = atoi(token);
      tokenCount++;
    }
    token = strtok(NULL,",");
    if (token != NULL) {
      oneline->original_mother = atoi(token);
      tokenCount++;
    }
    token = strtok(NULL,",");
    if (token != NULL) {
      oneline->original_father = atoi(token);
      tokenCount++;
    }
    token = strtok(NULL,",");
    if (token != NULL) {
      oneline->coefficient = atof(token);
      tokenCount++;
    }
    token = strtok(NULL,",");
    if (token != NULL) {
      oneline->upperbound = atof(token);
      tokenCount++;
    }
    if (flagLowerbound == 1) {
      token = strtok(NULL,",");
      if (token != NULL) {
	oneline->lowerbound = atof(token);
	#if 0
	if (oneline->lowerbound != 0.0) {
	  printf("Found !! Line = %d, lb = %e \n",
		 readLineCount, oneline->lowerbound);
	}
	#endif
	tokenCount++;
      }
    }
      
    if (flagLowerbound == 0 && tokenCount != 5) {
      printf("Skip Line %ld, anayalized faled, tokenCount = %ld (should be 5) \n", readLineCount, tokenCount);
      delete oneline;
      continue;
    }
    if (flagLowerbound == 1 && tokenCount != 6) {
      printf("Skip Line %ld, anayalized faled, tokenCount = %ld (should be 6) \n", readLineCount, tokenCount);
      delete oneline;
      continue;
    }

    #if 0
      printf("Readline %ld as --> \n", readLineCount);
      oneline->print();
    #endif

    readLines.push_back(oneline);
    readItemCount++;
    
  } // end of "while (fp_input != NULL)"
    
  fclose(fp_input);
  printf("Item Number = %ld\n", readLines.size());
  TimeStart(CONVERT_START);

  sort(readLines.begin(), readLines.end(), oneLine::compare);
  const idxint Z = readLines.size();

  // Assign the sorted id
  for (idxint index1=0; index1<Z; ++index1) {
    readLines[index1]->id = index1;
  }
  idxint p_mother    = -2;
  idxint p_mother_id = -2;
  idxint p_father    = -2;
  idxint p_father_id = -2;
  
  // Assign the sorted mother and father by binary search
  for (idxint index1=0; index1<Z; ++index1) {
    oneLine* oneline = readLines[index1];
    idxint mother = oneline->original_mother;
    idxint father = oneline->original_father;

    // mother must be larger than or equal to father
    // (more precisely, if mother == 0, then father == 0
    if (mother < father) {
      idxint tmp1 = father;
      father = mother;
      mother = tmp1;
    }

    if (mother == p_mother) {
      oneline->mother = p_mother_id;
    }
    else if (mother == 0) { // no mother assigned
      oneline->mother = -1;
      p_mother = mother;
      p_mother_id = -1;
    }
    else { // binary search
      p_mother = mother;
      idxint lower = 0;
      idxint upper = index1-1;
      idxint id = mother;
      while (1) {
	idxint mid = (lower+upper)/2;
	idxint mid_id = readLines[mid]->original_id;
	if (mid_id > id) {
	  upper = mid-1;
	}
	else if (mid_id < id) {
	  lower = mid+1;
	}
	else { // (mid_id == id)
	  oneline->mother = mid;
	  p_mother_id = mid;
	  break;
	}
      }
    }


    if (father == p_father) {
      oneline->father = p_father_id;
    }
    if (father == 0) { // no mother assigned
      oneline->father = -1;
      p_father = father;
      p_father_id = -1;
    }
    else { // binary search 
      p_father = father;
      idxint lower = 0;
      idxint upper = index1-1;
      idxint id = father;
      while (1) {
	idxint mid = (lower+upper)/2;
	idxint mid_id = readLines[mid]->original_id;
	if (mid_id > id)  {
	  upper = mid-1;
	}
	else if (mid_id < id) {
	  lower = mid+1;
	}
	else { // (mid_id == id) 
	  oneline->father = mid;
	  p_father_id = mid;
	  break;
	}
      }
    }
  }  // End of   Assign the sorted mother and father by binary search
 
  printf("Data consitency has been checked.\n");
  
  #if 0
  for (idxint index1=0; index1<Z; ++index1) {
    readLines[index1]->print();
  }
  #endif

  printf("Computing inbreeding coefficients ...\n");

  double* inbreeding = new double[Z];

  printf("This part is not available now.");

#if 0
  printf("inbreeding = \n");
  for (idxint index1 = 0; index1 < Z; ++index1) {
    printf("%e ", inbreeding[index1]);
  }
  printf("\n");
  idxint nnz_inbreeding = 0;
  for (idxint index1 = 0; index1 < Z; ++index1) {
    if (inbreeding[index1] != 0) {
      printf("inbreeding[%ld] = %e\n", index1, inbreeding[index1]);
      nnz_inbreeding++;
    }
  }
  printf("nnz = %d\n", nnz_inbreeding);
#endif

  idxint Ainv_nnz   = 0;
  idxint Ainv_index = 0;

  idxint* Ainv_jc = NULL;
  idxint* Ainv_ir = NULL;
  pfloat* Ainv_pr = NULL;

  TimeStart(BEFORE_AINV_SPARSE);
  printf("Constructing Ainv in the sparse format\n");

  vector<j_element*>* A_add;
  A_add = new vector<j_element*>[Z];
  for (idxint index1=0; index1<Z; ++index1) {
    oneLine* oneline = readLines[index1];
    idxint mother = oneline->mother;
    idxint father = oneline->father;
    
    if (mother >= 0 && father >= 0) {
      const double b = 4.0 / ((1-inbreeding[mother]) + (1-inbreeding[father]));
      /*
      Ainv_full[index1 + index1 * Z] += 1.0*b;
      Ainv_full[mother + index1 * Z] -= 0.5*b;
      Ainv_full[index1 + mother * Z] -= 0.5*b;
      Ainv_full[father + index1 * Z] -= 0.5*b;
      Ainv_full[index1 + father * Z] -= 0.5*b;
      Ainv_full[mother + mother * Z] += 0.25*b;
      Ainv_full[mother + father * Z] += 0.25*b;
      Ainv_full[father + mother * Z] += 0.25*b;
      Ainv_full[father + father * Z] += 0.25*b;
      */
      j_element* oneelement;
      oneelement = new j_element(index1,+1.0*b);
      A_add[index1].push_back(oneelement);
      oneelement = new j_element(mother,-0.5*b);
      A_add[index1].push_back(oneelement);
      oneelement = new j_element(index1,-0.5*b);
      A_add[mother].push_back(oneelement);
      oneelement = new j_element(father,-0.5*b);
      A_add[index1].push_back(oneelement);
      oneelement = new j_element(index1,-0.5*b);
      A_add[father].push_back(oneelement);
      oneelement = new j_element(mother,+0.25*b);
      A_add[mother].push_back(oneelement);
      oneelement = new j_element(mother,+0.25*b);
      A_add[father].push_back(oneelement);
      oneelement = new j_element(father,+0.25*b);
      A_add[mother].push_back(oneelement);
      oneelement = new j_element(father,+0.25*b);
      A_add[father].push_back(oneelement);
    }
    else if (mother >= 0 && father < 0) {
      // since mother >= father, 
      // the case (mother < 0 && father >= 0) cannot happen
      const double b = 4.0 / (1.0*(1-inbreeding[mother]) + 2.0*(1-0));
      /*
      Ainv_full[index1 + index1 * Z] += 1.0*b;
      Ainv_full[mother + index1 * Z] -= 0.5*b;
      Ainv_full[index1 + mother * Z] -= 0.5*b;
      Ainv_full[mother + mother * Z] += 0.25*b;
      */
      j_element* oneelement;
      oneelement = new j_element(index1,+1.0*b);
      A_add[index1].push_back(oneelement);
      oneelement = new j_element(mother,-0.5*b);
      A_add[index1].push_back(oneelement);
      oneelement = new j_element(index1,-0.5*b);
      A_add[mother].push_back(oneelement);
      oneelement = new j_element(mother,+0.25*b);
      A_add[mother].push_back(oneelement);
    }
    else {
      // the case (mother < 0 && father < 0)
      const double b = 4.0 / (2.0*(1-0) + 2.0*(1-0));
      /*
      Ainv_full[index1 + index1 * Z] += 1.0*b;
      */
      j_element* oneelement;
      oneelement = new j_element(index1,+1.0*b);
      A_add[index1].push_back(oneelement);
    }
  } 

  for (idxint index1=0; index1<Z; ++index1) {
    sort(A_add[index1].begin(), A_add[index1].end(), j_element::compare);
  }
  #if 0
  for (idxint index1=0; index1<Z; ++index1) {
    vector<j_element*>& A_j = A_add[index1];
    idxint length = A_j.size();
    for (idxint index2=0; index2<length;++index2) {
      j_element* oneelement=A_j[index2];
      printf("A_add[%ld]->i = %ld, value=%e\n",
	     index1, oneelement->i, oneelement->value);
    }
  }
  #endif

  idxint* j_diff_count = new idxint[Z];
  for (idxint index1=0; index1<Z; ++index1) {
    j_diff_count[index1] = 0;
  }

  for (idxint index1=0; index1<Z; ++index1) {
    idxint count1 = 1;
    vector<j_element*>& A_j = A_add[index1];
    idxint length = A_j.size();
    for (idxint index2=0; index2<length-1;++index2) {
      if (A_j[index2]->i != A_j[index2+1]->i) {
	count1++;
      }
    }
    j_diff_count[index1] = count1;
  }

  idxint* Ainv_sparse_jc = new idxint[Z+1];
  Ainv_sparse_jc[0] = 0;
  for (idxint index1=0; index1<Z; ++index1) {
    Ainv_sparse_jc[index1+1] = Ainv_sparse_jc[index1] + j_diff_count[index1];
  }

  Ainv_nnz = Ainv_sparse_jc[Z];
  idxint* Ainv_sparse_ir = new idxint[Ainv_nnz];
  pfloat* Ainv_sparse_pr = new pfloat[Ainv_nnz];
  idxint Ainv_sparse_index = 0;
  for (idxint index1=0; index1<Z; ++index1) {
    vector<j_element*>& A_j = A_add[index1];
    idxint length = A_j.size();
    idxint index2 = 0;
    while (index2 < length) {
      idxint add_i = A_j[index2]->i;
      double add_value = A_j[index2]->value;
      idxint index3 = index2+1;
      while (index3 < length && A_j[index3]->i == add_i) {
	add_value += A_j[index3]->value;
	index3++;
      }
      //  printf("Ainv_sparse_index = %ld, Ainv_nnz = %ld, j=%ld\n",
      //	 Ainv_sparse_index, Ainv_nnz, index1);
      Ainv_sparse_ir[Ainv_sparse_index] = add_i;
      Ainv_sparse_pr[Ainv_sparse_index] = add_value;
      Ainv_sparse_index++;
      index2 = index3;
    }
  }


  TimeEnd(AFTER_AINV_SPARSE);
  printf("Constructing Ainv is %.3lf seconds\n",
	 TimeCal(BEFORE_AINV_SPARSE, AFTER_AINV_SPARSE));
  Ainv_jc =   Ainv_sparse_jc;
  Ainv_ir =   Ainv_sparse_ir;
  Ainv_pr =   Ainv_sparse_pr;
  #if 0 
  printf("Ainv_nnz = %ld", Ainv_nnz);
  #endif

  #if 0
  for (idxint j=0; j<Z; ++j) {
    idxint Ainv_j_start = Ainv_jc[j];
    idxint Ainv_j_end   = Ainv_jc[j+1];
    for (idxint Ainv_index = Ainv_j_start; Ainv_index < Ainv_j_end; ++Ainv_index) {
      idxint i = Ainv_ir[Ainv_index];
      printf("A[%ld,%ld] = %e\n", i,j, Ainv_pr[Ainv_index]);
    }
  }
  #endif

  printf("Setting SOCP data ...\n");
  idxint ecos_n      = Z+1;
  idxint ecos_m      = 2+Z+Z+2+(1+Z);
  idxint ecos_p      = 0;
  idxint ecos_l      = 2+Z+Z+2;
  idxint ecos_ncones = 1;
  idxint ecos_q[1]   = {(idxint) (1+Z)};
  
  pfloat* ecos_c = NULL;
  ecos_c = new pfloat[Z+1];
  pfloat* minus_c = new pfloat[Z+1];
  for (idxint index1 = 0; index1 < Z; ++index1) {
    // printf("coefficient[%d] = %e\n", index1, readLines[index1]->coefficient);
    minus_c[index1] = -(readLines[index1]->coefficient);
  }
  // ecos_c[0:Z-1] = Ainv*minus_c[Z-1];
  multi_Ainv_x(Z, Ainv_jc, Ainv_ir, Ainv_pr, minus_c, ecos_c);
  ecos_c[Z] = 0; // for t

  pfloat* ecos_h = NULL;
  ecos_h = new pfloat[1 + 1 + Z+Z+1+1+1+Z];
  // ecos_h = [0; 0; u/N; -l/N; 1; -1; sqrt(2*theta); zeros(Z,1)];
  ecos_h[0] = 0;
  ecos_h[1] = 0;
  for (idxint index1=0; index1<Z; ++index1) {
    ecos_h[2+index1] = (readLines[index1]->upperbound)/N;
  }
  for (idxint index1=0; index1<Z; ++index1) {
    ecos_h[2+Z+index1] = -(readLines[index1]->lowerbound)/N;
  }
  ecos_h[2+Z+Z+0] = 1;
  ecos_h[2+Z+Z+1] = -1;
  ecos_h[2+Z+Z+2] = sqrt(2.0*theta);
  for (idxint index1=0; index1<Z; ++index1) {
    ecos_h[2+Z+Z+3+index1] = 0.0;
  }

  pfloat* e_vector = new pfloat[Z];
  for (idxint index1 = 0; index1 < Z; ++index1) {
    e_vector[index1] = 1.0;
  }
  // Ainv_e = Ainv*e_vector
  double* Ainv_e = new pfloat[Z];
  multi_Ainv_x(Z, Ainv_jc, Ainv_ir, Ainv_pr, e_vector, Ainv_e);

  // ecos_G = [Ainv; -Ainv; (Ainv*e)'; -(Ainv*e)'; zeros(1,Z); SOCP part]
  idxint* ecos_G_jc = NULL;
  idxint* ecos_G_ir = NULL;
  pfloat* ecos_G_pr = NULL;

  idxint* SOCP_counter = new idxint[Z];
  for (idxint j=0; j<Z; ++j) {
    SOCP_counter[j] = 0;
  }
  for (idxint j=0; j<Z; ++j) {
    oneLine* oneline = readLines[j];
    idxint id = oneline->id;
    idxint mother = oneline->mother;
    idxint father = oneline->father;
    SOCP_counter[id]++;
    if (mother >= 0) {
      SOCP_counter[mother]++;
    } 
    if (father >= 0 && father != mother) {
      SOCP_counter[father]++;
    } 
  }

  ecos_G_jc = new idxint[Z+1+1];
  ecos_G_jc[0] = 0;
  int zeroAinve = 0;
  for (idxint j=0; j<Z; ++j) {
    ecos_G_jc[j+1] = ecos_G_jc[j] + 2*(Ainv_jc[j+1] - Ainv_jc[j]) + SOCP_counter[j];
    if (Ainv_e[j] != 0) {
      zeroAinve++;
      ecos_G_jc[j+1] += 2;
    }
  }
  ecos_G_jc[Z+1] = ecos_G_jc[Z] + 2 + Z;
  #if 0
  for (idxint j=0; j<Z+1; ++j) {
    printf("A_[%ld] = %d\n", j+1, Ainv_jc[j+1]);
    printf("G_[%ld] = %d\n", j+1, ecos_G_jc[j+1]);
  }
  #endif
  idxint ecos_G_nnz = ecos_G_jc[Z+1];
  // printf("ecos_G_nnz = %ld\n", ecos_G_nnz);
  ecos_G_ir = new idxint[ecos_G_nnz];
  ecos_G_pr = new pfloat[ecos_G_nnz];

  for (idxint j=0; j<Z; ++j) {
    // Ainv
    idxint Ainv_j_start = Ainv_jc[j];
    idxint Ainv_j_end   = Ainv_jc[j+1];
    idxint ecos_G_index = ecos_G_jc[j];
    for (idxint Ainv_index = Ainv_j_start; Ainv_index < Ainv_j_end; ++Ainv_index) {
      ecos_G_ir[ecos_G_index] = 2+Ainv_ir[Ainv_index];
      ecos_G_pr[ecos_G_index] = Ainv_pr[Ainv_index];
      ecos_G_index++;
    }
    // -Ainv
    for (idxint Ainv_index = Ainv_j_start; Ainv_index < Ainv_j_end; ++Ainv_index) {
      ecos_G_ir[ecos_G_index] = 2+Z+Ainv_ir[Ainv_index];
      ecos_G_pr[ecos_G_index] = -Ainv_pr[Ainv_index];
      ecos_G_index++;
    }
  }
  for (idxint j=0; j<Z; ++j) {
    // Ainv*e, -Ainv*e
    if (Ainv_e[j] != 0) {
      idxint j_point = ecos_G_jc[j] + 2*(Ainv_jc[j+1] - Ainv_jc[j]) ;
      ecos_G_ir[j_point  ] = 2+2*Z;
      ecos_G_pr[j_point  ] = Ainv_e[j];
      ecos_G_ir[j_point+1] = 2+2*Z+1;
      ecos_G_pr[j_point+1] = -Ainv_e[j];
    }
  }
  idxint* j_counter = new idxint[Z];
  for (idxint j=0; j<Z; ++j) {
    j_counter[j] = ecos_G_jc[j] + 2*(Ainv_jc[j+1]-Ainv_jc[j]);
    if (Ainv_e[j] != 0) {
      j_counter[j] += 2;
    }
  }

  // SOCP part
  for (idxint index1=0; index1<Z; ++index1) {
    oneLine* oneline = readLines[index1];
    idxint id     = oneline->id;
    idxint mother = oneline->mother;
    idxint father = oneline->father;
    // oneline->print();
    idxint target_i = 2 + Z+Z+1+1+1+id;
    if (mother < 0 && father < 0) {
      double b = 4.0 / (2.0*(1-0) + 2.0*(1-0));
      ecos_G_ir[j_counter[id]] = target_i;
      ecos_G_pr[j_counter[id]] = 1.0*sqrt(b);
      j_counter[id]++;
    } 
    else if (father < 0) { // mother >= 0
      double b = 4.0 / (1.0*(1-inbreeding[mother]) + 2.0*(1-0));
      ecos_G_ir[j_counter[mother]] = target_i;
      ecos_G_pr[j_counter[mother]] = -0.5*sqrt(b);
      j_counter[mother]++;
      ecos_G_ir[j_counter[id]] = target_i;
      ecos_G_pr[j_counter[id]] = 1.0*sqrt(b);
      j_counter[id]++;
    } 
    else { // mother >= 0 && father >= 0
      // printf("mother = %d, father = %d\n", mother, father);
      double b = 4.0 / ((1-inbreeding[mother]) + (1-inbreeding[father]));
      if (mother != father) {
	ecos_G_ir[j_counter[father]] = target_i;
	ecos_G_pr[j_counter[father]] = -0.5*sqrt(b);
	j_counter[father]++;
	ecos_G_ir[j_counter[mother]] = target_i;
	ecos_G_pr[j_counter[mother]] = -0.5*sqrt(b);
	j_counter[mother]++;
      } 
      else {
	ecos_G_ir[j_counter[father]] = target_i;
	ecos_G_pr[j_counter[father]] = -2*0.5*sqrt(b);
	j_counter[father]++;
      }
      ecos_G_ir[j_counter[id]] = target_i;
      ecos_G_pr[j_counter[id]] = 1.0*sqrt(b);
      j_counter[id]++;
    }
  }

  // add a column for t
  idxint t_start = ecos_G_jc[Z];
  ecos_G_ir[t_start + 0] = 0;
  ecos_G_pr[t_start + 0] = -1.0;
  ecos_G_ir[t_start + 1] = 1;
  ecos_G_pr[t_start + 1] = 1.0;
  for (idxint index1 = 0; index1 < Z; ++index1) {
    ecos_G_ir[t_start + 2 + index1] = 2 + index1;
    ecos_G_pr[t_start + 2 + index1] = -1.0;
  }
  
  // printf("ecos_G_index = %ld\n", ecos_G_index);
  // printf("Ainv_nnz = %ld", Ainv_nnz);
  #if 0 
  for (idxint j=0; j<Z+1; ++j) {
    idxint G_j_start = ecos_G_jc[j];
    idxint G_j_end   = ecos_G_jc[j+1];
    for (idxint G_j_index = G_j_start; G_j_index < G_j_end; ++G_j_index) {
      if (ecos_G_pr[G_j_index] != 0) {
	idxint i = ecos_G_ir[G_j_index];
	printf("G[%ld,%ld] = %e\n", i,j, ecos_G_pr[G_j_index]);
      }
    }
  }
  #endif
  printf("n = %ld, m = %ld, p = %ld, l = %ld, ",
	 ecos_n, ecos_m, ecos_p, ecos_l);
  printf("ncones = %ld  {", ecos_ncones);
  for (idxint index1=0; index1<ecos_ncones; ++index1) {
    printf(" %ld", ecos_q[index1]);
  }
  printf("}\n");
  printf("Ainv_nnz = %ld, G_nnz = %ld\n", Ainv_nnz, ecos_G_nnz);
  printf("Writing SOCP data into ECOS format...\n");
  fflush(stdout);

  FILE* fp_socp = fopen(ecos_input, "w");
  if (fp_socp == NULL) {
    printf("Cannot open %s\n", ecos_input);
    exit(-1);
  }

  fprintf(fp_socp, "# n\n");
  fprintf(fp_socp, "%ld\n", ecos_n);
  fprintf(fp_socp, "# m\n");
  fprintf(fp_socp, "%ld\n", ecos_m);
  fprintf(fp_socp, "# p\n");
  fprintf(fp_socp, "%ld\n", ecos_p);
  fprintf(fp_socp, "# l\n");
  fprintf(fp_socp, "%ld\n", ecos_l);
  fprintf(fp_socp, "# ncones\n");
  fprintf(fp_socp, "%ld\n", ecos_ncones);
  fprintf(fp_socp, "# q\n");
  fprintf(fp_socp, "%ld\n", ecos_q[0]); // only one cone
  fprintf(fp_socp, "# c %ld\n", ecos_n);
  for (idxint i=0; i<ecos_n; ++i) {
    fprintf(fp_socp, "%.16e\n", ecos_c[i]);
  }
  fprintf(fp_socp, "# h\n");
  for (idxint i=0; i<ecos_m; ++i) {
    fprintf(fp_socp, "%.16e\n", ecos_h[i]);
  }
  fprintf(fp_socp, "# Gjc\n");
  for (idxint i=0; i<ecos_n+1; ++i) {
    fprintf(fp_socp, "%ld\n", ecos_G_jc[i]);
  }
  fprintf(fp_socp, "# Gir\n");
  for (idxint i=0; i<ecos_G_nnz; ++i) {
    fprintf(fp_socp, "%ld\n", ecos_G_ir[i]);
  }
  fprintf(fp_socp, "# Gpr\n");
  for (idxint i=0; i<ecos_G_nnz; ++i) {
    fprintf(fp_socp, "%.16e\n", ecos_G_pr[i]);
  }
  fprintf(fp_socp, "# Ajc\n"); // empty
  fprintf(fp_socp, "# Air\n"); // empty
  fprintf(fp_socp, "# Apr\n"); // empty
  fprintf(fp_socp, "# b\n");  // empty
#if 0
  for (int i=0; i<ecos_p; ++i) {
    fprintf(fp_socp, "%.16e\n", ecos_b[i]);
  }
#endif
  fclose(fp_socp);
  TimeEnd(CONVERT_END);
  printf("Time for conversion = %.3f\n", TimeCal(CONVERT_START, CONVERT_END));


  TimeStart(ECOS_START);
  char system_command[16384];
  sprintf(system_command, "\"%s\" \"%s\" \"%s\" \"%s\" \"%s\" \"%s\"", 
	  ecos_wrapper, ecos_input, ecos_info, ecos_x, ecos_y, ecos_z);
  printf("--->Calling SOCP solver<---\n");
  printf("-->command<-- %s\n", system_command);
  fflush(stdout);
  system(system_command);

  TimeEnd(ECOS_END);

  printf("Time for solver = %.3f\n", TimeCal(ECOS_START, ECOS_END));

  FILE* fp_x = fopen(ecos_x, "r");
  pfloat* result_x = NULL;
  result_x = new double[Z+1];
  for (idxint index1=0; index1<Z+1; ++index1) {
    char readline[MAX_CHARACTER];
    fgets(readline, MAX_CHARACTER, fp_x);
    sscanf(readline, "%lf", &result_x[index1]);
  }

  FILE* fp_info = fopen(ecos_info, "r");
  char readlineInfo[MAX_CHARACTER];
  fgets(readlineInfo, MAX_CHARACTER, fp_info);
  char* token_info = strtok(readlineInfo,"=");
  token_info = strtok(NULL, "=");
  int exitflag = atoi(token_info);
  printf("===> Ecos Status = %d : (", exitflag);
  switch(exitflag) {
  case 0: printf("optimal"); break;
  case 1: printf("primal-infeasible"); break;
  case 2: printf("dual-infeasible"); break;
  case -1: printf("maximum-iteration"); break;
  case -2: printf("numerical-problem[search direction]"); break;
  case -3: printf("numerical-problem[exterior]"); break;
  case -7: printf("unknwon problem"); break;
  default: printf("Cannot read Ecos Status \n"); break;
  }
  printf(") <===\n");

  fclose(fp_info);

  double* contributions = new double[Z];
  double objective_value = 0.0;
  if (exitflag >= -3 && exitflag <= 2) {
    printf("Converting solutions ...\n");
    multi_Ainv_x(Z, Ainv_jc, Ainv_ir, Ainv_pr, result_x, contributions);

    for (idxint index1=0; index1<Z; ++index1) {
      readLines[index1]->contribution = contributions[index1];
#if 0
      if (contributions[index1] > 1.0e-5) {
	printf("contributions[%ld] = %e\n", index1, contributions[index1]);
      }
#endif
    }

    sort(readLines.begin(), readLines.end(), oneLine::compare_back);
    printf("Writing solutions to %s ...\n", filename_output);
    for (idxint index1=0; index1 < ecos_n; ++index1) {
      objective_value += (-ecos_c[index1])*result_x[index1];
    }
    fprintf(fp_output, "# Objective Value = %e\n", objective_value);
    for (idxint index1=0; index1<Z; ++index1) {
      double adjust = N * readLines[index1]->contribution;
      if (adjust >= -1.0e-7 && adjust < 1.0e-7) {
	adjust = 0.0;
      }
      if (adjust != 0) {
	fprintf(fp_output, "%ld, %e\n", 
		readLines[index1]->original_id, adjust);
      }
      else {
	fprintf(fp_output, "%ld, 0.0\n", 
		readLines[index1]->original_id);
      }
    }
    fclose(fp_output);
  } // if (exitflag >= -3 && exitflag <= 2)

  printf("Freeing all memory space ... \n");
  // Free all memories

  for (idxint index1=0; index1<Z; ++index1) {
    delete readLines[index1];
  }
  delete[] result_x;
  delete[] inbreeding;
  delete[] point;
  delete[] D;
  delete[] L;
  for (idxint index1=0; index1<Z; ++index1) {
    vector<j_element*>& A_j = A_add[index1];
    idxint length = A_j.size();
    for (idxint index2=0; index2<length; ++index2) {
      delete A_j[index2];
    }
  }
  delete[] A_add;
  delete[] j_diff_count;

  delete[] Ainv_jc;
  delete[] Ainv_ir;
  delete[] Ainv_pr;
  delete[] ecos_c;
  delete[] minus_c;
  delete[] ecos_h;
  delete[] e_vector;
  delete[] Ainv_e;
  delete[] SOCP_counter;
  delete[] ecos_G_jc;
  delete[] ecos_G_ir;
  delete[] ecos_G_pr;
  delete[] j_counter;
  delete[] contributions;

  printf("--->Finish<----- \n");

  printf("Objective value is %e\n", objective_value);

  TimeEnd(MAIN_END);
  printf("Time for conversion = %.3f\n", TimeCal(CONVERT_START, CONVERT_END));
  printf("Time for solver = %.3f\n", TimeCal(ECOS_START, ECOS_END));
  printf("Total time is %.3lf seconds\n", TimeCal(MAIN_START, MAIN_END));
  
  return 0;
}
