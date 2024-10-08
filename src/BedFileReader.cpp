// header file
#include <fstream>  
#include <iostream> 
#include <stdio.h>
#include <bitset>
#include <vector>
#include <cstring>
#include <chrono>

#include <Rcpp.h>
#include "BedFileReader.h"

using namespace std;
using namespace Rcpp;

/// @brief 
/// @param famName 
/// @param bimName 
/// @param bedName  

// BedFileReader::BedFileReader() : famName_temp(""), bimName_temp(""), bedName_temp(""), m_line_counter(0), m_num_of_snps(0), m_size_of_esi(0) {
// }

BedFileReader::BedFileReader(string famName, string bimName, string bedName)
    : famName_temp(famName), bimName_temp(bimName), bedName_temp(bedName) { 
    /* ---------- initialization ----------*/
    this->famName_temp = famName;
    this->bimName_temp = bimName;
    this->bedName_temp = bedName;
    this->m_line_counter = -1; // total number of individual 
    this->m_num_of_snps = -1; // total number of snps

    /* ---------- fam file reading ----------*/
    this->fam.open(famName);    // open file 

    // error messege
    if (!this->fam){
        error.open("reading_file.log");
        error << "CANT_OPEN_FAM_FILE4READ\n";
        error.close();

    }

    if (this->fam.is_open()){        // check if the file is open
        while(!this->fam.eof()){    // while end of the file
            string str;
            getline(this->fam, str);
            this->m_line_counter++; // read fam file to ge number of individuals
        }
        this->fam.close();  
        
    }
    
    cout << "Number of individuals: " << this->m_line_counter << endl;
    this->m_size_of_esi = (this->m_line_counter + 3)/4; //+3 for magic number & mode
    cout << "Number of bytes: " << this->m_size_of_esi << endl;

    /* ---------- bim file reading ----------*/

    this->bim.open(bimName);    

    // error messege
    if (!this->bim){
        error.open("reading_file.log", ios_base::app); // append error messege to log file
        error << "CANT_OPEN_BIM_FILE4READ\n";
        error.close();
    }

   if (this->bim.is_open()){
        while(!this->bim.eof()){
            string str;
            getline(this->bim, str);
            this->m_num_of_snps++;
        }
        this->bim.close();  
        
    }

    cout << "Number of snps: " << this->m_num_of_snps << endl;

    /* ---------- bed file reading ----------*/
    this->bed.open(bedName, ios::binary); // read it as binary

    //error messege
    if (!this->bed){
        error.open("reading_file.log", ios_base::app); 
        error << "CANT_OPEN_BED_FILE4READ\n";
        error.close();

    }

    char tmpbuf[3]; 
    memset(tmpbuf, '\0', sizeof(tmpbuf));
    this->bed.read(tmpbuf, sizeof(tmpbuf)); // three first bytes - permanent in bed file

}

// void BedFileReader::snp_index_func(){
    
//     this->bim.open(this->bimName_temp);

//     // error messege
//     if (!this->bim.is_open()){
//     error.open("reading_file.log", ios_base::app);
//     error << "BIM_FILE_READ_ERROR: cannot open the file\n";
//     error.close();
//     } 

//     size_t curLine = 0;
//     string line;

//     while(getline(this->bim, line)) { 
//         curLine++;
//         std::istringstream iss(line);

//         int column1;
//         string column2;

//         iss >> column1 >> column2; 
//         this->snp_index.insert({column2, curLine});
//     }

//     this->bim.close(); 
// }

void BedFileReader::snp_index_func(){
  this->bim.open(this->bimName_temp);
  // error messege
  if (!this->bim.is_open()){
  error.open("reading_file.log", ios_base::app);
  error << "BIM_FILE_READ_ERROR: cannot open the file\n";
  error.close();
  }
  size_t curLine = 0;
  string line;
  snp_index.clear();

  while(getline(this->bim, line)) {
    curLine++;
    std::istringstream iss(line);
    int column1;
    string column2;
    iss >> column1 >> column2;
    this->snp_index.insert({column2, curLine});
    //cout<<curLine<<column2<<endl;
  }
  //<<"map generation completed, number of snp: "<<curLine<<endl;
  this->bim.close();
  //cout<<"bim file closed"<<endl;
}

// int BedFileReader::findSnpIndex(string snpName){
    
//     this->bim.open(this->bimName_temp); //dictionary (map) : snpid(std::string) as key, snpidx as value

//     // error messege
//     if (!this->bim.is_open()){
//     error.open("reading_file.log", ios_base::app);
//     error << "BIM_FILE_READ_ERROR: cannot open the file\n";
//     error.close();
//     } 

//     size_t curLine = 0;
//     string line;
    
//     while(getline(this->bim, line)) { 
//         curLine++;
        
//         if (line.find(snpName, 0) != string::npos) { //only 2nd col
//             //cout << "found: " << this->snpName << " line: " << curLine << endl;
//             this->snpIndex = curLine;
//         }
//     }
//     this->bim.close();  

//     //cout << this->snpIndex << endl;
//     return this->snpIndex;
// }

int BedFileReader::findSnpIndex(string snpName){
  snpIndex = this->snp_index.find(snpName)->second;
  //cout << snpIndex << endl;
  return snpIndex;
}

vector<int> BedFileReader::readOneSnp(int snpIndex){

    this->snpIndex = snpIndex;
    /* ---------- initialization ----------*/
    int bits_val[MY_CHAR_BIT];
    vector<int> temp_snp_info0(m_line_counter, 0);
    vector<int> temp_snp_info1(m_line_counter, 0); 
    vector<char> encoded_snp_info(m_size_of_esi, 0); 
    vector<char> buff(m_size_of_esi, 0); 
    
    size_t individuals_counter = 0;
    size_t pos = (snpIndex-1) * this->m_size_of_esi + 3;
    size_t pos_cur = this->bed.tellg();

    // error messege
    if (!this->bed.is_open()){
    error.open("reading_file.log", ios_base::app);
    error << "1. BED_FILE_READ_ERROR: cannot open the file\n";
    error.close();
    } 

    if(!pos_cur){
        error.open("reading_file.log", ios_base::app); 
        error << "2. BED_FILE_READ_ERROR: cannot get the position\n";
        error.close();
    }

    if(!this->bed.seekg(pos, ios::beg)){
        error.open("reading_file.log", ios_base::app); 
        error << "3. BED_FILE_READ_ERROR: cannot get seek the position\n";
        error.close();
    }

    this->bed.read((char*) &buff[0],this->m_size_of_esi);

    for (size_t i = 0; i < this->m_size_of_esi; i++) {	
    //===============================================================					
    //=== This part converts Byte "buff[i]" to bits values "bits_val"
    //=== for example byte buff[0] = "w" ->  bits_val = 11101011

        memset(bits_val, 0, sizeof(bits_val));

        int k = MY_CHAR_BIT;  //8
            while (k > 0){
                -- k;
                bits_val[k] = (buff[i]&(1 << k) ? 1 : 0);
            }
    
        //=== here interpret Bit information "bits_val" to snps and count it - decode it
        decode_byte(bits_val, &individuals_counter, &temp_snp_info0[0], &temp_snp_info1[0]);
    
    }

  return temp_snp_info0;
}
// temp_snp_info0 - vector of integers that holds counting for first character 
//					of genotype(in plink format, alt allele) per individual 
// temp_snp_info1 - vector of integers that holds counting for second character 
//					of genotype(in plink format, ref allele) per individual


void BedFileReader::readAllSnp(string fileName){

    chrono::steady_clock::time_point begin = chrono::steady_clock::now();
    ofstream matrix;
    matrix.open(fileName, ios_base::app);
    for (size_t i= 1; i <= this->m_num_of_snps; i++){
        vector<int> a0 = readOneSnp(i); 
        for (size_t j= 0; j < this->m_line_counter; j++){
            matrix << a0[j] << " ";
        }
        matrix << "\n";
        
    }
    matrix.close();
    chrono::steady_clock::time_point end = chrono::steady_clock::now();
    cout << "Time take to read all snps is " << chrono::duration_cast<chrono::milliseconds>(end - begin).count() << "[ms]" << endl;
}

// ADDED
// void BedFileReader::readSomeSnp(string fileName, vector<string> snpNames){

//     chrono::steady_clock::time_point begin = chrono::steady_clock::now();
//  //   ofstream matrix2;
//  //   matrix2.open(fileName, ios_base::app);

//     for (const string& snpName : snpNames) {
// //        cout << snpName << endl;
//         int i = findSnpIndex(snpName);
//         vector<int> a0 = readOneSnp(i);
// /*
//         for (int j = 0; j < this->m_line_counter; j++) {
//             if(j==1){
//                 cout<<a0[j]<<endl;
//             }
//             matrix2 << a0[j] << " ";
//         }
//         matrix2 << "\n";
//  */
//     }
//  //   matrix2.close();
//     chrono::steady_clock::time_point end = chrono::steady_clock::now();
//     cout << "Time take to read some snps is " << chrono::duration_cast<chrono::milliseconds>(end - begin).count() << "[ms]" << endl;
// }

// vector<vector<int>> BedFileReader::readSomeSnp(vector<string> snpNames, vector<int> sampleList = vector<int>()){

//     cout << "snpNames size: " << snpNames.size() << ", sampleList size: " << sampleList.size() << endl;

//     chrono::steady_clock::time_point begin = chrono::steady_clock::now();
//     vector<vector<int>> genomat;

//     for (const string& snpName : snpNames) {
//         //cout << snpName << endl;
//         int i = findSnpIndex(snpName);
//         vector<int> a0 = readOneSnp(i);
//         vector<int> a1 = {};
//         // cout << a0.size() << endl;
        
//         int cnt = 0;
        
//         if (!sampleList.empty()) {
//             for (const int& sampleIdx : sampleList){
//                 if (sampleIdx!=-1){
//                     a1.push_back(a0[sampleIdx]);
//                 }
//                 else{
//                     a1.push_back(9);
//                 }
//                 cnt++;
//             }
//         }
//         else{
//             a1 = a0;
//         }

//         // cout<<a1[1]<<endl;
//         genomat.push_back(a1);
//         // cout << "sample size: " << cnt << endl;
//     }

//     // cout << "genomat row: " << genomat.size() << ", genomat col: " << genomat[1].size() << endl;

//     chrono::steady_clock::time_point end = chrono::steady_clock::now();
//     cout << "Time take to read some snps is " << chrono::duration_cast<chrono::milliseconds>(end - begin).count() << "[ms]" << endl;
    
//     return genomat;
// }

vector<vector<int>> BedFileReader::readSomeSnp(vector<string> snpNames, vector<int> sampleList = vector<int>()){

    cout << "snpNames size: " << snpNames.size() << ", sampleList size: " << sampleList.size() << endl;

    chrono::steady_clock::time_point begin = chrono::steady_clock::now();
    vector<vector<int>> genomat;

    for (const string& snpName : snpNames) {
        int i = findSnpIndex(snpName);
        vector<int> a0 = readOneSnp(i);
        
        if (!sampleList.empty()) {
            vector<int> a1;
            a1.reserve(sampleList.size());
            for (const int& sampleIdx : sampleList){
                if (sampleIdx != -1 && sampleIdx < a0.size()){
                    a1.push_back(a0[sampleIdx]);
                }
                else{
                    a1.push_back(9);
                }
            }
            genomat.push_back(a1);
        }
        else {
            genomat.push_back(a0);
        }
    }

    chrono::steady_clock::time_point end = chrono::steady_clock::now();
    cout << "Time taken to read some snps is " << chrono::duration_cast<chrono::milliseconds>(end - begin).count() << "[ms]" << endl;
    
    return genomat;
}

vector<float> BedFileReader::calculatePRS(vector<string> snpList, vector<float> betaList, vector<int> sampleList = vector<int>()){
    chrono::steady_clock::time_point begin = chrono::steady_clock::now();
    
    cout << "snpList size: " << snpList.size() << ", sampleList size: " << sampleList.size() << endl;

    int snpIndex;
    vector<int> oneSnp; 
    vector<float> PRS;

    // if (!sampleList.empty()) {
    //     PRS.resize(sampleList.size(), 0);
    // } else {
    //     PRS.resize(m_line_counter, 0);
    // }

    for (int i = 0; i < snpList.size(); i++){
        //snpIndex = findSnpIndex(snpList[i]);
        
        snpIndex = this->snp_index.find(snpList[i])->second;
        oneSnp = readOneSnp(snpIndex);
        
        if (!sampleList.empty()) {
            if (i==0) {
                for (const int& sampleIdx : sampleList){
                    if (sampleIdx!=-1){
                        //PRS[j] += oneSnp[sampleIdx] * betaList[i];
                        cout << "oneSNP size: " << oneSnp.size() << ", sampleIdx: " << sampleIdx << "Isin: " << (sampleIdx <= oneSnp.size()) << endl;
                        //PRS.push_back(oneSnp[sampleIdx] * betaList[i]);
                    }
                    else{
                        PRS.push_back(0);
                    }
                }   
            }
            else {
                int j = 0;
                for (const int& sampleIdx : sampleList){
                    if (sampleIdx!=-1){
                        cout << "oneSNP size: " << oneSnp.size() << ", sampleIdx: " << sampleIdx << "Isin: " << (sampleIdx <= oneSnp.size()) << endl;
                        //PRS[j] += oneSnp[sampleIdx] * betaList[i];
                    }
                    // PRS of NA remains 0
                    j++;
                }
            }
        }
        // calculate only sampleList 
        // if (!sampleList.empty()) {
        //     //cout<<"got sample list"<<endl;
        //     int j = 0;
        //     vector<float> a0;
        //     for (const int& sampleIdx : sampleList){
        //         if (sampleIdx!=-1){
        //             PRS[j] += oneSnp[sampleIdx] * betaList[i];
        //         }
        //         // PRS of NA remains 0
        //         j++;
        //     }
        //     //cout<<j<<endl;
        //     //2586 samples for train
        // }
        else {
            if (i==0) {
                for (int j = 0; j < oneSnp.size(); j++){
                    // PRS.push_back(oneSnp[j] * betaList[i]); //sample x 1
                }
            }
            else{
                for (int j = 0; j < oneSnp.size(); j++){
                    // PRS[j] += oneSnp[j] * betaList[i];
                }
            }
        }
    }
    chrono::steady_clock::time_point end = chrono::steady_clock::now();
    cout << "Time take to read some snps is " << chrono::duration_cast<chrono::milliseconds>(end - begin).count() << "[ms]" << endl;

    return PRS;
}

vector<vector<float>> BedFileReader::calculatePRS_mat(vector<string> snpList, vector<vector<float>> betaAll, vector<int> sampleList = vector<int>()){
    chrono::steady_clock::time_point begin = chrono::steady_clock::now();
    
    vector<int> oneSnp; 
    int models = betaAll.size();
    int sampleCount = sampleList.empty() ? m_line_counter : sampleList.size();

    vector<vector<float>> PRS(models, vector<float> (sampleCount, 0.0f));
     // models * samples

    // if (!sampleList.empty()) {
    //     PRS.resize(models, vector<float> (sampleList.size(), 0));
    // } 
    // else {
    //     PRS.resize(m_line_counter, 0);
    // }

    for (int k = 0; k < models; ++k){
        for (int i = 0; i < snpList.size(); ++i){
            //snpIndex = findSnpIndex(snpList[i]);
            
            snpIndex = this->snp_index.find(snpList[i])->second;
            oneSnp = readOneSnp(snpIndex);
            
            // calculate only sampleList 
            if (!sampleList.empty()) {
                //cout<<"got sample list"<<endl;
                int j = 0;
                for (const int& sampleIdx : sampleList){
                    if (sampleIdx!=-1){
                        if (sampleIdx >= oneSnp.size() || j >= PRS[k].size()) {
                            cerr << "Error: Index out of bounds in sampleList processing." << endl;
                            continue;
                        }
                        PRS[k][j] += oneSnp[sampleIdx] * betaAll[k][i];
                    }
                    // PRS of NA remains 0
                    j++;
                }
                //cout<<j<<endl;
                //2586 samples for train
            }
            else{
                for (int j = 0; j < oneSnp.size(); ++j){
                    if (j >= PRS[k].size()) {
                        cerr << "Error: Index out of bounds in full sample processing." << endl;
                        continue;
                    }
                    PRS[k][j] += oneSnp[j] * betaAll[k][i];
                }
            }
        }
    }
    chrono::steady_clock::time_point end = chrono::steady_clock::now();
    cout << "Time take to read some snps is " << chrono::duration_cast<chrono::milliseconds>(end - begin).count() << "[ms]" << endl;

    return PRS;
}

void BedFileReader::decode_byte(int* bits_val, size_t * individuals_counter, int* temp_snp_info0, int* temp_snp_info1){
	//int flag = 0;
	for (int i = 0; i < 4; ++i)
	{
		if(bits_val[i*2] == 0 && bits_val[i*2+1] == 0)
		{
			*individuals_counter += 1;
			if (*individuals_counter > m_line_counter)
				return;
			temp_snp_info0[*individuals_counter - 1] = 2;
			temp_snp_info1[*individuals_counter - 1] = 0;
			//flag = 1 ;//Homozegote 1 for example GG      write 20 ; 00 - will point to [0] letter   + 2 to [0]

		}
		else if(bits_val[i*2] == 1 && bits_val[i*2+1] == 1)
		{
			*individuals_counter += 1;
			if (*individuals_counter > m_line_counter)
				return;
			temp_snp_info0[*individuals_counter - 1] = 0;
			temp_snp_info1[*individuals_counter - 1] = 2;
			//flag = 2 ;//Homozegote 2 for example AA      write 02 ; 11 - will point to [1] letter   + 2 to [1]
		}
		else if(bits_val[i*2] == 1 && bits_val[i*2+1] == 0)
		{
			*individuals_counter += 1;		
			if (*individuals_counter > m_line_counter)
				return;
			temp_snp_info0[*individuals_counter - 1] = 9;
			temp_snp_info1[*individuals_counter - 1] = 9;
			//flag = 3 ; //Missing value                   nothing to add - write 99 ;
		}
		else if(bits_val[i*2] == 0 && bits_val[i*2+1] == 1)
		{
			*individuals_counter += 1;
			if (*individuals_counter > m_line_counter)
				return;
			temp_snp_info0[*individuals_counter - 1] = 1;
			temp_snp_info1[*individuals_counter - 1] = 1;

			//flag = 4 ; //Heterozegote for example AG  or GA     write 11 ; 01 - will point to [0] and [1] letter   +1 +1 to [1]
		}
		// else
		// 	flag = 5 ; //Error
	}
}

void BedFileReader::close_bed(){
  if(this->bed.is_open()){
    this->bed.close();
  }
}
   
// // [[Rcpp::export]]
// RcppExport SEXP BedFileReader__new(SEXP famName_, SEXP bimName_, SEXP bedName_) {
//     // create an external pointer to a Uniform object
//   // convert inputs to appropriate C++ types
//   string famName = Rcpp::as<string>(famName_);
//   string bimName = Rcpp::as<string>(bimName_);
//   string bedName = Rcpp::as<string>(bedName_);

//   // create a pointer to an Uniform object and wrap it as an external pointer
//   Rcpp::XPtr<BedFileReader> reader( new BedFileReader(famName, bimName, bedName), true );

//   // return the external pointer to the R side
//   return reader;
// }

// // [[Rcpp::export]]
// RcppExport SEXP BedFileReader__snp_index_func(SEXP xp){
//     Rcpp::XPtr<BedFileReader> reader(xp);
//     reader->snp_index_func();
// }

// // [[Rcpp::export]]
// RcppExport SEXP BedFileReader__readOneSnp(SEXP xp, SEXP snpIndex_){
//     Rcpp::XPtr<BedFileReader> reader(xp);

//     int snpIndex = Rcpp::as<int>(snpIndex_);

//     vector<int> oneSnp = reader->readOneSnp(snpIndex);
    
//     return Rcpp::wrap(oneSnp);
// }

// // RcppExport SEXP BedFileReader__calculatePRS(SEXP xp, SEXP snpList_, SEXP betaList_, SEXP sampleList_){
// //     Rcpp::XPtr<BedFileReader> reader(xp);

// //     vector<string> snpList = Rcpp::as< vector<string> >(snpList_);
// //     vector<float> betaList = Rcpp::as< std::vector<float> >(betaList_);
// //     vector<int> sampleList = Rcpp::as< std::vector<int> >(sampleList_);

// //     vector<float> PRS = reader->calculatePRS(snpList, betaList, sampleList);

// //     return Rcpp::wrap(PRS);
// // }

// // RcppExport SEXP BedFileReader__calculatePRS_mat(SEXP xp, SEXP snpList_, SEXP betaAll_, SEXP sampleList_){
// //     Rcpp::XPtr<BedFileReader> reader(xp);

// //     vector<string> snpList = Rcpp::as< vector<string> >(snpList_);
// //     vector<vector<float>> betaAll = Rcpp::as< std::vector<std::vector<float>> >(betaAll_);
// //     vector<int> sampleList = Rcpp::as< std::vector<int> >(sampleList_);

// //     vector<vector<float>> PRS = reader->calculatePRS_mat(snpList, betaAll, sampleList);

// //     return Rcpp::wrap(PRS);
// // }

// // [[Rcpp::export]]
// RcppExport SEXP BedFileReader__readAllSnp(SEXP xp, SEXP filename_){
//     Rcpp::XPtr<BedFileReader> reader(xp);

//     string filename = Rcpp::as<string>(filename_);

//     reader->readAllSnp(filename);
// }

// // [[Rcpp::export]]
// RcppExport SEXP BedFileReader__findSnpIndex(SEXP xp, SEXP snpName_){
//     Rcpp::XPtr<BedFileReader> reader(xp);

//     string snpName = Rcpp::as<string>(snpName_);

//     int snpIndex = reader->findSnpIndex(snpName);
//     return Rcpp::wrap(snpIndex);
// }

RCPP_MODULE(BedFileReader_module) {
    class_<BedFileReader>("BedFileReader")
        .constructor<std::string, std::string, std::string>()
        .method("snp_index_func", &BedFileReader::snp_index_func, "Mapping a SNP name to its index in bim")
        .method("findSnpIndex", &BedFileReader::findSnpIndex, "Find index of a single SNP")
        .method("readOneSnp", &BedFileReader::readOneSnp, "Read a single SNP based on index");
}



