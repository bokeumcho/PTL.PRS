// header file
#include <fstream>  
#include <iostream> 
#include <stdio.h>
#include <bitset>
#include <vector>
#include <cstring>
#include <chrono>

#include <Rcpp.h>
#include <RcppParallel.h>
#include "BedFileReader.h"

using namespace std;
using namespace Rcpp;
using namespace RcppParallel;

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
    
    Rcout << "Number of individuals: " << this->m_line_counter << endl;
    this->m_size_of_esi = (this->m_line_counter + 3)/4; //+3 for magic number & mode
    Rcout << "Number of bytes: " << this->m_size_of_esi << endl;

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

    Rcout << "Number of snps: " << this->m_num_of_snps << endl;

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

void BedFileReader::snp_index_func() {
    // Open the file
    bim.open(bimName_temp);
    if (!bim.is_open()) {
        std::ofstream error("reading_file.log", std::ios_base::app);
        if (error.is_open()) {
            error << "BIM_FILE_READ_ERROR: Cannot open file: " << bimName_temp << "\n";
        }
        return;
    }

    size_t curLine = 0;
    std::string line;
    snp_index.clear();

    // Process each line of the file
    while (std::getline(bim, line)) {
        ++curLine;
        std::istringstream iss(line);
        int col1;
        std::string col2;
        if (!(iss >> col1 >> col2)) {
            // Log parse error but continue reading the file
            std::ofstream error("reading_file.log", std::ios_base::app);
            if (error.is_open()) {
                error << "BIM_FILE_PARSE_ERROR: Could not parse line " << curLine << "\n";
            }
            continue;
        }
        // Insert or update the SNP index mapping
        snp_index[col2] = curLine;
    }

    bim.close();
}

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

    // chrono::steady_clock::time_point begin = chrono::steady_clock::now();
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
    // chrono::steady_clock::time_point end = chrono::steady_clock::now();
    // cout << "Time take to read all snps is " << chrono::duration_cast<chrono::milliseconds>(end - begin).count() << "[ms]" << endl;
}

// ADDED
// // The parallel worker for processing SNPs
// struct SNPWorker : public Worker {
//     // Input data
//     const vector<string>& snpNames;
//     const vector<int>& sampleList;
//     // Output data (preallocated)
//     vector<vector<double>>& genomat;
//     // Pointer to the parent object to call its methods
//     BedFileReader* reader;

//     // Constructor
//     SNPWorker(const vector<string>& snpNames,
//                 const vector<int>& sampleList,
//                 vector<vector<double>>& genomat,
//                 BedFileReader* reader)
//         : snpNames(snpNames), sampleList(sampleList), genomat(genomat), reader(reader) { }

//     // Overloaded operator() to process a range of SNP names
//     void operator()(std::size_t begin, std::size_t end) {
//         for (std::size_t i = begin; i < end; i++) {
//             const string &snpName = snpNames[i];
//             int idx = reader->findSnpIndex(snpName);
//             vector<int> a0 = reader->readOneSnp(idx);
//             vector<double> a1;

//             if (!sampleList.empty()) {
//                 a1.reserve(sampleList.size());
//                 for (const int &sampleIdx : sampleList) {
//                     if (sampleIdx != -1 && static_cast<size_t>(sampleIdx) < a0.size()) {
//                         a1.push_back(static_cast<double>(a0[sampleIdx]));
//                     } else {
//                         a1.push_back(9.0);
//                     }
//                 }
//             } else {
//                 a1.resize(a0.size());
//                 for (size_t j = 0; j < a0.size(); j++) {
//                     a1[j] = static_cast<double>(a0[j]);
//                 }
//             }
//             genomat[i] = a1;
//         }
//     }
// };

// // Parallelized function using RcppParallel
// vector<vector<double>> BedFileReader::readSomeSnpParallel(const vector<string>& snpNames, const vector<int>& sampleList = vector<int>()) {
//     // Rcout << "currently using" << defaultNumThreads() << "threads in parallel" << endl;
    
//     // Preallocate output vector maintaining order (one per SNP)
//     vector<vector<double>> genomat(snpNames.size());

//     // Create our worker instance
//     SNPWorker worker(snpNames, sampleList, genomat, this);

//     // Run in parallel using parallelFor from RcppParallel
//     parallelFor(0, snpNames.size(), worker);

//     // Optional: Perform mean imputation per row as in your original code.
//     for (auto &row : genomat) {
//         double sum = 0.0;
//         int count = 0;
//         for (double val : row) {
//             if (val != 9.0) {
//                 sum += val;
//                 count++;
//             }
//         }
//         if (count > 0) {
//             double mean_val = sum / count;
//             for (double &val : row) {
//                 if (val == 9.0) {
//                     val = mean_val;
//                 }
//             }
//         }
//     }

//     return genomat;
// }

vector<vector<double>> BedFileReader::readSomeSnp(vector<string> snpNames, vector<int> sampleList = vector<int>()){
    // cout << "snpNames size: " << snpNames.size() << ", sampleList size: " << sampleList.size() << endl;

    // auto begin = chrono::steady_clock::now();
    vector<vector<double>> genomat;

    for (const string& snpName : snpNames) {
        int idx = findSnpIndex(snpName);
        vector<int> a0 = readOneSnp(idx);
        vector<double> a1;

        if (!sampleList.empty()) {
            a1.reserve(sampleList.size());
            for (const int& sampleIdx : sampleList) {
                if (sampleIdx != -1 && static_cast<size_t>(sampleIdx) < a0.size()){
                    a1.push_back(static_cast<double>(a0[sampleIdx]));
                }
                else {
                    a1.push_back(9.0);
                }
            }
        }
        else {
            a1.resize(a0.size());
            for (size_t j = 0; j < a0.size(); j++){
                a1[j] = static_cast<double>(a0[j]);
            }
        }
        genomat.push_back(a1);
    }

    // Perform mean imputation for each SNP (row):
    // Replace each 9.0 with the mean of the non-9 values in that row.
    for (auto& row : genomat) {
        double sum = 0.0;
        int count = 0;
        for (double val : row) {
            if (val != 9.0) {
                sum += val;
                count++;
            }
        }
        if (count > 0) {
            double mean_val = sum / count;
            for (double& val : row) {
                if (val == 9.0) {
                    val = mean_val;
                }
            }
        }
    }

    // auto end = chrono::steady_clock::now();
    // cout << "Time taken to read some snps is " 
    //      << chrono::duration_cast<chrono::milliseconds>(end - begin).count() 
    //      << "[ms]" << endl;
    
    return genomat;
}

vector<float> BedFileReader::calculatePRS(vector<string> snpList, vector<float> betaList, vector<int> sampleList = vector<int>()){
    // chrono::steady_clock::time_point begin = chrono::steady_clock::now();
    
    // cout << "snpList size: " << snpList.size() << ", sampleList size: " << sampleList.size() << endl;

    int snpIndex;
    vector<int> oneSnp; 
    vector<float> PRS;

    // if (!sampleList.empty()) {
    //     PRS.resize(sampleList.size(), 0);
    // } else {
    //     PRS.resize(m_line_counter, 0);
    // }

    for (size_t i = 0; i < snpList.size(); i++){
        //snpIndex = findSnpIndex(snpList[i]);
        
        snpIndex = this->snp_index.find(snpList[i])->second;
        oneSnp = readOneSnp(snpIndex);
        
        if (!sampleList.empty()) {
            if (i==0) {
                for (const int& sampleIdx : sampleList){
                    if (sampleIdx!=-1){
                        //PRS[j] += oneSnp[sampleIdx] * betaList[i];
                        // cout << "oneSNP size: " << oneSnp.size() << ", sampleIdx: " << sampleIdx << "Isin: " << (static_cast<size_t>(sampleIdx) <= oneSnp.size()) << endl;
                        // PRS.push_back(oneSnp[sampleIdx] * betaList[i]);
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
                        // cout << "oneSNP size: " << oneSnp.size() << ", sampleIdx: " << sampleIdx << "Isin: " << (static_cast<size_t>(sampleIdx) <= oneSnp.size()) << endl;
                        //PRS[j] += oneSnp[sampleIdx] * betaList[i];
                    }
                    // PRS of NA remains 0
                    j++;
                }
            }
        }
        else {
            if (i==0) {
                for (size_t j = 0; j < oneSnp.size(); j++){
                    // PRS.push_back(oneSnp[j] * betaList[i]); //sample x 1
                }
            }
            else{
                for (size_t j = 0; j < oneSnp.size(); j++){
                    // PRS[j] += oneSnp[j] * betaList[i];
                }
            }
        }
    }
    // chrono::steady_clock::time_point end = chrono::steady_clock::now();
    // cout << "Time take to read some snps is " << chrono::duration_cast<chrono::milliseconds>(end - begin).count() << "[ms]" << endl;

    return PRS;
}

vector<vector<float>> BedFileReader::calculatePRS_mat(vector<string> snpList, vector<vector<float>> betaAll, vector<int> sampleList = vector<int>()){
    // chrono::steady_clock::time_point begin = chrono::steady_clock::now();
    
    vector<int> oneSnp; 
    int models = betaAll.size();
    int sampleCount = sampleList.empty() ? m_line_counter : sampleList.size();

    vector<vector<float>> PRS(models, vector<float> (sampleCount, 0.0f));
     // models * samples

    for (int k = 0; k < models; ++k){
        for (size_t i = 0; i < snpList.size(); ++i){
            //snpIndex = findSnpIndex(snpList[i]);
            
            snpIndex = this->snp_index.find(snpList[i])->second;
            oneSnp = readOneSnp(snpIndex);
            
            // calculate only sampleList 
            if (!sampleList.empty()) {
                //cout<<"got sample list"<<endl;
                size_t j = 0;
                for (const int& sampleIdx : sampleList){
                    if (sampleIdx!=-1){
                        if (static_cast<size_t>(sampleIdx) >= oneSnp.size() || j >= PRS[k].size()) {
                            Rcout << "Error: Index out of bounds in sampleList processing." << endl;
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
                for (size_t j = 0; j < oneSnp.size(); ++j){
                    if (j >= PRS[k].size()) {
                        Rcout << "Error: Index out of bounds in full sample processing." << endl;
                        continue;
                    }
                    PRS[k][j] += oneSnp[j] * betaAll[k][i];
                }
            }
        }
    }
    // chrono::steady_clock::time_point end = chrono::steady_clock::now();
    // cout << "Time taken to calculate PRS matrix is " << chrono::duration_cast<chrono::milliseconds>(end - begin).count() << "[ms]" << endl;

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
   

RCPP_MODULE(BedFileReader_module) {
    class_<BedFileReader>("BedFileReader")
        .constructor<std::string, std::string, std::string>()
        .method("snp_index_func", &BedFileReader::snp_index_func, "Mapping a SNP name to its index in bim")
        .method("findSnpIndex", &BedFileReader::findSnpIndex, "Find index of a single SNP")
        .method("readOneSnp", &BedFileReader::readOneSnp, "Read a single SNP based on index")
        .method("readSomeSnp", &BedFileReader::readSomeSnp, "Read multiple SNPs based on index")
        // .method("readSomeSnpParallel", &BedFileReader::readSomeSnpParallel, "Read multiple SNPs based on index parallelly")
        .method("calculatePRS_mat", &BedFileReader::calculatePRS_mat, "Calculate PRS of multiple SNPs & betas");
    }



