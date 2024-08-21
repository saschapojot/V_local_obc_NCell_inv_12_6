//
// Created by polya on 8/17/24.
//


#include "mc_read_load_compute.hpp"



///
/// @param x
/// @param leftEnd
/// @param rightEnd
/// @param eps
/// @return return a value within distance eps from x, on the open interval (leftEnd, rightEnd)
double mc_computation::generate_uni_open_interval(const double &x, const double &leftEnd, const double &rightEnd, const double &eps){


    double xMinusEps=x-eps;
    double xPlusEps=x+eps;

    double unif_left_end=xMinusEps<leftEnd?leftEnd:xMinusEps;
    double unif_right_end=xPlusEps>rightEnd?rightEnd:xPlusEps;

//    std::random_device rd;
//    std::ranlux24_base e2(rd());

    double unif_left_end_double_on_the_right=std::nextafter(unif_left_end, std::numeric_limits<double>::infinity());



    std::uniform_real_distribution<> distUnif(unif_left_end_double_on_the_right,unif_right_end); //[unif_left_end_double_on_the_right, unif_right_end)

    double xNext=distUnif(e2);
    return xNext;



}

///
/// @param x proposed value
/// @param y current value
/// @param a left end of interval
/// @param b right end of interval
/// @param epsilon half length
/// @return proposal probability S(x|y)
double mc_computation::S_uni(const double &x, const double &y,const double &a, const double &b, const double &epsilon){

    if (a<y and y<a+epsilon){
        return 1.0/(y-a+epsilon);
    } else if( a+epsilon<=y and y<b+epsilon){
        return 1.0/(2.0*epsilon);
    }else if(b-epsilon<=y and y<b){
        return 1.0/(b-y+epsilon);
    } else{

        std::cerr<<"value out of range."<<std::endl;
        std::exit(10);


    }


}


///
/// @param xVecCurr
/// @param j
/// @param xVecNext
void mc_computation::proposal_uni(const std::shared_ptr<double[]> & xVecCurr,const int&pos,std::shared_ptr<double[]>&xVecNext){
    double lm=potFuncPtr->getLm();

    double xVec_pos_new= generate_uni_open_interval(xVecCurr[pos],0,lm,h);

    std::memcpy(xVecNext.get(),xVecCurr.get(),2*N*sizeof (double ));

    xVecNext[pos]=xVec_pos_new;

}



double mc_computation::acceptanceRatio_uni(const std::shared_ptr<double[]> & xVecCurr,const std::shared_ptr<double[]> & xVecNext,const int &pos,const double &UCurr, double &UNext){

    double lm=potFuncPtr->getLm();

    UNext=(*potFuncPtr)(xVecNext.get(),pos);

    double numerator = -this->beta*UNext;
    double denominator=-this->beta*UCurr;
    double R=std::exp(numerator - denominator);

    double S_xVecCurrNext= S_uni(xVecCurr[pos],xVecNext[pos],0,lm,h);
    double S_xVecNextCurr= S_uni(xVecNext[pos],xVecCurr[pos],0,lm,h);

    double ratio_pos=S_xVecCurrNext/S_xVecNextCurr;
    if (std::fetestexcept(FE_DIVBYZERO)) {
        std::cout << "Division by zero exception caught." << std::endl;
        std::exit(15);
    }
    if (std::isnan(ratio_pos)) {
        std::cout << "The result is NaN." << std::endl;
        std::exit(15);
    }
    R*=ratio_pos;

    return std::min(1.0,R);

}


void mc_computation::saveLastData2Csv(const std::shared_ptr<double[]>& array, const  int& arraySize, const std::string& filename, const int& numbersPerRow){
    //saves last row to csv
    std::ofstream outFile(filename);

    if (!outFile.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return;
    }
    outFile << std::setprecision(std::numeric_limits<double>::digits10 + 1) << std::fixed;
    outFile<<"U";
    for (int i=0;i<2*N;i++){
        if(i%2==0){
            //A, even position
            int j=i/2;
            outFile<<",xA"+std::to_string(j);
        }//end even
        else{
            //B, odd position
            int j=i/2;
            outFile<<",xB"+std::to_string(j);
        }//end odd

    }//end for


    outFile<<"\n";
    for(int i=arraySize-numbersPerRow;i<arraySize;i++){
        outFile<<array[i];
        if ((i + 1) % numbersPerRow == 0) {
            outFile << '\n';
        } else {
            outFile << ',';
        }


    }
    outFile.close();

}


void mc_computation::save_array_to_pickle_one_column(double *ptr, const int& startingInd, std::size_t size,const int & numbersPerRow, const std::string& filename){
    using namespace boost::python;
    try {
        Py_Initialize();  // Initialize the Python interpreter
        if (!Py_IsInitialized()) {
            throw std::runtime_error("Failed to initialize Python interpreter");
        }

        // Debug output
//        std::cout << "Python interpreter initialized successfully." << std::endl;

        // Import the pickle module
        object pickle = import("pickle");
        object pickle_dumps = pickle.attr("dumps");

        // Create a Python list from the C++ array
        list py_list;
        for (std::size_t i = startingInd; i < size; i+=numbersPerRow) {
            py_list.append(ptr[i]);
        }

        // Serialize the list using pickle.dumps
        object serialized_array = pickle_dumps(py_list);

        // Extract the serialized data as a string
        std::string serialized_str = extract<std::string>(serialized_array);

        // Write the serialized data to a file
        std::ofstream file(filename, std::ios::binary);
        if (!file) {
            throw std::runtime_error("Failed to open file for writing");
        }
        file.write(serialized_str.data(), serialized_str.size());
        file.close();

        // Debug output
//        std::cout << "Array serialized and written to file successfully." << std::endl;
    } catch (const error_already_set&) {
        PyErr_Print();
        std::cerr << "Boost.Python error occurred." << std::endl;
    } catch (const std::exception& e) {
        std::cerr << "Exception: " << e.what() << std::endl;
    }

    if (Py_IsInitialized()) {
        Py_Finalize();  // Finalize the Python interpreter
    }




}



std::string mc_computation::generate_varName(const int &ind,const int &numbersPerRow){
    if (ind==0){
        return "U";
    }//end ind=0
    else {
        if (ind % 2 == 1) {
            //A, even, ind is odd
            int j=ind/2;
            return "xA" + std::to_string(j);
        } else {
            //B, odd, ind is even
            int j=ind/2-1;
            return "xB" + std::to_string(j);
        }

    }//end ind !=0

}



void mc_computation::execute_mc_one_sweep(std::shared_ptr<double[]>&xVecCurr,std::shared_ptr<double[]>& xVecNext, const int &fls, const int& swp) {
    double UFull=potFuncPtr->potentialFull(xVecCurr.get());

    //next U
    double UNext;
    double UCurr;
    for (int j = 0; j < 2 * N; j++) {
        //one mc in sweep
        int pos = dist0_2N_minus1(e2);// position to change the value
        //propose next
        this->proposal_uni(xVecCurr, pos, xVecNext);

        //next U


         UCurr = (*potFuncPtr)(xVecCurr.get(), pos);
        // double UCurr=potFuncPtr->potentialFull(xVecCurr.get());
        //accept reject

        double r = this->acceptanceRatio_uni(xVecCurr, xVecNext, pos, UCurr, UNext);

        double u = distUnif01(e2);
        double UCurrCpy=UCurr;
        if (u <= r) {
            std::memcpy(xVecCurr.get(), xVecNext.get(), 2 * N * sizeof(double));
            UCurr = UNext;

        }//end of accept-reject
        UFull+=UCurr-UCurrCpy;
        // U_dist_ptr[swp * 2 * N * varNum + j * varNum+0] = UFull;
        // std::memcpy(U_dist_ptr.get() + swp * 2 * N * varNum + j * varNum + 1, xVecCurr.get(), 2 * N * sizeof(double));

    }//end sweep for
    U_dist_ptr[swp  * varNum +0] = UFull;
    std::memcpy(U_dist_ptr.get() + swp * varNum  + 1, xVecCurr.get(), 2 * N * sizeof(double));


}



void mc_computation::execute_mc(const std::shared_ptr<double[]> &xVec, const int & sweepInit, const int & flushNum) {

    std::shared_ptr<double[]> xVecCurr = std::shared_ptr<double[]>(new double[2 * N], std::default_delete<double[]>());
    std::shared_ptr<double[]> xVecNext = std::shared_ptr<double[]>(new double[2 * N], std::default_delete<double[]>());


    std::memcpy(xVecCurr.get(), xVec.get(), 2 * N * sizeof(double));

    int sweepStart = sweepInit;
    for (int fls = 0; fls < flushNum; fls++) {
        const auto tMCStart{std::chrono::steady_clock::now()};
        for (int swp = 0; swp < sweepToWrite; swp++) {
            execute_mc_one_sweep(xVecCurr, xVecNext, fls, swp);

        }//end sweep for

        int sweepEnd = sweepStart + sweepToWrite - 1;

        std::string fileNameMiddle = "sweepStart" + std::to_string(sweepStart) + "sweepEnd" + std::to_string(sweepEnd);
        std::string out_U_distPickleFileName_pkl = this->U_dist_dataDir + "/" + fileNameMiddle + ".U_dist.pkl";

        std::string out_U_distPickleFileName_csv = this->U_dist_dataDir + "/" + fileNameMiddle + ".U_dist.csv";
        saveLastData2Csv(U_dist_ptr, sweepToWrite  * varNum, out_U_distPickleFileName_csv, varNum);

        for (int startingInd = 0; startingInd < varNum; startingInd++) {
            std::string varName = generate_varName(startingInd, varNum);
            std::string outVarPath = this->U_dist_dataDir + "/" + varName + "/";
            if (!fs::is_directory(outVarPath) || !fs::exists(outVarPath)) {
                fs::create_directories(outVarPath);
            }
            std::string outVarFile = outVarPath + "/" + fileNameMiddle + "." + varName + ".pkl";
            save_array_to_pickle_one_column(U_dist_ptr.get(), startingInd, sweepToWrite  * varNum, varNum,
                                            outVarFile);

        }
        const auto tMCEnd{std::chrono::steady_clock::now()};
        const std::chrono::duration<double> elapsed_secondsAll{tMCEnd - tMCStart};
        std::cout << "sweep " + std::to_string(sweepStart) + " to sweep " + std::to_string(sweepEnd) + ": "
                  << elapsed_secondsAll.count() / 3600.0 << " h" << std::endl;
        sweepStart = sweepEnd + 1;
    }//end flush for loop
    std::cout << "mc executed for " << flushNum << " flushes." << std::endl;


}


void mc_computation::init_and_run(){
    this->execute_mc(xVecInit,sweepLastFile+1,newFlushNum);


}
