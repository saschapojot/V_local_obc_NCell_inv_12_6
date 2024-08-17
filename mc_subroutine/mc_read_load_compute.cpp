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
void mc_computation::proposal_uni(const std::shared_ptr<double[]> & xVecCurr,const int&j,std::shared_ptr<double[]>&xVecNext){
    double lm=potFuncPtr->getLm();

    double xVecj_new= generate_uni_open_interval(xVecCurr[j],0,lm,h);

    std::memcpy(xVecNext.get(),xVecCurr.get(),2*N*sizeof (double ));

    xVecNext[j]=xVecj_new;

}



double mc_computation::acceptanceRatio_uni(const std::shared_ptr<double[]> & xVecCurr,const std::shared_ptr<double[]> & xVecNext,const int &j,const double &UCurr, double &UNext){

    double lm=potFuncPtr->getLm();

    UNext=(*potFuncPtr)(xVecNext.get(),j);

    double numerator = -this->beta*UNext;
    double denominator=-this->beta*UCurr;
    double R=std::exp(numerator - denominator);

    double S_xVecCurrNext= S_uni(xVecCurr[j],xVecNext[j],0,lm,h);
    double S_xVecNextCurr= S_uni(xVecNext[j],xVecCurr[j],0,lm,h);

    double ratio_j=S_xVecCurrNext/S_xVecNextCurr;
    if (std::fetestexcept(FE_DIVBYZERO)) {
        std::cout << "Division by zero exception caught." << std::endl;
        std::exit(15);
    }
    if (std::isnan(ratio_j)) {
        std::cout << "The result is NaN." << std::endl;
        std::exit(15);
    }
    R*=ratio_j;

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
            outFile<<",x"+std::to_string(i)+"A";
        }//end even
        else{
            //B, odd position
            outFile<<",x"+std::to_string(i)+"B";
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
        if (ind % 2 == 0) {
            //A, even
            return "xA" + std::to_string(ind - 1);
        } else {
            //B, odd
            return "xB" + std::to_string(ind - 1);
        }

    }//end ind !=0

}
