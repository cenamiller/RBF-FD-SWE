#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>

using namespace std;

void throwout_ss(istringstream&, int);
void throwout_file(ifstream&, int);

// Read in results as .csv file
int main() {
    ifstream inFile;
    ofstream outFile;
    string x;
    string read_val;
    int nskip;
    int nattempts;
    bool openCL;
    string fnames[4] = {"output/output_cfdl.txt","output/output_cfdl_sfdl.txt","output/output_sfdl.txt","output/output_default.txt"};
    outFile.open("out.csv");
    
    for (int j = 0; j < 4; j++) {
        
        inFile.open(fnames[j]);
        if (!inFile) {
            cout << "Unable to open " << fnames[j] << endl;
            continue;
        }
        
        outFile << fnames[j] << endl;
        
        while (!inFile.eof()) {
            
            // Read in script name
            throwout_file(inFile, 1);
            if (!getline(inFile, x)) break;
            istringstream script_name(x);
            throwout_ss(script_name, 3);
            script_name >> read_val;
            openCL = (read_val == "runCL.sh");
            
            if (!openCL) {
               
                // Read in the number of threads
                if (!getline(inFile, x)) break;
                istringstream nunits_line(x);
                throwout_ss(nunits_line, 2);
                nunits_line >> read_val;
                outFile << read_val << ',';
            }
                
            // Read in the number of attempts
            throwout_file(inFile, 8);
            if (!getline(inFile, x)) break;
            istringstream attempts_line(x);
            throwout_ss(attempts_line, 3);
            attempts_line >> nattempts;
            
            if (openCL) {
                
                // Read in the number of subdevices
                throwout_file(inFile, 26);
                if (!getline(inFile, x)) break;
                istringstream nsub_line(x);
                throwout_ss(nsub_line, 3);
                nsub_line >> read_val;
                outFile << read_val << ',';
            }
            
            // Read in the eval_rhs time
            if (openCL) nskip = 9;
            else nskip = 17;
            throwout_file(inFile, nskip + nattempts);
            if (!getline(inFile, x)) break;
            istringstream eval_rhs_line(x);
            throwout_ss(eval_rhs_line, 4);
            eval_rhs_line >> read_val;
            outFile << read_val << endl;
            
            throwout_file(inFile, 6);
        }
        inFile.close();
        outFile << endl;
    }
    outFile.close();
    return 0;
}

void throwout_ss(istringstream& xss, int num_throw) {
    string throwout;
    for (int i = 0; i < num_throw; i++) xss >> throwout;
    return;
}

void throwout_file(ifstream& inFile, int num_throw) {
    string throwout;
    for (int i = 0; i < num_throw; i++) getline(inFile, throwout);
    return;
}
