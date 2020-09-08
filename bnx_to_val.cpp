#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
using namespace std;

vector<float> reff;

//converts the simulated Rmaps from bnx format to the valouev format (we use this format for processing)

int main(int argc, char *argv[])
{
    // read_ref();

    // check for correct number of arguments
    if (argc < 3)
    {
        // bnx is output of bionano genomics - give path and filename
        // valouev format - has .val extension
        // first line is identifier string
        // second line has name of enzyme
        // after that is rmap--in fragment representation, as opposed to position representation (like in bnx)
        // rmm format - just give output file name
        // enzyme restriction - use BspQI
        cout << "Usage: bnxcnv <bnx_filename> <rmm_filename> <enzyme_restriction_pattern>" << endl;
        return 1;
    }

    int min_size = 0;
    string line;

    ifstream bnxfile(argv[1]); // input stream class to operate on bnx file

    string fname(argv[2]); // output file name

    string enz(argv[3]); // restriction enzyme name

    ofstream valfile((fname + ".val").c_str()); // .val output (easier to parse than bnx for this experiment)

    ofstream cometfile((fname + ".com").c_str()); // .com format file (dont worry about this, simpler valouev format)

    ofstream rmmfile((fname).c_str()); // rmm format output

    vector<vector<int>> rmap_values;

    vector<vector<int>> pos_values;

    vector<vector<float>> intensity_values;
    vector<vector<float>> snr_values;

    int count = 0; // holds total number of Rmaps
    long tot_frag_size = 0;
    long tot_frags = 0;

    long frag_size = 0;
    int number_frags = 0;

    if (bnxfile.is_open()) // check if bnxfile is able to be opened, throws error otherwise
    {
        while (getline(bnxfile, line)) // get current line of bnxfile until failbit is set
        {
            istringstream iss(line); // creates stream from current line
            string word;
            iss >> word; // puts first "word" in iss into word
            if (word.compare("1") != 0)
            {   // if the current word at the start of the stream is 0,
                // indicating that line is the molecule header,
                // run loop again and get next line
                continue;
            }
            else // if it is indeed the data, proceed...
            {
                vector<int> rmap;

                vector<int> pos;
                int old = 0;
                while (iss >> word) // while line stream has another white-space separated word
                {
                    int temp;       // holds current 
                    istringstream(word) >> temp;
                    pos.push_back(temp);
                    if (old == 0)
                    {
                        old = temp;
                        continue;
                    }
                    else
                    {
                        rmap.push_back(temp - old);
                    }

                    tot_frag_size += temp - old;
                    tot_frags++;
                    old = temp;
                }

                if (rmap.size() < min_size)
                    continue;

                string snr_line;
                getline(bnxfile, snr_line);

                string temst;

                istringstream iss_snr(snr_line);

                float temfl;

                iss_snr >> temst;

                vector<float> snr;

                while (iss_snr >> temfl)
                {
                    snr.push_back(temfl);
                }

                string int_line;
                getline(bnxfile, int_line);

                istringstream iss_int(int_line);

                iss_int >> temst;

                vector<float> intensity;

                while (iss_int >> temfl)
                {
                    intensity.push_back(temfl);
                }

                rmap_values.push_back(rmap);
                pos_values.push_back(pos);
                intensity_values.push_back(intensity);
                snr_values.push_back(snr);
                count++;
            }
        }
        bnxfile.close(); // when no more lines are found, close the ifstream
    }

    else
        cout << "Unable to open file";

    rmmfile << "r\t" << rmap_values.size() << "\t 1 \t" << enz.size() << "\t" << enz << endl;

    for (int i = 0; i < pos_values.size(); i++)
    {

        rmmfile << "R \t" << pos_values[i].back() << "\t" << pos_values[i].size() << "\t";

        for (int j = 0; j < pos_values[i].size(); j++)
        {
            rmmfile << pos_values[i][j] << "\t";
        }

        rmmfile << endl;

        rmmfile << "I \t" << intensity_values[i].size() << "\t";

        for (int j = 0; j < intensity_values[i].size(); j++)
        {
            rmmfile << intensity_values[i][j] << "\t";
        }
        rmmfile << endl;

        rmmfile << "N \t" << snr_values[i].size() << "\t";

        for (int j = 0; j < snr_values[i].size(); j++)
        {
            rmmfile << snr_values[i][j] << "\t";
        }
        rmmfile << endl;

        cometfile << "Rmap_" << i;
        valfile << "Rmap_" << i << "\n";
        valfile << enz << "\t" << enz << "\t";

        for (int j = 0; j < rmap_values[i].size(); j++)
        {
            valfile << " " << (float)rmap_values[i][j] / 1000;
            cometfile << " " << (float)(float)rmap_values[i][j] / 1000;
        }

        valfile << "\n \n";
        cometfile << "\n";
    }

    cout << "Total number of Rmaps: " << count << endl;
    //  cout<<"Average fragment size: "<<(float)frag_size/number_frags;;
    return 0;
}
