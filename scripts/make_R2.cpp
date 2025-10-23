#include <iostream>
#include <fstream>
#include <string>
#include <unordered_set>
#include <vector>
#include <map>

using namespace std;

struct Read {
    string header;
    string sequence;
    string plus;
    string quality;
};

string extractName(const string& header) {
    size_t slashPos = header.find('/');
    if (slashPos != string::npos) {
        return header.substr(1, slashPos - 1);
    }
    return "";
}


map<string, Read> readFastq2Map(const string& filename) {
    cout << "reading Fastq file to MAP: " << filename << endl;
    map <string, Read> reads;
    ifstream file(filename);
    string line;

    while (getline(file, line)) {
        if (line[0] == '@') {
            Read read;
            read.header = line;
            getline(file, read.sequence);
            getline(file, read.plus);
            getline(file, read.quality);
            reads[extractName(read.header)] = read;
        }
    }
    return reads;
}


vector<Read> readFastq(const string& filename) {
    vector<Read> reads;
    ifstream file(filename);
    string line;
    while (getline(file, line)) {
        if (line[0] == '@') {
            Read read;
            read.header = line;
            getline(file, read.sequence);
            getline(file, read.plus);
            getline(file, read.quality);
            reads.push_back(read);
        }
    }
    return reads;
}



int main(int argc, char* argv[]) {
    if (argc < 3) {
        cerr << "Usage: program <merged_R1_file> <list_of_R2_files>" << endl;
        return 1;
    }
    string mergedR1File = argv[1];
    map<string,Read> mergedR1Reads = readFastq2Map(mergedR1File);


    cout << "Running processing, iterating over R2 fastq files..." << endl;
    for (int i = 2; i < argc; ++i) {
        string r2File = argv[i];
        cout << "working with fastq file: " << r2File << endl;
        vector<Read> r2Reads = readFastq(r2File);

        string r1OutputFile = r2File + ".R1.fastq";
        ofstream out(r1OutputFile);

        for (const auto& R2read : r2Reads) {
            string readName = extractName(R2read.header);
            Read R1read = mergedR1Reads[readName];
            out << R1read.header << endl;
            out << R1read.sequence << endl;
            out << R1read.plus << endl;
            out << R1read.quality << endl;

        }
        cout << "Created file: " << r1OutputFile << endl;
    }
    return 0;
}

