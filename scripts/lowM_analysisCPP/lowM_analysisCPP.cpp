#include <iostream>
#include <fstream>
#include <unordered_set>
#include "api/BamReader.h"
#include "api/BamWriter.h"
#include "api/BamAlignment.h"
#include <vector>
#include "api/api_global.h"


//iterate over bam file
//for each readName get barcode
//if barcode in unordered set of 1000 barcode, write 


std::unordered_set<std::string> loadBarcodeList(std::string barcodeList){
    std::unordered_set<std::string> barcodeUS;
    std::ifstream file(barcodeList);
    std::string line;
    while(std::getline(file, line)){
        barcodeUS.insert(line);
    }
    return barcodeUS;
}


std::string chr2spName(std::string chr){
    std::string sep_char = "_";
    size_t pos = chr.find_last_of(sep_char);
    std::string sp_name = chr.substr(0, pos);
    return sp_name;
}

std::string trimReadName(std::string readName){
    std::string sep_char = "/";
    size_t pos = readName.find_last_of(sep_char);
    return readName.substr(0, pos);
}

struct Read {
    std::string header;
    std::string sequence;
    std::string plus;
    std::string quality;
};

std::map<std::string, std::string> readName2barcode(const std::string& barcodeFastqFile){
    std::string readPrefix = "E25";
    std::map<std::string, std::string> readName2barcodeMAP;
    std::ifstream barcodeFastqFileIS(barcodeFastqFile);
    std::string line;
    while (std::getline(barcodeFastqFileIS, line)){
        if((line[0] == '@') && (line.find(readPrefix) != std::string::npos)){
            Read read;
            read.header = line;
            std::getline(barcodeFastqFileIS, read.sequence);
            std::getline(barcodeFastqFileIS, read.plus);
            std::getline(barcodeFastqFileIS, read.quality);
            readName2barcodeMAP[trimReadName(read.header)] = read.sequence;
        }
    }
    barcodeFastqFileIS.close();
    return readName2barcodeMAP;
}


void writeAl(BamTools::BamAlignment al, BamTools::SamHeader header, BamTools::RefVector references, std::string sp_name, std::string outDIR){
    std::string outputFilename = outDIR + "/" + sp_name + ".bam";
    BamTools::BamWriter writer;
    if (!writer.Open(outputFilename, header, references)){
        std::cerr << "Could not open output BAM file" << std::endl;
    }
    writer.SaveAlignment(al);
}


int main(int argc, char* argv[]){
    std::string barcodeFastqFile = argv[1];
    std::string barcodes = argv[2];
    std::string inputBam = argv[3];
    std::string outDIR = argv[4];

    std::unordered_set<std::string> barcodeUS = loadBarcodeList(barcodes);
    std::map<std::string, std::string> readName2barcodeMAP = readName2barcode(barcodeFastqFile);

    {
        using namespace BamTools;
        BamReader reader;
        if (!reader.Open(inputBam)){
            std::cerr << "Could not open input BAM files." << std::endl;
            return 0;
        }
        BamAlignment al;
        while ( reader.GetNextAlignmentCore(al) ) {
            std::string readName = al.Name;
            SamHeader header = reader.GetHeader();
            std::string referenceName;
            reader.GetReferenceData()[al.RefID];
            std::string barcode = readName2barcodeMAP[readName];
            if (barcodeUS.find(barcode) != barcodeUS.end()){
                RefVector references = reader.GetReferenceData();
                std::string sp_name = chr2spName(referenceName);
                writeAl(al, header, references, sp_name, outDIR);
            }
        }
    }
    return 1;
}