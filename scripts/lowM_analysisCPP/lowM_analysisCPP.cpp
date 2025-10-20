#include <iostream>
#include <fstream>
#include <unordered_set>
#include "api/BamReader.h"
#include "api/BamWriter.h"
#include "api/BamAlignment.h"
#include <vector>
#include "api/api_global.h"
#include <memory>


//iterate over bam file
//for each readName get barcode
//if barcode in unordered set of 1000 barcode, write 


std::unordered_set<std::string> loadBarcodeList(std::string barcodeList){
    std::cout << "loading barcodes into unordered set..." << std::endl;
    std::unordered_set<std::string> barcodeUS;
    std::ifstream file(barcodeList);
    std::string line;
    while(std::getline(file, line)){
        barcodeUS.insert(line);
    }
    std::cout << "DONE" << std::endl;
    return barcodeUS;
}

std::string chr2spName(std::string chr){
    char sep_char = '_';
    std::string sp_name = "";
    int sep_count = 0;
    for (int idx = 0; idx<chr.length(); idx++){
        if((char)chr[idx] != sep_char){
            sp_name += chr[idx];
        }
        else if((char)chr[idx] == sep_char && sep_count==0){
            sp_name += sep_char;
            sep_count++;
        }
        else if((char)chr[idx] == sep_char && sep_count==1){
            break;
        }
    }
    if (sp_name.find_first_of(sep_char) ==  sp_name.find_last_of(sep_char)){
        return sp_name;
    }
    else{
        std::cerr << "ERROR: too many sep chars in sp name." << std::endl;;
        throw std::runtime_error("ERROR: too many sep chars in sp name.");
    }
}

std::string trimReadName(std::string readName){
    std::string sep_char = "/";
    size_t pos = readName.find_last_of(sep_char);
    if(pos == std::string::npos){
        return readName;
    }
    if(pos > readName.length()){
        return readName;
    }
    return readName.substr(0, pos);
}

struct Read {
    std::string header;
    std::string sequence;
    std::string plus;
    std::string quality;
};

std::map<std::string, std::string> readName2barcode(std::string barcodeFastqFile){
    std::cout << "Making readName2barcode MAP..." << std::endl;
    std::string readPrefix = "E25";
    std::map<std::string, std::string> readName2barcodeMAP;
    std::ifstream barcodeFastqFileIS(barcodeFastqFile);
    std::cout << "barcodeFastqFile was loaded..." << std::endl;
    std::string line;
    while (std::getline(barcodeFastqFileIS, line)){
        if(!line.empty() && line[0] == '@' && line.find(readPrefix) != std::string::npos){
            Read read;
            line.erase(0, 1);
            read.header = line;
            std::getline(barcodeFastqFileIS, read.sequence);
            std::getline(barcodeFastqFileIS, read.plus);
            std::getline(barcodeFastqFileIS, read.quality);
            readName2barcodeMAP[trimReadName(read.header)] = read.sequence.substr(0, 20);
        }
    }
    barcodeFastqFileIS.close();
    std::cout << "DONE" << std::endl;
    return readName2barcodeMAP;
}

std::map<std::string, std::unique_ptr<BamTools::BamWriter>> openHandles(
    const std::vector<std::string>& species,
    const std::string& outDIR,
    const BamTools::SamHeader& header,
    const BamTools::RefVector& references)
{
    std::map<std::string, std::unique_ptr<BamTools::BamWriter>> handles;
    for (const std::string& sp : species) {
        std::string outputFilename = outDIR + "/" + sp + ".bam";
        auto writer = std::make_unique<BamTools::BamWriter>();
        if (!writer->Open(outputFilename, header, references)) {
            std::cerr << "Could not open output BAM file: " << outputFilename << std::endl;
            std::exit(1);
        }
        handles[sp] = std::move(writer);
    }
    return handles;
}


int main(int argc, char* argv[]){
    std::cout << "Processing arguments..." << std::endl;
    std::string barcodeFastqFile = argv[1];
    std::string barcodes = argv[2];
    std::string inputBam = argv[3];
    std::string outDIR = argv[4];
    std::vector<std::string> species = {"Acomys_dimidiatus", "Carassius_gibelio", "Dicrostonyx_torquatus", "Ellobius_talpinus", "Larus_michahellis", "Myotis_brandtii", "Panthera_tigris", "Sander_lucioperca", "Ursus_maritimus"};
    if (!outDIR.empty() && outDIR.back() == '/'){outDIR.pop_back();}
    std::cout << "Running with arguments:" << std::endl;
    std::cout << "barcodeFastqFile:    " << barcodeFastqFile << std::endl;
    std::cout << "barcodes:    " << barcodes << std::endl;
    std::cout << "inputBam:    " << inputBam << std::endl;
    std::cout << "outDIR:    " << outDIR << std::endl;

    std::unordered_set<std::string> barcodeUS = loadBarcodeList(barcodes);
    std::map<std::string, std::string> readName2barcodeMAP = readName2barcode(barcodeFastqFile);

    {
        using namespace BamTools;
        BamReader reader;
        if (!reader.Open(inputBam)){
            std::cerr << "Could not open input BAM files." << std::endl;
            return 1;
        }
        SamHeader header = reader.GetHeader();
        RefVector references = reader.GetReferenceData();
        std::map<std::string, std::unique_ptr<BamTools::BamWriter>> handles = openHandles(species, outDIR, header, references);
        BamAlignment al;
        long int progress = 0;
        long int barcode_match = 0;
        long int barcode_not_found = 0;
        long int sp_not_found = 0;

        std::cout << "Starting itearations over bam file..." << std::endl;
        while ( reader.GetNextAlignment(al) ) {
            progress++;
            if(progress % 1000000 == 0){std::cout << "Current progress: " << progress << "; barcode_match: " << barcode_match << "; barcode_not_found: " << barcode_not_found << std::endl;}
            std::string readName = al.Name;
            
            if (al.RefID < 0) continue;
            if (al.RefID >= static_cast<int32_t>(references.size())) continue;
            auto it = readName2barcodeMAP.find(trimReadName(readName));
            if (it == readName2barcodeMAP.end()) {
                barcode_not_found++;
                continue;
            }
            std::string barcode = it->second;

            if (barcodeUS.find(barcode) != barcodeUS.end()){
                barcode_match++;
                std::string referenceName = references[al.RefID].RefName;

                std::string sp_name = chr2spName(referenceName);
                if (handles.contains(sp_name)){
                    handles[sp_name]->SaveAlignment(al);
                }
                else{
                    sp_not_found++;
                    std::cout << sp_name << " was not found! " << "Created from " << referenceName << std::endl;
                    continue;
                }
            }
        }
        reader.Close();
        for(const auto& pair : handles){
            pair.second->Close();
        } 
    }
    return 0;
}