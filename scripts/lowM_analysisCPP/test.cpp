#include "api/BamMultiReader.h"
#include "api/BamWriter.h"
#include <vector>
#include <string>

int main(){
    using namespace std;
    using namespace BamTools;
    // at some point, start our merge operation
    vector<string> inputFilenames;
    string outputFilename;
    // provide some input & output filenames
    // attempt to open our BamMultiReader
    BamMultiReader reader;
    if (!reader.Open(inputFilenames)) {
        cerr << "Could not open input BAM files." << endl;
    }
    // retrieve 'metadata' from BAM files, these are required by BamWriter
    const SamHeader header = reader.GetHeader();
    const RefVector references = reader.GetReferenceData();
    // attempt to open our BamWriter
    BamWriter writer;
    if ( !writer.Open(outputFilename, header, references) ) {
        cerr << "Could not open output BAM file" << endl;
    }
    // iterate through all alignments, only keeping ones with high map quality
    BamAlignment al;
    while ( reader.GetNextAlignmentCore(al) ) {
        if ( al.MapQuality >= 90 )
        writer.SaveAlignment(al);
    }
    // close the reader & writer
    reader.Close();
    writer.Close();
    // merge is now complete, continue whatever we were doing
    return 1;
}