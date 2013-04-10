#include <tclap/CmdLine.h>
#include "NGS++.h"

using namespace std;
using namespace NGS;

/**< This sample requires you to download and make available the TCLAP library */
/**< This requires the BAM file to be indexed and will crash if it is not. No validation is done to insure it is indexed */
int main(int argc, char* argv[])
{
 try
    {
        TCLAP::CmdLine cmd("Signal Extraction from BAM", ' ', "0.2");
        /**< Declare and add arguments */
        TCLAP::ValueArg<std::string> regPath("i","interval","File containing genomic interval regions",true,"null","Filepath",cmd);
        TCLAP::ValueArg<std::string> bamPath("b","bam","Bam file to extract signal from",true,"null","string",cmd);

        /**< Parse */
        cmd.parse( argc, argv );
        ifstream bedFile;

        utility::loadStream(regPath.getValue(),bedFile);
        uParser bedParser(&bedFile,"BED");
        uBasicNGSExperiment intervalExp;
        /**< Load the bedFile */
        intervalExp.loadWithParser(bedParser);
          std::cout << "Loaded "<< intervalExp.count()<<"regions\n";
        /**< Declare BamReader */
        BamTools::BamReader bamParser;
        /**< Open indexed BAM file with BamReader */
        bamParser.Open(bamPath.getValue());
        std::cout << "Opened BAM file\n";
        std::ostream & FILE = (std::cout);

        /**<  Function that will load the data corresponding to the region, load it tu an Experiment, genereate the Signal and write it to standard output. */
        auto writeSignal = [&](const uBasicNGS & oneInterval)
        {
            /**< Reset the BamReader */
            bamParser.Rewind();
            /**< Declare condition for find_if, check for chromosome */
            /**< Note, it would be more efficient to declare once outside scope and std::bind the parameter */
            auto isChrRef=[&](const BamTools::RefData &  item){
                return ( oneInterval.getChr()==item.RefName);
            };
            /**< Check if the scaffold name is found in the Bam data */
            if (std::find_if(bamParser.GetReferenceData().begin(), bamParser.GetReferenceData().end(), isChrRef)!=bamParser.GetReferenceData().end())
            {
                int refID = ( ( std::find_if(bamParser.GetReferenceData().begin(), bamParser.GetReferenceData().end(), isChrRef) - bamParser.GetReferenceData().begin())  );

                /**< Filter the BAMReader to only supply data in the region internal. */
                bamParser.SetRegion(refID,oneInterval.getStart(),refID,oneInterval.getEnd());
                uTagsExperiment bamExperiment;
                bamExperiment.loadWithBamTools(bamParser,0);
                /**< Generate the signal from the internal. */
                vector<float> myRegion = bamExperiment.getRegionSignal(oneInterval.getChr(), oneInterval.getStart(),oneInterval.getEnd(),true);
                std::ostream_iterator<float> output_iterator(FILE);
                /**< Write the signal. */
                cout <<oneInterval.getChr()<<"\t"<<oneInterval.getStart()<<oneInterval.getEnd()<<"\t";
                std::copy(
                          myRegion.begin(),
                          myRegion.end(),
                          output_iterator);
                cout <<"\n";
            }
            else
                {

            }
        };
        /**< Operator the above on every loaded site from the bed file */
        intervalExp.applyOnSites(writeSignal);
    }
    catch (TCLAP::ArgException &e)  // catch any exceptions
    {
        std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;
    }
}








