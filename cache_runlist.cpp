#include "TChain.h"
#include "TTreeReader.h"
#include "TNtuple.h"

#include "VASkyMap.h"
#include "VSOptions.hpp"
#include "VSDataConverter.hpp"
#include "VATime.h"
#include "VASimpleCuts.h"

#include <boost/algorithm/string.hpp>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <string>


double ZABINS[] = {0.00,14.1,25.46,32.59,37.57,42.56,47.55};
double EBINS[] = {316,630,1259,2512,5012};
int TBINS[] = {3,4};
double AZBINS[] = {0,45,90,135,180,225,270,315,360};
double OBINS[] = {0.0, 1.0, 1.4, 1.7, 2.0};
typedef struct Cut_Params{
    double eventRa;
    double eventDec;
    double zenith;
    double energy;
    double ntels;
    double azimuth;
    double offset;
    double msl;
    double height;
    int zbin;
    int ebin;
    int abin;
    int obin;
} Cut_Params_t;

int cuts_za = 0;
int cuts_e = 0;
int cuts_tel = 0;
int cuts_az = 0;
int cuts_off = 0;

typedef struct SourceCut
{
  VACoordinatePair coords;
  double cut_rad;
  SourceCut(std::vector<std::string> line_in) :
  coords(std::stod(line_in[1]),std::stod(line_in[2]),VACoordinates::J2000,VACoordinates::Deg),
  cut_rad(std::stod(line_in[3]))
  {}

  double DistFromSource(VACoordinatePair evt_coords){
    return coords.angularSeparation_Deg(evt_coords);
  }

  bool InsideExclRadius(VACoordinatePair evt_coords) {
    double dist = coords.angularSeparation_Deg(evt_coords);
    return (dist < cut_rad);
  }
} SourceCut_t;

int bin_event(Cut_Params_t *params){
    bool fail = false;
    /*
    //ZA cut
    int z_bin = -42;
    for(int z = 0; z < 6; z++){
        if(params->zenith < ZABINS[z+1] && params->zenith > ZABINS[z]){
            z_bin = z;
        }
    }
    if(z_bin == -42){
        fail = true;
        cuts_za++;
    }


    //Energy Cut
    int e_bin = -42;
    for(int e = 0; e < 4; e++){
        if(params->energy < EBINS[e+1] && params->energy > EBINS[e]){
            e_bin = e;
        }
    }
    if(e_bin == -42){
        fail = true;
        cuts_e++;
    }

    //Tel Cut
    if(params->ntels == 2){
        fail = true;
        cuts_tel++;
    }

    //AZ Cut
    int a_bin = -42;
    for(int a = 0; a < 8; a++){
        if(params->azimuth < AZBINS[a+1] && params->azimuth > AZBINS[a]){
            a_bin = a;
        }
    }
    if(a_bin == -42){
        fail = true;
        cuts_az++;
    }

    //Offset Cut
    int o_bin = -42;
    for(int o = 0; o < 8; o++){
        if(params->offset < OBINS[o+1] && params->offset > OBINS[o]){
            o_bin = o;
        }
    }
    if(o_bin == -42){
        fail = true;
        cuts_off++;
    }
    */

    //Tel Cut
    if(params->ntels < 3){
      fail = true;
      cuts_tel++;
    }

    //Offset Cut
    if(params->offset > 1.75){
      fail = true;
      cuts_off++;
    }

    return fail;
}





void cacheData_vegas(std::string pathbase, bool binned_output){
    std::stringstream excl_file;
    excl_file << pathbase << "_src.txt";
    //Check for source exclusion file
    if(access(excl_file.str().c_str(), F_OK)){
        std::cerr << "No source exclusion file found at "
        << excl_file.str()
        << ". This may cause problems."
        << std::endl;
    }

    std::vector<SourceCut_t> sourcecuts;
    std::string line;
    std::vector<std::string> line_fields;
    std::ifstream srclist(excl_file.str().c_str());
    while(std::getline(srclist, line)){
        boost::split(line_fields, line, boost::is_any_of(","));
        sourcecuts.emplace_back(line_fields);
        line_fields.clear();
    }



    //Setup flist
    std::stringstream infile;
    infile << pathbase << ".list";
    std::ifstream flist(infile.str().c_str());

    //Cache events
    TFile *f;
    std::stringstream outpath;
    outpath << "Data_" << pathbase << ".root";
    TFile g(outpath.str().c_str(), "RECREATE");

    TNtuple *ntuple = new TNtuple("ntuple", "Passed Events",
      "MSW:MSL:Zenith:Energy_GeV:NTels:Azimuth:Offset:RA:DEC:BDT:Src_Excl:Src_Dist:ShowerHeightMax_KM");
    while(std::getline(flist, line)){
        f = TFile::Open(line.c_str());
        if(!f){
            std::cout << "Cannot open " << line << std::endl;
            continue;
        }
        else{
            std::cout << "Reading " << line << std::endl;
        }

        TTreeReader reader("SelectedEvents/CombinedEventsTree", f);
        TTreeReaderValue<VAShowerData> shower(reader, "S");
        TTreeReaderValue<double> bdt(reader, "BDTScore");

        int nRead = 0;
        int nPass = 0;

        while(reader.Next()){
            nRead++;
            Cut_Params_t *params = new Cut_Params_t();
            params->eventRa = shower->fDirectionRA_J2000_Rad * TMath::RadToDeg();
            params->eventDec = shower->fDirectionDec_J2000_Rad * TMath::RadToDeg();
            params->zenith = 90.0 - (shower->fDirectionElevation_Rad * TMath::RadToDeg());
            params->energy = shower->fEnergy_GeV;
            auto tels_used = shower->fTelUsedInReconstruction;
            params->ntels = count(tels_used.begin(), tels_used.end(), 1);
            params->azimuth = shower->fDirectionAzimuth_Rad * TMath::RadToDeg();
            VACoordinatePair shower_coords = VACoordinatePair(
                shower->fDirectionRA_J2000_Rad,
                shower->fDirectionDec_J2000_Rad,
                VACoordinates::J2000,
                VACoordinates::Rad
            );
            VACoordinatePair tracking_coords = VACoordinatePair(
                shower->fArrayTrackingRA_J2000_Rad,
                shower->fArrayTrackingDec_J2000_Rad,
                VACoordinates::J2000,
                VACoordinates::Rad
            );
            params->offset = tracking_coords.angularSeparation_Deg(shower_coords);
            params->msl = shower->fMSL;
            params->height = shower->fShowerMaxHeight_KM;

            if(bin_event(params)) continue;
            nPass++;

            //Check Source Distance
            double in_excl = false;
            double src_dist = 200;
            VACoordinatePair eventCoord = VACoordinatePair(
              params->eventRa,
              params->eventDec,
              VACoordinates::J2000,
              VACoordinates::Deg);
            for(auto &it : sourcecuts){
              if(it.InsideExclRadius(eventCoord)) in_excl = true;
              double dist = it.DistFromSource(eventCoord);
              if(dist < src_dist) src_dist = dist;
            }

            ntuple->Fill(shower->fMSW,
                         shower->fMSL,
                         params->zenith,
                         params->energy,
                         params->ntels,
                         params->azimuth,
                         params->offset,
                         params->eventRa,
                         params->eventDec,
                         *bdt,
                         in_excl,
                         src_dist,
                         shower->fShowerMaxHeight_KM);

            delete params;
        }
        f->Close();
    }

    std::cout << "Zenith  cuts " << cuts_za << std::endl;
    std::cout << "Energy  cuts " << cuts_e << std::endl;
    std::cout << "Tel     cuts " << cuts_tel << std::endl;
    std::cout << "Azimuth cuts " << cuts_az << std::endl;
    std::cout << "Offset  cuts " << cuts_off << std::endl;
    std::cout << "Writing " << ntuple->GetEntries() << " entries." << std::endl;
    std::cout << "Writing to " << outpath.str() << "." << std::endl;
    ntuple->GetCurrentFile()->Write();
    g.Close();
}

int main(int argc, char** argv){
  std::string pathbase = argv[1];
  bool binned = 1;
  if(argc == 3) binned = 0;
  cacheData_vegas(pathbase, binned);
  return 0;
}
