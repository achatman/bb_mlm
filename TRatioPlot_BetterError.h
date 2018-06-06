#ifndef TRATIOPLOT_BETTERERROR_H
#define TRATIOPLOT_BETTERERROR_H


#include "TGraphAsymmErrors.h"
#include "TFitResult.h"

//Author: Paul Gessinger   25/08/2016
//Altered to account for errors on both histograms when using diffsig
class TRatioPlot_BetterError : public TRatioPlot{

private:
  enum CalculationMode {
    kDivideHist = 1, ///< Use `TH1::Divide` to create the ratio.
    kDivideGraph = 2, ///< Use `TGraphAsymmErrors::Divide` to create the ratio.
    kDifference = 3, ///< Calculate the difference between the histograms.
    kFitResidual = 4, ///< Calculate the fit residual between the histogram and a fit stored within it.
    kDifferenceSign = 5 ///< Calculate the difference divided by the error.
  };

  enum ErrorMode {
    kErrorSymmetric = 1, ///< Use the regular `TH1::GetBinError` as the error
    kErrorAsymmetric = 2, ///< Use `TH1::GetBinErrorUp` and `TH1::GetBinErrorLow` for the error, depending on y values.
    kErrorFunc = 3 ///< Use the square root of the function value as the error.
  };

  enum HideLabelMode {
    kHideUp = 1, ///< Hide the first label of the upper y axis when there is low space.
    kHideLow = 2, ///< Hide the last label of the lower y axis when there is low space.
    kNoHide = 3, ///< Do not hide labels when there is low space.
    kForceHideUp = 4, ///< Always hide the first label of the upper y axis
    kForceHideLow = 5 ///< Always hide the last label of the lower y axis
  };

protected:

  void Init(TH1* h1, TH1* h2, Option_t *option)
  {

    fH1 = h1;
    fH2 = h2;

    SetupPads();

    TString optionString = TString(option);

    if (optionString.Contains("divsym")) {
      optionString.ReplaceAll("divsym", "");
      fMode = CalculationMode::kDivideHist;
    } else if (optionString.Contains("diffsig")) {
      optionString.ReplaceAll("diffsig", "");
      fMode = CalculationMode::kDifferenceSign;

      // determine which error style
      if (optionString.Contains("errasym")) {
        fErrorMode = ErrorMode::kErrorAsymmetric;
        optionString.ReplaceAll("errasym", "");
      }

      if (optionString.Contains("errfunc")) {
        fErrorMode = ErrorMode::kErrorFunc;
        optionString.ReplaceAll("errfunc", "");
      }
    } else if (optionString.Contains("diff")) {
      optionString.ReplaceAll("diff", "");
      fMode = CalculationMode::kDifference;
    } else {
      fMode = CalculationMode::kDivideGraph; // <- default
    }

    fOption = optionString;


    fH1DrawOpt = "hist";
    fH2DrawOpt = "E";
    fGraphDrawOpt = "AP";


    // build ratio, everything is ready
    if (!BuildLowerPlot()) return;

    // taking x axis information from h1 by cloning it x axis
    fSharedXAxis = (TAxis*)(fH1->GetXaxis()->Clone());
    fUpYaxis = (TAxis*)(fH1->GetYaxis()->Clone());
    fLowYaxis = (TAxis*)(fRatioGraph->GetYaxis()->Clone());
  }


  ////////////////////////////////////////////////////////////////////////////////
  /// Build the lower plot according to which constructor was called, and
  /// which options were passed.
  Int_t BuildLowerPlot()
  {
    std::cout << "POINT_BUILD" << std::endl;
    std::cout << "Calculation mode: " << fMode << std::endl;
    // Clear and delete the graph if not exists
    if (fRatioGraph != 0) {
      fRatioGraph->IsA()->Destructor(fRatioGraph);
      fRatioGraph = 0;
    }

    if (fConfidenceInterval1 == 0) {
      fConfidenceInterval1 = new TGraphErrors();
    }

    if (fConfidenceInterval2 == 0) {
      fConfidenceInterval2 = new TGraphErrors();
    }

    static Double_t divideGridlines[] = {0.7, 1.0, 1.3};
    static Double_t diffGridlines[] = {0.0};
    static Double_t signGridlines[] = {1.0, 0.0, -1.0};

    // Determine the divide mode and create the lower graph accordingly
    // Pass divide options given in constructor
    if (fMode == CalculationMode::kDivideGraph) {
      // use TGraphAsymmErrors Divide method to create

      SetGridlines(divideGridlines, 3);

      TH1 *tmpH1 = (TH1*)fH1->Clone();
      TH1 *tmpH2 = (TH1*)fH2->Clone();

      tmpH1->Scale(fC1);
      tmpH2->Scale(fC2);

      TGraphAsymmErrors *ratioGraph = new TGraphAsymmErrors();
      ratioGraph->Divide(tmpH1, tmpH2, fOption.Data());
      fRatioGraph = ratioGraph;

      delete tmpH1;
      delete tmpH2;

    } else if (fMode == CalculationMode::kDifference) {
      SetGridlines(diffGridlines, 3);

      TH1 *tmpHist = (TH1*)fH1->Clone();

      tmpHist->Reset();

      tmpHist->Add(fH1, fH2, fC1, -1*fC2);
      fRatioGraph = new TGraphErrors(tmpHist);

      delete tmpHist;
    } else if (fMode == CalculationMode::kDifferenceSign) {

      SetGridlines(signGridlines, 3);

      fRatioGraph = new TGraphAsymmErrors();
      Int_t ipoint = 0;
      Double_t res;
      Double_t error;

      Double_t val;
      Double_t val2;

      for (Int_t i=0; i<=fH1->GetNbinsX();++i) {
        val = fH1->GetBinContent(i);
        val2 = fH2->GetBinContent(i);

        if (fErrorMode == ErrorMode::kErrorAsymmetric) {

          Double_t errUp = fH1->GetBinErrorUp(i);
          Double_t errLow = fH1->GetBinErrorLow(i);

          if (val - val2 > 0) {
            // h1 > h2
            error = errLow;
          } else {
            // h1 < h2
            error = errUp;
          }

        } else if (fErrorMode == ErrorMode::kErrorSymmetric) {
          error = TMath::Sqrt( TMath::Power(fH1->GetBinError(i), 2) + TMath::Power(fH2->GetBinError(i), 2) );
          std::cout << "POINT_ERROR" << std::endl;
        } else {
          Warning("BuildLowerPlot", "error mode is invalid");
          error = 0;
        }

        if (error != 0) {

          res = (val - val2) / error;

          ((TGraphAsymmErrors*)fRatioGraph)->SetPoint(ipoint, fH1->GetBinCenter(i), res);
          ((TGraphAsymmErrors*)fRatioGraph)->SetPointError(ipoint,  fH1->GetBinWidth(i)/2., fH1->GetBinWidth(i)/2., 0.5, 0.5);

          ++ipoint;

        }
      }

    } else if (fMode == CalculationMode::kFitResidual) {

      SetGridlines(signGridlines, 3);

      TF1 *func = dynamic_cast<TF1*>(fH1->GetListOfFunctions()->At(0));

      if (func == 0) {
        // this is checked in constructor and should thus not occur
        Error("BuildLowerPlot", "h1 does not have a fit function");
        return 0;
      }

      fRatioGraph = new TGraphAsymmErrors();
      Int_t ipoint = 0;

      Double_t res;
      Double_t error;

      std::vector<double> ci1;
      std::vector<double> ci2;

      Double_t x_arr[fH1->GetNbinsX()];
      std::fill_n(x_arr, fH1->GetNbinsX(), 0);
      Double_t ci_arr1[fH1->GetNbinsX()];
      std::fill_n(ci_arr1, fH1->GetNbinsX(), 0);
      Double_t ci_arr2[fH1->GetNbinsX()];
      std::fill_n(ci_arr2, fH1->GetNbinsX(), 0);
      for (Int_t i=0; i<fH1->GetNbinsX();++i) {
        x_arr[i] = fH1->GetBinCenter(i+1);
      }

      Double_t cl1 = fCl1;
      Double_t cl2 = fCl2;

      if (fFitResult != 0) {
        // use this to get conf int

        fFitResult->GetConfidenceIntervals(fH1->GetNbinsX(), 1, 1, x_arr, ci_arr1, cl1);
        for (Int_t i=1; i<=fH1->GetNbinsX();++i) {
          ci1.push_back(ci_arr1[i-1]);
        }

        fFitResult->GetConfidenceIntervals(fH1->GetNbinsX(), 1, 1, x_arr, ci_arr2, cl2);
        for (Int_t i=1; i<=fH1->GetNbinsX();++i) {
          ci2.push_back(ci_arr2[i-1]);
        }
      } else {
        (TVirtualFitter::GetFitter())->GetConfidenceIntervals(fH1->GetNbinsX(), 1, x_arr, ci_arr1, cl1);
        for (Int_t i=1; i<=fH1->GetNbinsX();++i) {
          ci1.push_back(ci_arr1[i-1]);
        }
        (TVirtualFitter::GetFitter())->GetConfidenceIntervals(fH1->GetNbinsX(), 1, x_arr, ci_arr2, cl2);
        for (Int_t i=1; i<=fH1->GetNbinsX();++i) {
          ci2.push_back(ci_arr2[i-1]);
        }

      }

      Double_t x;
      Double_t val;

      for (Int_t i=0; i<=fH1->GetNbinsX();++i) {
        val = fH1->GetBinContent(i);
        x = fH1->GetBinCenter(i+1);

        if (fErrorMode == ErrorMode::kErrorAsymmetric) {

          Double_t errUp = fH1->GetBinErrorUp(i);
          Double_t errLow = fH1->GetBinErrorLow(i);

          if (val - func->Eval(fH1->GetBinCenter(i)) > 0) {
            // h1 > fit
            error = errLow;
          } else {
            // h1 < fit
            error = errUp;
          }

        } else if (fErrorMode == ErrorMode::kErrorSymmetric) {
          error = fH1->GetBinError(i);
        } else if (fErrorMode == ErrorMode::kErrorFunc) {

          error = sqrt(func->Eval(x));

        } else {
          Warning("BuildLowerPlot", "error mode is invalid");
          error = 0;
        }

        if (error != 0) {

          res = (fH1->GetBinContent(i)- func->Eval(fH1->GetBinCenter(i) ) ) / error;
          //__("x="<< x << " y=" << res << " err=" << error);

          ((TGraphAsymmErrors*)fRatioGraph)->SetPoint(ipoint, fH1->GetBinCenter(i), res);
          ((TGraphAsymmErrors*)fRatioGraph)->SetPointError(ipoint,  fH1->GetBinWidth(i)/2., fH1->GetBinWidth(i)/2., 0.5, 0.5);

          fConfidenceInterval1->SetPoint(ipoint, x, 0);
          fConfidenceInterval1->SetPointError(ipoint, x, ci1[i] / error);
          fConfidenceInterval2->SetPoint(ipoint, x, 0);
          fConfidenceInterval2->SetPointError(ipoint, x, ci2[i] / error);

          ++ipoint;

        }

      }

    } else if (fMode == CalculationMode::kDivideHist){
      SetGridlines(divideGridlines, 3);

      // Use TH1's Divide method
      TH1 *tmpHist = (TH1*)fH1->Clone();
      tmpHist->Reset();

      tmpHist->Divide(fH1, fH2, fC1, fC2, fOption.Data());
      fRatioGraph = new TGraphErrors(tmpHist);

      delete tmpHist;
    } else {
      // this should not occur
      Error("BuildLowerPlot", "Invalid fMode value");
      return 0;
    }

    // need to set back to "" since recreation. we don't ever want
    // title on lower graph

    if (fRatioGraph == 0) {
      Error("BuildLowerPlot", "Error creating lower graph");
      return 0;
    }

    fRatioGraph->SetTitle("");
    fConfidenceInterval1->SetTitle("");
    fConfidenceInterval2->SetTitle("");

    return 1;
  }

public:
  TRatioPlot_BetterError(TH1* h1, TH1* h2, Option_t *option)
  {
    std::cout << "POINT_CREATE" << std::endl;
    gROOT->GetListOfCleanups()->Add(this);

    if (!h1 || !h2) {
      Warning("TRatioPlot", "Need two histograms.");
      return;
    }

    Bool_t h1IsTH1=h1->IsA()->InheritsFrom(TH1::Class());
    Bool_t h2IsTH1=h2->IsA()->InheritsFrom(TH1::Class());

    if (!h1IsTH1 && !h2IsTH1) {
      Warning("TRatioPlot", "Need two histograms deriving from TH2 or TH3.");
      return;
    }

    fHistDrawProxy = h1;

    Init(h1, h2, option);

  }

};


#endif