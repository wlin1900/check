import sys
import argparse
import pathlib
import ROOT
import json
import re

toolspath = pathlib.Path(__file__).parent.joinpath("../../tools/").resolve()
sys.path.append(str(toolspath))
import b2genvalues
import categorydef
import dmiddef
import pdgvalues
import read_ntuples


def arg_parser():

    _myUsage = " \n"
    _myUsage += " \n"
    _myUsage += "python3 " + sys.argv[0] + " \n"
    _myUsage += "        --inputROOT    [(str)   input_histogram_pdf.root] \n"
    _myUsage += "        --outputROOT   [(str)   output_histogram_pdf.root] \n"
    _myUsage += "        --fitLumi      [(float) LLL.L] \n"
    _myUsage += "        --x1VarName    [(str)   X1-VARIABLE-OBJECT-NAME] _(optional)_ \n"
    _myUsage += "        --x2VarName    [(str)   X2-VARIABLE-OBJECT-NAME] _(optional)_ \n"
    _myUsage += "        --wsName       [(str)   WORKSPACE-OBJECT-NAME]   _(optional)_ \n"
    _myUsage += "        --dataSetLabel [(str)   LABEL]                   _(optional)_ \n"
    _myUsage += "        --modelNum     [(int)   KKK]                     _(optional)_ \n"
    _myUsage += "        --inputUnfoldFactorJsonName_R    [(str)  inputUnfoldFactorJsonName_R]"
    _myUsage += "        --inputFactorJsonName [(str) inputFactorJsonName]                   \n"
    _myUsage += " \n"

    _parser = argparse.ArgumentParser(usage=_myUsage)
    _parser.add_argument("--inputROOT", default="",
                         help="Input ROOT file name containing histogram pdfs classified by candidate categories.")
    _parser.add_argument("--outputROOT", default="",
                         help="Output ROOT file name to save 1D histograms pdfs for the fit.")
    _parser.add_argument("--fitLumi", default=0., type=float,
                         help="Luminosity of (pseudo) data used for the fit in fb-1.")
    _parser.add_argument("--x1VarName", default="x1Var",
                         help="Object name of X1 variable in the two-dimensional fit.")
    _parser.add_argument("--x2VarName", default="x2Var",
                         help="Object name of X2 variable in the two-dimensional fit.")
    _parser.add_argument("--wsName", default="wsHistPdf",
                         help="Object name of RooWorkSpace")
    _parser.add_argument("--dataSetLabel", default="",
                         help="Label for name of RooDataSet.")
    _parser.add_argument("--modelNum", default=None, type=int,
                         help="Number for identification of the model. The input value will overwrite the original value.")
    _parser.add_argument("--inputUnfoldFactorJsonName_R", default="", 
                         help="Input the name of the file which store T_migration factors")
    _parser.add_argument("--inputFactorJsonName", default="",
                         help="inputFactorJsonName")
    input_checker(parser=_parser)

    return _parser.parse_args()


def input_checker(parser):

    _inROOTPath = pathlib.Path(parser.parse_args().inputROOT).resolve()

    print(" ")
    print("=" * 100)
    print(" * Input ROOT    |", parser.parse_args().inputROOT)
    print(" * Output ROOT   |", parser.parse_args().outputROOT)
    print("." * 100)
    print(" * RooWorkspace  |", parser.parse_args().wsName)
    print(" * DataSet label |", parser.parse_args().dataSetLabel)
    print(" * Model number  |", parser.parse_args().modelNum)
    print("=" * 100)
    print(" ")

    print("*" * 50)
    print("*****   Luminosity for fit =  {:^7.2f} fb-1   *****".format(parser.parse_args().fitLumi))
    print("*" * 50)
    print(" ")

    print("-" * 100)
    print(" Fit variables ")
    print("-" * 100)
    print(" * X1 .................... {} ".format(parser.parse_args().x1VarName))
    print(" * X2 .................... {} ".format(parser.parse_args().x2VarName))
    print("-" * 100)
    print(" ")

    _checker = not _inROOTPath.is_file()
    _checker += len(parser.parse_args().outputROOT) < 5
    _checker += parser.parse_args().fitLumi < 1.0
    _checker += len(parser.parse_args().x1VarName) < 1
    _checker += len(parser.parse_args().x2VarName) < 1
    _checker += len(parser.parse_args().wsName) < 1
    _checker += len(parser.parse_args().inputUnfoldFactorJsonName_R) < 1
    _checker += len(parser.parse_args().inputFactorJsonName) < 1
    if _checker:
        parser.error("\n[ERROR] Treminated by argument errors! \n")

    return

def check_dependencies(var):

    print(f"\n>>> Checking {var.GetName()} ({var.ClassName()}) = {var.getVal()}")
    servers = var.getParameters(ROOT.RooArgSet())
    for s in servers:
        try:
            val = s.getVal()
        except Exception:
            val = "N/A"
        print(f"  - {s.GetName():40s} {s.ClassName():20s} value = {val}")


class CombinePDFs:

    def __init__(self):

        self.define_empty_dictionaries()

        return

    def define_empty_dictionaries(self):

        self.inHistPdfDict = {}
        self.inHistPdfSigMCDict = {}
        self.sumEntriesVarDict = {}
        self.pdfDict = {}
        self.numVarDict = {}
        self.ratioVarDict = {}
        self.RDstVar = []
        self.RDVar = []
        self.unfold_numVarDict = {}
        self._temp_var_parts = [] 
        self._temp_var_parts1 = [] 

        return

    def define_fit_luminosity_variable(
                self,
                fitLumi=0.0,
            ):

        self.fitLumiVar = ROOT.RooConstVar(
            "fitLumiVar",
            "Fit luminosity",
            float(fitLumi),
        )
        self.fitLumiVar.setUnit("fb^{-1}")

        return

    def read_roofit_variables(
                self,
                inROOTName="",
                x1VarName="",
                x2VarName="",
                wsName="",
            ):

        _inROOTPath = pathlib.Path(inROOTName).resolve()
        self.inROOT = ROOT.TFile(str(_inROOTPath))
        self.inTree = self.inROOT.Get("tree")
        self.wsHistPdf = self.inROOT.Get(wsName)

        self.xVars = []
        self.xVars.append(self.wsHistPdf.var(x1VarName))
        self.xVars.append(self.wsHistPdf.var(x2VarName))
        self.xVarsList = ROOT.RooArgList(*self.xVars)
        self.xVarsSet = ROOT.RooArgSet(*self.xVars)
        self.wgtVar = self.wsHistPdf.var("wgtVar")
        self.q2Var = self.wsHistPdf.var("q2Var")
        self.isChargedBEvent = self.wsHistPdf.var("isChargedBEvent")

        self.gnrMCLumiVar = self.wsHistPdf.arg("gnrMCLumiVar")
        self.nGenSigMCVar = self.wsHistPdf.arg("nGenSigMCVar")
        self.fitCatNumVar = self.wsHistPdf.arg("fitCatNumVar")
        self.dataSetNumVar = self.wsHistPdf.arg("dataSetNumVar")
        self.modelNumVar = self.wsHistPdf.arg("modelNumVar")

        return

    def define_fit_parameter_variables(
                self,
            ):
        
        initRDsts = pdgvalues.PDG.initRDsts
        initRDs = pdgvalues.PDG.initRDs
        init_Vals = zip(initRDsts, initRDs)
        for i, (initRDst, initRD) in enumerate(init_Vals):
            RDst_name = f"RDstVar_{i+1}"
            RD_name = f"RDVar_{i+1}"

            RDstVar = ROOT.RooRealVar(
                RDst_name,
                "#it{R(D*)}",
                initRDst,
                0.0,
                10.0,
            )

            self.RDstVar.append(RDstVar)

            RDVar = ROOT.RooRealVar(
                RD_name,
                "#it{R(D)}",
                initRD,
                0.0,
                10.0,
            )

            self.RDVar.append(RDVar)

        return 

    def read_MC_pdfs(
                self,
                dataSetLabel="",
            ):

        print(" ")
        print("[INFO] Reading RooHistPdf objects ...... ")

        reader = read_ntuples.Read()
        n_bins = reader._get_bin_info()

        for _candCatName in categorydef.candCatNumDict.keys():
            self.inHistPdfDict[_candCatName] = {}
            for i in range(n_bins):
                _pdfName = "rooHistPdf"
                _pdfName += f"_bin{i+1}_MC"
                _pdfName += "_{}".format(_candCatName)
                _pdfName += "_{}".format(dataSetLabel) if len(dataSetLabel) > 0 else ""
                _inHistPdf = self.wsHistPdf.pdf(_pdfName)
                if _inHistPdf != None:
                    self.inHistPdfDict[_candCatName][i] = _inHistPdf

        self.inHistPdfDict["Total"] = {} 
        for i in range(n_bins):
            _candCatName = "Total"
            _pdfName = "rooHistPdf"
            _pdfName += f"_bin{i+1}_MC"
            _pdfName += "_{}".format(_candCatName)
            _pdfName += "_{}".format(dataSetLabel) if len(dataSetLabel) > 0 else ""
            _inHistPdf = self.wsHistPdf.pdf(_pdfName)
            if _inHistPdf != None:
                self.inHistPdfDict[_candCatName][i] = _inHistPdf

            for _i, _candCatName in enumerate(self.inHistPdfDict.keys()):
                print(" {:>2d}. {} ".format(_i + 1, _candCatName))
            print(" ")

        return

    def read_sum_entries_variables(self):

        for _candCatName, _histPdfs in self.inHistPdfDict.items():
            self.sumEntriesVarDict[_candCatName] = {}
            for i, _histPdf in _histPdfs.items():
                _varName = _histPdf.GetName().replace("rooHistPdf", "sumEntries")
                _sumEntriesVar = self.wsHistPdf.arg(_varName)
                self.sumEntriesVarDict[_candCatName][i] = _sumEntriesVar

        return

    def combine_pdf_objects(self):

        for _candCatName, _histPdfs in self.inHistPdfDict.items():
            self.pdfDict[_candCatName] = {}
            if _candCatName == "Total":
                continue
            for i, _histPdf in _histPdfs.items():
                _pdfName = _histPdf.GetName().replace("rooHistPdf", "pdf")
                self.pdfDict[_candCatName][i] = _histPdf.Clone(_pdfName)

        return

    def define_reconstruction_efficiency_variables(self):

        _sigCandCatName = "SignalDstTv"
        _normCandCatName = "NormDstlv"
        self.effVarSig_Dst = self.define_efficiency_variable(candCatName=_sigCandCatName)
        self.effVarNorm_Dst = self.define_efficiency_variable(candCatName=_normCandCatName)

        _sigCandCatName = "SignalDTv"
        _normCandCatName = "NormDlv"
        self.effVarSig_D = self.define_efficiency_variable(candCatName=_sigCandCatName)
        self.effVarNorm_D = self.define_efficiency_variable(candCatName=_normCandCatName)

        return

    def define_efficiency_variable(
                self,
                candCatName="",
                nominalVal=1.0 * (10.0 ** -4.0),
            ):

        #_name = self.pdfDict[candCatName].GetName().replace("pdf_MC", "effVar")
        num = len(self.pdfDict[candCatName])
        _effVar = {}
        for i in range(num): 
            _name = self.pdfDict[candCatName][i].GetName().replace("pdf", "effVar")
            _candCatNameIndex = categorydef.convert_candidate_category_to_legend_string(candCatName=candCatName)
            _fitCatNum = int(self.fitCatNumVar.getVal())
            _title = "#it{#epsilon}_{CANDCAT}^{DECAYMODE}"
            _title = _title.replace("CANDCAT", _candCatNameIndex)
            _title = _title.replace("DECAYMODE", dmiddef.strXcDecDict[_fitCatNum])
            _effVar[i+1] = ROOT.RooRealVar(
                _name,
                _title,
                nominalVal,
                0.0,
                1.0,
            )
        #_effVar = ROOT.RooConstVar(
        #    _name,
        #    _title,
        #    nominalVal,
        #)
        return _effVar

    def define_floated_fit_parameters(self):

        for _candCatName in categorydef.floatedCandCatNames:
            _pdfs = self.pdfDict.get(_candCatName)
            if _pdfs is None:
                continue
            self.numVarDict[_candCatName] = {}
            for i, _pdf in _pdfs.items():
                _name = _pdf.GetName().replace("pdf", "numVar")
                _candCatNameIndex = categorydef.convert_candidate_category_to_legend_string(candCatName=_candCatName)
                _fitCatNum = int(self.fitCatNumVar.getVal())
                _title = "#it{N}_{CANDCAT}^{DECAYMODE}"
                _title = _title.replace("CANDCAT", _candCatNameIndex)
                _title = _title.replace("DECAYMODE", dmiddef.strXcDecDict[_fitCatNum])
                _numVal = self.sumEntriesVarDict[_candCatName][i].getVal()
                _numVal *= (self.fitLumiVar.getVal() / self.gnrMCLumiVar.getVal())
                _totalNumVal = self.sumEntriesVarDict["Total"][i].getVal()
                if _numVal == 0.0:
                    continue
                _numVar = ROOT.RooRealVar(
                    _name,
                    _title,
                    _numVal,
                    _totalNumVal * -10.0,
                    _totalNumVal * 10.0,
                )
                self.numVarDict[_candCatName][i] = _numVar

        return

    def define_fixed_parameters(self):

        for _candCatName in categorydef.fixedCandCatNames:
            _pdfs = self.pdfDict.get(_candCatName)
            if _pdfs is None:
                continue
            self.numVarDict[_candCatName] = {}
            for i, _pdf in _pdfs.items():
                _name = _pdf.GetName().replace("pdf", "numVar")
                _candCatNameIndex = categorydef.convert_candidate_category_to_legend_string(candCatName=_candCatName)
                _fitCatNum = int(self.fitCatNumVar.getVal())
                _title = "#it{N}_{CANDCAT}^{DECAYMODE}"
                _title = _title.replace("CANDCAT", _candCatNameIndex)
                _title = _title.replace("DECAYMODE", dmiddef.strXcDecDict[_fitCatNum])
                _numVal = self.sumEntriesVarDict[_candCatName][i].getVal()
                _numVal *= (self.fitLumiVar.getVal() / self.gnrMCLumiVar.getVal())
                _numVar = ROOT.RooConstVar(
                    _name,
                    _title,
                    _numVal,
                )
                self.numVarDict[_candCatName][i] = _numVar

        return

    #Step 1. N(B->D(*)lnv)_{gen}
    def define_branching_ratio_fit_parameters(self, _inputFactorJsonName=""):

        self.f00Var = ROOT.RooRealVar(
                "f00Var",
                "#it{f}_{00}",
                pdgvalues.PDG.f00,
                0.0,
                1.0,
            )

        self.fpmVar = ROOT.RooRealVar(
                "fpmVar",
                "#it{f}_{pm}",
                pdgvalues.PDG.fpm,
                0.0,
                1.0,
            )

        reader = read_ntuples.Read()
        corr_factors = reader.read_factors(_inputFactorJsonName)

        _NBBVal = b2genvalues.nGenBBbar1batch / b2genvalues.lumi1batch
        _NBBVal *= self.fitLumiVar.getVal()
        self.NBBbarVar = ROOT.RooRealVar(
                "NBBbarVar",
                "Number of #it{B#bar{Bbar}} events",
                _NBBVal,
                0.0,
                _NBBVal * 10.0,
            )

        _frac_Dst_v0_Val = b2genvalues.brB02Dstlnu / (b2genvalues.brB02Dstlnu + b2genvalues.brBp2Dstlnu)
        self.frac_Dst_b0 = ROOT.RooConstVar(
            "frac_Dst_b0",
            "fraction of B0 in Dstlv",
            _frac_Dst_v0_Val
        )
        _frac_D_v0_Val = b2genvalues.brB02Dlnu / (b2genvalues.brB02Dlnu + b2genvalues.brBp2Dlnu)
        self.frac_D_b0 = ROOT.RooConstVar(
            "frac_D_b0",
            "fraction of B0 in Dlv",
            _frac_D_v0_Val
        )
        
        self.f_bin_NormDlv_b0 = {}
        self.f_bin_NormDlv_bp = {}
        self.f_bin_NormDstlv_b0 = {}
        self.f_bin_NormDstlv_bp = {}

        _candCatName = "NormDstlv"
        _name = "brVar_Norm_BDstlv"
        _title = "#it{Br}(#it{#bar{B}}#rightarrow#it{D*l{#nu}})"
        _brVal = (b2genvalues.brB02Dstlnu + b2genvalues.brBp2Dstlnu) / 2.0
        _strFormula = "(@0 * 2.0) * (@2 * @1 * @5 + (1 - @2) * @3 * @6) * (@4 * 2.0) * @7"

        self.brVarNorm_Dst = ROOT.RooRealVar(
                _name,
                _title,
                _brVal,
                0.0,
                1.0,
            )

        _pdfs = self.pdfDict.get(_candCatName)  # measured bins dict
        self.unfold_numVarDict[_candCatName] = {}

        for i, _pdf in _pdfs.items():

            truth_bin_index = i + 1

            f_bin_val = corr_factors[_candCatName][truth_bin_index]
            self.f_bin_NormDstlv_b0[truth_bin_index] = ROOT.RooConstVar(
                    f"f_bin_{_candCatName}_{truth_bin_index}_b0", "bin correction factor B0", f_bin_val["0"]
                )
            self.f_bin_NormDstlv_bp[truth_bin_index] = ROOT.RooConstVar(
                    f"f_bin_{_candCatName}_{truth_bin_index}_bp", "bin correction factor Bp", f_bin_val["1"]
                )

            _name = _pdf.GetName().replace("pdf", "unfold_numVar")
            _title = "#it{unfold_numVar}(#it{#bar{B}}#rightarrow#it{D*l{#nu}})"

            _argList = ROOT.RooArgList(
                    self.brVarNorm_Dst,
                    self.f00Var,
                    self.frac_Dst_b0,
                    self.fpmVar,
                    self.NBBbarVar,
                    self.f_bin_NormDstlv_b0[truth_bin_index],
                    self.f_bin_NormDstlv_bp[truth_bin_index],
                    self.effVarNorm_Dst[truth_bin_index],
                )

            _numVar = ROOT.RooFormulaVar(
                    _name,
                    _title,
                    _strFormula,
                    _argList,
                )
            self.unfold_numVarDict[_candCatName][truth_bin_index] = _numVar
            
        _name = self.unfold_numVarDict[_candCatName][1].GetName()
        match = re.search(r"(Xc.*)$", _name)
        unfold_numVar_NormDstlv_bin0 = 2223.75
        self.unfold_numVarDict[_candCatName][0] = ROOT.RooConstVar(
                f"unfold_numVar_bin0_MC_NormDstlv_{match.group(1)}",
                "fixed MC template yield for truth first bin",
                unfold_numVar_NormDstlv_bin0
            )

        _candCatName = "NormDlv"
        _name = "brVar_Norm_BDlv"
        _title = "#it{Br}(#it{#bar{B}#rightarrow#it{D}l{#nu}})"
        _brVal = (b2genvalues.brB02Dlnu + b2genvalues.brBp2Dlnu) / 2.0
        _strFormula = "(@0 * 2.0) * (@2 * @1 * @5 + (1 - @2) * @3 * @6) * (@4 * 2.0) * @7"

        self.brVarNorm_D = ROOT.RooRealVar(
                _name,
                _title,
                _brVal,
                0.0,
                1.0,
            )

        _pdfs = self.pdfDict.get(_candCatName)
        self.unfold_numVarDict[_candCatName] = {}

        for i, _pdf in _pdfs.items():
            truth_bin_index = i + 1  # measured bin i ↔ truth bin i+1

            f_bin_val = corr_factors[_candCatName][truth_bin_index]
            self.f_bin_NormDlv_b0[truth_bin_index] = ROOT.RooConstVar(
                    f"f_bin_{_candCatName}_{truth_bin_index}_b0", "bin correction factor B0", f_bin_val["0"]
                )
            self.f_bin_NormDlv_bp[truth_bin_index] = ROOT.RooConstVar(
                    f"f_bin_{_candCatName}_{truth_bin_index}_bp", "bin correction factor Bp", f_bin_val["1"]
                )

            _name = _pdf.GetName().replace("pdf", "unfold_numVar")
            _title = "#it{unfold_numVar}(#it{#bar{B}}#rightarrow#it{Dl{#nu}})"

            _argList = ROOT.RooArgList(
                    self.brVarNorm_D,
                    self.f00Var,
                    self.frac_D_b0,
                    self.fpmVar,
                    self.NBBbarVar,
                    self.f_bin_NormDlv_b0[truth_bin_index],
                    self.f_bin_NormDlv_bp[truth_bin_index],
                    self.effVarNorm_D[truth_bin_index],
                )

            _numVar = ROOT.RooFormulaVar(
                    _name,
                    _title,
                    _strFormula,
                    _argList,
                )
            self.unfold_numVarDict[_candCatName][truth_bin_index] = _numVar

        _name = self.unfold_numVarDict[_candCatName][1].GetName()
        match = re.search(r"(Xc.*)$", _name)
        unfold_numVar_NormDlv_bin0 = 113.08
        self.unfold_numVarDict[_candCatName][0] = ROOT.RooConstVar(
                f"unfold_numVar_bin0_MC_NormDlv_{match.group(1)}",
                "fixed MC template yield for truth first bin",
                unfold_numVar_NormDlv_bin0
            )

        return

    # Step 2. N(B->D(*)lnv)_{measured}
    def define_norm_fit_parameters(self, R_Migrate=""):

        with open(R_Migrate, 'r') as f:
            R_Migrate = json.load(f)

        self.s_bin0 = ROOT.RooRealVar(
            "s_bin0",
            "scale factor for truth bin0 template",
            1.0,
            0.5,
            2
        )

        for _candCatName in ['NormDstlv', 'NormDlv']:

            self.numVarDict[_candCatName] = {}

            for meas_bin_index in self.unfold_numVarDict[_candCatName].keys():

                if meas_bin_index == 3:
                    continue
                _argList = ROOT.RooArgList()
                var_parts = []

                _name = self.unfold_numVarDict[_candCatName][meas_bin_index].GetName()
                match = re.search(r"(Xc.*)$", _name)
                Migrate_factor_name = "{}_{}".format(_candCatName, match.group(1))

                for truth_bin_index in self.unfold_numVarDict[_candCatName].keys():

                    if truth_bin_index == 0:
                        continue

                    Migrate_factor = R_Migrate[Migrate_factor_name][str(meas_bin_index)][str(truth_bin_index)]
                    _N_truth_i = self.unfold_numVarDict[_candCatName][truth_bin_index]

                    _name_part = _N_truth_i.GetName().replace(
                        f"unfold_numVar_{truth_bin_index}",
                        f"numVar_part_meas{meas_bin_index}_truth{truth_bin_index}"
                    )

                    _var_part = ROOT.RooFormulaVar(
                        _name_part,
                        f"Migrated contribution truth bin {truth_bin_index} → measured bin {meas_bin_index}",
                        f"@0 * {float(Migrate_factor)}",
                        ROOT.RooArgList(_N_truth_i),
                    )

                    var_parts.append(_var_part)
                    _argList.add(_var_part)


                N0_const = self.unfold_numVarDict[_candCatName][0]
                Migrate_factor_bin0 = R_Migrate[Migrate_factor_name][str(meas_bin_index)]["0"]

                T_j_val = float(Migrate_factor_bin0) * N0_const.getVal()

                T_j = ROOT.RooConstVar(
                        f"T_bin0_to_meas{meas_bin_index}_{_candCatName}",
                        f"template contribution of bin0 to measured bin{meas_bin_index}",
                        T_j_val
                    )

                sT_var = ROOT.RooFormulaVar(
                        f"sT_bin0_to_meas{meas_bin_index}_{_candCatName}",
                        f"s x T_j for measured bin {meas_bin_index}",
                        "@0*@1",
                        ROOT.RooArgList(self.s_bin0, T_j)
                    )

                var_parts.append(sT_var)
                _argList.add(sT_var)

                # ------------- Step 3.μ_j(θ, s) -------------
                sum_expr = " + ".join([f"@{k}" for k in range(_argList.getSize())])

                _name_total = _name.replace("unfold_numVar", f"numVar")
                _title_total = f"Total expected yield in measured bin {meas_bin_index}"

                _totalVar = ROOT.RooFormulaVar(
                    _name_total,
                    _title_total,
                    sum_expr,
                    _argList
                )

                print(f"Test---------------------------------->{_candCatName}---->{meas_bin_index}")
                print(_totalVar)

                self.numVarDict[_candCatName][meas_bin_index] = _totalVar

        return

    def define_q2_integration(self):

            reader = read_ntuples.Read()
            #n_bins = reader._get_bin_info()
            q2_bins = [3.2, 4, 6, 7.5, 12]
            n_bins = len(q2_bins)

            self.const_C2 = []

            for i in range(n_bins-1):

                q2_min = q2_bins[i]
                q2_max = q2_bins[i+1]

                q2_lower_limit = self.q2Var.getMin()
                q2_upper_limit = self.q2Var.getMax()

                actual_q2_min = max(q2_min, q2_lower_limit)
                actual_q2_max = min(q2_max, q2_upper_limit)

                range_name = f"q2_bin{i}"
                self.q2Var.setRange(range_name, actual_q2_min, actual_q2_max)

                #_numVal1 = actual_q2_max - actual_q2_min if actual_q2_max > actual_q2_min else 0.
                _numVal1 = 1

                const_C2 = ROOT.RooConstVar(
                    f"C2_bin{i}",
                    f"Integration of RD_C2",
                    _numVal1
                )

                self.const_C2.append(const_C2)

            return

    def define_q2_integration_RDst(self):
            
            self.mtau = ROOT.RooConstVar("mtau", "mtau", pdgvalues.PDG.getMass(15))

            reader = read_ntuples.Read()
            #n_bins = reader._get_bin_info()
            q2_bins = [3.2, 4, 6, 7.5, 12]
            n_bins = len(q2_bins)

            self.const_C1_RDst = []

            mtau = self.mtau.getVal() 

            for i in range(n_bins-1):
                q2_min, q2_max = q2_bins[i], q2_bins[i + 1]
                q2_center = 0.5 * (q2_min + q2_max)

                factor = (1 - mtau * mtau / q2_center) ** 2

                const_C1 = ROOT.RooConstVar(
                    f"C1_RDst_bin{i}",
                    f"Integration of 1/RDst_C1 (bin {i})",
                    factor
                )
                self.const_C1_RDst.append(const_C1)

    def define_q2_integration_RD(self):
            
            mass_Bplus = ROOT.RooConstVar("mass_Bplus", "mass_Bplus", pdgvalues.PDG.getMass(521))
            mass_Bzero = ROOT.RooConstVar("mass_Bzero", "mass_Bzero", pdgvalues.PDG.getMass(511))
            mass_Dzero = ROOT.RooConstVar("mass_Dzero", "mass_Dzero", pdgvalues.PDG.getMass(421))
            mass_Dplus = ROOT.RooConstVar("mass_Dplus", "mass_Dplus", pdgvalues.PDG.getMass(411))
            self.mtau = ROOT.RooConstVar("mtau", "mtau", pdgvalues.PDG.getMass(15))

            self.mB = ROOT.RooFormulaVar(
                "mB", "mB",
                "@0 * @1 + (1 - @0) * @2",
                ROOT.RooArgList(self.isChargedBEvent, mass_Bplus, mass_Bzero)
            )
            self.mD = ROOT.RooFormulaVar(
                "mD", "mD",
                "@0 * @1 + (1 - @0) * @2",
                ROOT.RooArgList(self.isChargedBEvent, mass_Dzero, mass_Dplus)
            )

            self.denom = ROOT.RooFormulaVar(
                "denom", "denom",
                "pow((@0*@0 - @1*@1), 2)",
                ROOT.RooArgList(self.mB, self.mD)
            )

            reader = read_ntuples.Read()
            #n_bins = reader._get_bin_info()
            q2_bins = [3.2, 4, 6, 7.5, 12]
            n_bins = len(q2_bins)

            mtau = self.mtau.getVal() 
            mB = self.mB.getVal() 
            mD = self.mD.getVal() 
            denom = self.denom.getVal()

            self.const_C1_RD = []

            for i in range(n_bins-1):
                q2_min, q2_max = q2_bins[i], q2_bins[i + 1]
                q2_center = 0.5 * (q2_min + q2_max)

                kinematic_factor = (1 - mtau * mtau / q2_center) * (1 - mtau * mtau / q2_center)

                denom_factor = (
                    (mB - mD) * (mB - mD) - q2_center
                ) * (
                    (mB + mD) * (mB + mD) - q2_center
                )

                numVal1 = kinematic_factor * denom / denom_factor

                const_C1 = ROOT.RooConstVar(
                    f"C1_RD_bin{i}",
                    f"Integration of 1/RD_C1 (bin {i})",
                    numVal1
                )
                self.const_C1_RD.append(const_C1)

    #Step 1. N(B->D(*)Tnv)_{gen}
    def define_unfold_signal_fit_parameters(self):

        _candCatName = "SignalDstTv"
        self.unfold_numVarDict[_candCatName] = {}

        _name = self.unfold_numVarDict["NormDstlv"][1].GetName()
        match = re.search(r"(Xc.*)$", _name)

        for i in self.unfold_numVarDict["NormDstlv"].keys():
            if i == 0:  
                continue

            _name = f"unfold_numVar_bin{i}_MC_{_candCatName}_{match.group(1)}"
            _title = f"Unfolded expected signal yield for bin {i} ({_candCatName})"

            _argList = ROOT.RooArgList(
                self.RDstVar[i-1],                        
                self.unfold_numVarDict["NormDstlv"][i],   
                self.effVarSig_Dst[i],                   
                self.effVarNorm_Dst[i],                   
                self.const_C1_RDst[i],                    
                self.const_C2[i],                         
            )

            _numVar = ROOT.RooFormulaVar(
                _name,
                _title,
                "2.0 * @0 * @1 * (@2 / @3) * @4 * @5",
                _argList,
            )
            self.unfold_numVarDict[_candCatName][i] = _numVar

        unfold_numVar_SignalDstTv_bin0 = 6.52
        self.unfold_numVarDict[_candCatName][0] = ROOT.RooConstVar(
                f"unfold_numVar_bin0_MC_SignalDstTv_{match.group(1)}",
                "fixed MC template yield for truth first bin",
                unfold_numVar_SignalDstTv_bin0
            )


        _candCatName = "SignalDTv"
        self.unfold_numVarDict[_candCatName] = {}

        for i in self.unfold_numVarDict["NormDlv"].keys():
            if i == 0:
                continue

            _name = f"unfold_numVar_bin{i}_MC_{_candCatName}_{match.group(1)}"
            _title = f"Unfolded expected signal yield for bin {i} ({_candCatName})"

            _argList = ROOT.RooArgList(
                self.RDVar[i-1],
                self.unfold_numVarDict["NormDlv"][i],
                self.effVarSig_D[i],
                self.effVarNorm_D[i],
                self.const_C1_RD[i],
                self.const_C2[i],
            )

            _numVar = ROOT.RooFormulaVar(
                _name,
                _title,
                "2.0 * @0 * @1 * (@2 / @3) * @4 * @5",
                _argList,
            )
            self.unfold_numVarDict[_candCatName][i] = _numVar

        unfold_numVar_SignalDTv_bin0 = 0
        self.unfold_numVarDict[_candCatName][0] = ROOT.RooConstVar(
                f"unfold_numVar_bin0_MC_SignalDTv_{match.group(1)}",
                "fixed MC template yield for truth first bin",
                unfold_numVar_SignalDTv_bin0
            )

        return

    # Step 2. N(B->D(*)Tnv)_{measured}
    def define_signal_fit_parameters(self, R_Migrate=""):

        with open(R_Migrate, 'r') as f:
            R_Migrate = json.load(f)

        # bin0 缩放参数
        self.s_bin0 = ROOT.RooRealVar(
            "s_bin0_signal",
            "scale factor for truth bin0 template",
            1.0,
            0.5,
            2.0
        )

        for _candCatName in ['SignalDstTv', 'SignalDTv']:
            self.numVarDict[_candCatName] = {}

            for meas_bin_index in self.unfold_numVarDict[_candCatName].keys():

                if meas_bin_index == 3:
                    continue

                _argList = ROOT.RooArgList()
                var_parts = []

                _name = self.unfold_numVarDict[_candCatName][meas_bin_index].GetName()
                match = re.search(r"(Xc.*)$", _name)
                Migrate_factor_name = "{}_{}".format(_candCatName, match.group(1))

                for truth_bin_index in self.unfold_numVarDict[_candCatName].keys():

                    if truth_bin_index == 0:
                        continue

                    Migrate_factor = R_Migrate[Migrate_factor_name][str(meas_bin_index)][str(truth_bin_index)]
                    _N_truth_i = self.unfold_numVarDict[_candCatName][truth_bin_index]

                    _name_part = _N_truth_i.GetName().replace(
                        f"unfold_numVar_{truth_bin_index}",
                        f"numVar_part_meas{meas_bin_index}_truth{truth_bin_index}"
                    )

                    _var_part = ROOT.RooFormulaVar(
                        _name_part,
                        f"Migrated contribution truth bin {truth_bin_index} → measured bin {meas_bin_index}",
                        f"@0 * {float(Migrate_factor)}",
                        ROOT.RooArgList(_N_truth_i)
                    )

                    var_parts.append(_var_part)
                    _argList.add(_var_part)


                N0_const = self.unfold_numVarDict[_candCatName][0]
                Migrate_factor_bin0 = R_Migrate[Migrate_factor_name][str(meas_bin_index)]["0"]

                T_j_val = float(Migrate_factor_bin0) * N0_const.getVal()
                T_j = ROOT.RooConstVar(
                        f"T_bin0_to_meas{meas_bin_index}_{_candCatName}",
                        f"template contribution of bin0 to measured bin {meas_bin_index}",
                        T_j_val
                    )

                sT_var = ROOT.RooFormulaVar(
                        f"sT_bin0_to_meas{meas_bin_index}_{_candCatName}",
                        f"s x T_j for measured bin {meas_bin_index}",
                        "@0*@1",
                        ROOT.RooArgList(self.s_bin0, T_j)
                    )

                var_parts.append(sT_var)
                _argList.add(sT_var)

                # μ_j(θ, s)
                sum_expr = " + ".join([f"@{k}" for k in range(_argList.getSize())])
                _name_total = _name.replace("unfold_numVar", f"numVar")
                _title_total = f"Total expected yield in measured bin {meas_bin_index}"

                _totalVar = ROOT.RooFormulaVar(
                    _name_total,
                    _title_total,
                    sum_expr,
                    _argList
                )
                print(f"Test---------------------------------->{_candCatName}---->{meas_bin_index}")
                print(_totalVar)

                self.numVarDict[_candCatName][meas_bin_index] = _totalVar

        return

    def define_yield_ratio_parameter(
                    self,
                    numeratorCandCatName="",
                    denominatorCandCatName="",
                ):
            
            if self.numVarDict.get(numeratorCandCatName) is None:
                return
            
            _name = self.numVarDict[numeratorCandCatName].GetName().replace("numVar", "ratioVar")
            _title = self.numVarDict[numeratorCandCatName].GetTitle().replace("#it{N}", "#it{r}")
            _ratioVal = self.sumEntriesVarDict[numeratorCandCatName].getVal()
            _ratioVal /= self.sumEntriesVarDict[denominatorCandCatName].getVal()
            _ratioVar = ROOT.RooConstVar(
                _name,
                _title,
                float(_ratioVal),
            )
            self.ratioVarDict[numeratorCandCatName] = _ratioVar

            _name = self.numVarDict[numeratorCandCatName].GetName()
            _title = self.numVarDict[numeratorCandCatName].GetTitle()
            _numVar = ROOT.RooFormulaVar(
                _name,
                _title,
                "@0 * @1",
                ROOT.RooArgList(self.numVarDict[denominatorCandCatName], _ratioVar),
            )
            self.numVarDict[numeratorCandCatName] = _numVar

            return

    def create_fit_model(self, dataSetLabel=""):
            _name_prefix = f"model_{dataSetLabel}" if dataSetLabel else "model"
            self.models = {}

            all_bins = set()
            for candCatName, pdf_bins in self.pdfDict.items():
                if isinstance(pdf_bins, dict):
                    all_bins.update(pdf_bins.keys())

            for i in sorted(all_bins):
                _pdfs = []
                _numVars = []

                for candCatName in self.pdfDict:
                    pdf_bins = self.pdfDict.get(candCatName, {})
                    numVar_bins = self.numVarDict.get(candCatName, {})

                    if not isinstance(pdf_bins, dict) or not isinstance(numVar_bins, dict):
                        continue

                    _pdf = pdf_bins.get(i)
                    _numVar = numVar_bins.get(i)

                    if _pdf is None or _numVar is None:
                        continue

                    _pdfs.append(_pdf)
                    _numVars.append(_numVar)

                if not _pdfs or not _numVars:
                    continue

                model_name = f"{_name_prefix}_bin{i+1}"
                self.models[i] = ROOT.RooAddPdf(
                    model_name,
                    model_name,
                    ROOT.RooArgList(*_pdfs),
                    ROOT.RooArgList(*_numVars),
                )

            return

    def overwrite_model_number(
                    self,
                    modelNum=None,
                ):

            if self.modelNumVar == None and modelNum is None:
                self.modelNumVar = ROOT.RooConstVar(
                    "modelNumVar",
                    "#Model",
                    0.0,
                )
            elif self.modelNumVar == None and modelNum is not None:
                self.modelNumVar = ROOT.RooConstVar(
                    "modelNumVar",
                    "#Model",
                    float(modelNum),
                )
            elif self.modelNumVar != None and modelNum is not None:
                self.modelNumVar.changeVal(float(modelNum))

            return

    def output_roofit_objects(
                    self,
                    outROOTPath="",
                ):

            self.wsPdf = ROOT.RooWorkspace(
                "wsPdf",
                "Workspace for pdf objects of RooHistPdf",
            )
            for _xVar in self.xVars:
                self.wsPdf.Import(_xVar)
            for RDstVar in self.RDstVar:
                self.wsPdf.Import(RDstVar)
            for RDVar in self.RDVar:
                self.wsPdf.Import(RDVar)
            self.wsPdf.Import(self.wgtVar)
            self.wsPdf.Import(self.q2Var)
            self.wsPdf.Import(self.isChargedBEvent)
            self.wsPdf.Import(self.fitLumiVar)
            self.wsPdf.Import(self.gnrMCLumiVar)
            self.wsPdf.Import(self.nGenSigMCVar)
            self.wsPdf.Import(self.fitCatNumVar)
            if self.dataSetNumVar != None:
                self.wsPdf.Import(self.dataSetNumVar)
            self.wsPdf.Import(self.modelNumVar)

            for i, effVarSig_Dst in self.effVarSig_Dst.items():
                self.wsPdf.Import(effVarSig_Dst)
            for i, effVarNorm_Dst in self.effVarNorm_Dst.items():
                self.wsPdf.Import(effVarNorm_Dst)
            for i, effVarSig_D in self.effVarSig_D.items():
                self.wsPdf.Import(effVarSig_D)
            for i, effVarNorm_D in self.effVarNorm_D.items():
                self.wsPdf.Import(effVarNorm_D)
            for _candCatName, _ratioVar in self.ratioVarDict.items():
                self.wsPdf.Import(_ratioVar)
            self.wsPdf.Import(self.f00Var)
            self.wsPdf.Import(self.fpmVar)
            self.wsPdf.Import(self.NBBbarVar)
            self.wsPdf.Import(self.brVarNorm_Dst)
            self.wsPdf.Import(self.brVarNorm_D)
            for f_bin_var in self.f_bin_NormDlv_b0.values():
                self.wsPdf.Import(f_bin_var)
            for f_bin_var in self.f_bin_NormDlv_bp.values():
                self.wsPdf.Import(f_bin_var)
            for f_bin_var in self.f_bin_NormDstlv_b0.values():
                self.wsPdf.Import(f_bin_var)
            for f_bin_var in self.f_bin_NormDstlv_bp.values():
                self.wsPdf.Import(f_bin_var)

            for const_C1_RD in self.const_C1_RD:
                self.wsPdf.Import(const_C1_RD)
            for const_C1_RDst in self.const_C1_RDst:
                self.wsPdf.Import(const_C1_RDst)
            for const_C2 in self.const_C2:
                self.wsPdf.Import(const_C2)

            for i, model in self.models.items():
                self.wsPdf.Import(model)
            for _candCatName, _numVarDict in self.numVarDict.items():
                for i, numVarDict in _numVarDict.items():
                    self.wsPdf.Import(numVarDict)

            for _candCatName, _unfold_numVarDict in self.unfold_numVarDict.items():
                for i, unfold_numVarDict in _unfold_numVarDict.items():
                    self.wsPdf.Import(unfold_numVarDict)
            
                    
            self.wsPdf.Print()
            self.wsPdf.writeToFile(str(outROOTPath))

            return

    def setup_output_ROOT_file(
                    self,
                    outROOTPath="",
                ):

            self.outROOT = ROOT.TFile(str(outROOTPath), "update")
            self.outTree = self.inTree.CopyTree("1")

            return

    def output_objects(
                    self,
                    outROOTName="",
                ):

            _outROOTPath = pathlib.Path(outROOTName).resolve()
            self.output_roofit_objects(outROOTPath=_outROOTPath)
            self.setup_output_ROOT_file(outROOTPath=_outROOTPath)
            self.outROOT.cd()
            self.outTree.Write()
            self.outROOT.Close()

            print(" ")
            print("=" * 100)
            print(" * Output ROOT file |", str(_outROOTPath))
            print("=" * 100)
            print(" * File size ........ {:.3f} MB".format(_outROOTPath.stat().st_size / (1024 ** 2)))
            print(" ")

            return

def combine_pdfs(
            inROOTName="",
            outROOTName="",
            fitLumi=0.0,
            x1VarName="",
            x2VarName="",
            wsName="",
            dataSetLabel="",
            modelNum=None,
            inputUnfoldFactorJsonName_R="",
            _inputFactorJsonName="",
        ):

    pdfCombiner = CombinePDFs()
    pdfCombiner.define_fit_luminosity_variable(fitLumi=fitLumi)
    pdfCombiner.read_roofit_variables(
        inROOTName=inROOTName,
        x1VarName=x1VarName,
        x2VarName=x2VarName,
        wsName=wsName,
    )
    pdfCombiner.define_fit_parameter_variables()
    pdfCombiner.read_MC_pdfs(dataSetLabel=dataSetLabel)
    pdfCombiner.read_sum_entries_variables()
    pdfCombiner.combine_pdf_objects()
    pdfCombiner.define_reconstruction_efficiency_variables()
    pdfCombiner.define_q2_integration()
    pdfCombiner.define_q2_integration_RD()
    pdfCombiner.define_q2_integration_RDst()
    pdfCombiner.define_floated_fit_parameters()
    pdfCombiner.define_fixed_parameters()
    pdfCombiner.define_branching_ratio_fit_parameters(_inputFactorJsonName=_inputFactorJsonName)
    pdfCombiner.define_norm_fit_parameters(R_Migrate=inputUnfoldFactorJsonName_R)
    pdfCombiner.define_unfold_signal_fit_parameters()
    pdfCombiner.define_signal_fit_parameters(R_Migrate=inputUnfoldFactorJsonName_R)
    pdfCombiner.define_yield_ratio_parameter(numeratorCandCatName="SignalDstTvLepMisid",denominatorCandCatName="SignalDstTv",)
    pdfCombiner.define_yield_ratio_parameter(numeratorCandCatName="NormDstlvLepMisid",denominatorCandCatName="NormDstlv",)
    pdfCombiner.define_yield_ratio_parameter(numeratorCandCatName="SignalDTvLepMisid",denominatorCandCatName="SignalDTv",)
    pdfCombiner.define_yield_ratio_parameter(numeratorCandCatName="NormDlvLepMisid",denominatorCandCatName="NormDlv",)

    pdfCombiner.create_fit_model(dataSetLabel=dataSetLabel)
    pdfCombiner.overwrite_model_number()
    pdfCombiner.output_objects(outROOTName=outROOTName)

    return


if __name__ == "__main__":

    args = arg_parser()

    in_args = {}
    in_args["inROOTName"] = args.inputROOT
    in_args["outROOTName"] = args.outputROOT
    in_args["fitLumi"] = args.fitLumi
    in_args["x1VarName"] = args.x1VarName
    in_args["x2VarName"] = args.x2VarName
    in_args["wsName"] = args.wsName
    in_args["dataSetLabel"] = args.dataSetLabel
    in_args["modelNum"] = args.modelNum
    in_args["inputUnfoldFactorJsonName_R"] = args.inputUnfoldFactorJsonName_R
    in_args["_inputFactorJsonName"] = args.inputFactorJsonName

    combine_pdfs(**in_args)
