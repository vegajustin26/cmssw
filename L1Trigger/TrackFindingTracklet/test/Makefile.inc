#common definitions go in here
CXX = g++ CXXFLAGS = -O - Wall - fPIC - g - ansi - Wextra - Wno - deprecated - std =
                         c++ 11 - I.- I../../..LD = g++ LDFLAGS = -O LIBDIR = $(WORKDIR) AR = ar ARFLAGS = cr

                             OS
    : = $(shell uname - s) ifeq($(OS), Darwin) SOFLAGS = -dynamiclib else SOFLAGS =
          -shared endif

              ifndef ROOTSYS
              $(error *
                **Please set up Root) endif

                  ifeq("$(wildcard lib)", "") $(shell mkdir lib) endif

#
#Add system headers here(hack around problem with makedepend
#if you added file remove the dummyIncl directory to rewrite files
#
                      ifeq("$(wildcard dummyIncl)",
                           "") $(shell mkdir
                                     dummyIncl) $(shell mkdir dummyIncl /
                                                  sys) $(shell mkdir dummyIncl / ext)
                          $(shell touch dummyIncl / TBranch.h) $(shell touch dummyIncl / TCanvas.h) $(
                              shell touch
                                  dummyIncl /
                              TChain
                                  .h) $(shell touch dummyIncl / TColor.h)
                              $(shell touch dummyIncl /
                                TEfficiency.h) $(shell touch dummyIncl / TFile.h)
                                  $(shell touch dummyIncl / TH1D.h) $(shell touch dummyIncl / TH1F.h) $(
                                      shell touch
                                          dummyIncl /
                                      TH2F.h) $(shell touch dummyIncl / TLatex.h)
                                      $(shell touch dummyIncl /
                                        TLegend.h) $(shell touch dummyIncl / TMath.h)
                                          $(shell touch dummyIncl /
                                            TROOT.h) $(shell touch dummyIncl / TStyle.h)
                                              $(shell touch dummyIncl /
                                                TSystem.h) $(shell touch dummyIncl / TTree.h)
                                                  $(shell touch dummyIncl /
                                                    algorithm) $(shell touch dummyIncl / array)
                                                      $(shell touch dummyIncl / bitset) $(shell touch dummyIncl / cassert) $(
                                                          shell
                                                              touch
                                                                  dummyIncl /
                                                          cctype) $(shell touch dummyIncl / cmath)
                                                          $(shell touch dummyIncl /
                                                            cstdlib) $(shell touch dummyIncl / fstream)
                                                              $(shell touch dummyIncl /
                                                                iomanip) $(shell touch dummyIncl / iostream)
                                                                  $(shell touch dummyIncl / limits) $(
                                                                      shell
                                                                          touch
                                                                              dummyIncl /
                                                                      map) $(shell touch dummyIncl /
                                                                             set) $(shell touch dummyIncl / sstream)
                                                                      $(shell touch dummyIncl / string) $(
                                                                          shell touch
                                                                              dummyIncl /
                                                                          vector) $(shell touch dummyIncl / algorithm)
                                                                          $(shell touch dummyIncl / ext / hash_set) $(
                                                                              shell touch dummyIncl /
                                                                              sys /
                                                                              time.h) endif

                                                                              ifeq("$(wildcard FWCore)",
                                                                                   "") $(shell mkdir FWCore)
                                                                                  $(shell mkdir FWCore /
                                                                                    MessageLogger) $(shell mkdir FWCore /
                                                                                                     MessageLogger /
                                                                                                     interface) endif

                                                                                      ifeq(
                                                                                          "$(wildcard "
                                                                                          "FWCore/MessageLogger/"
                                                                                          "interface/MessageLogger.h)",
                                                                                          "") $(shell echo
                                                                                                "#ifndef "
                                                                                                "L1Trigger_"
                                                                                                "TrackFindingTracklet_"
                                                                                                "interface_Logger_h" >
                                                                                                FWCore / MessageLogger /
                                                                                                    interface /
                                                                                                    MessageLogger.h)
                                                                                          $(shell echo
                                                                                            "#define "
                                                                                            "L1Trigger_"
                                                                                            "TrackFindingTracklet_"
                                                                                            "interface_Logger_h" >>
                                                                                            FWCore / MessageLogger /
                                                                                                interface /
                                                                                                MessageLogger
                                                                                                    .h) $(shell echo
                                                                                                          "#include<"
                                                                                                          "iostream>" >>
                                                                                                          FWCore /
                                                                                                              MessageLogger /
                                                                                                              interface /
                                                                                                              MessageLogger
                                                                                                                  .h) $(shell echo
                                                                                                                        "namespace edm {" >>
                                                                                                                        FWCore /
                                                                                                                            MessageLogger / interface /
                                                                                                                            MessageLogger
                                                                                                                                .h)
#INFO statement
                                                                                              $(shell echo
                                                                                                "  class "
                                                                                                "LogVerbatim{" >>
                                                                                                FWCore / MessageLogger /
                                                                                                    interface /
                                                                                                    MessageLogger.h) $(shell echo
                                                                                                                       "  public:" >>
                                                                                                                       FWCore /
                                                                                                                           MessageLogger /
                                                                                                                           interface /
                                                                                                                           MessageLogger
                                                                                                                               .h)
                                                                                                  $(shell echo
                                                                                                    "    "
                                                                                                    "LogVerbatim(std::"
                                                                                                    "string type) { "
                                                                                                    "if(0) "
                                                                                                    "{std::cout<<type<<"
                                                                                                    "std::endl;} }" >>
                                                                                                    FWCore /
                                                                                                        MessageLogger /
                                                                                                        interface /
                                                                                                        MessageLogger.h)
                                                                                                      $(shell echo
                                                                                                        "    "
                                                                                                        "~LogVerbatim"
                                                                                                        "() { "
                                                                                                        "std::cout "
                                                                                                        "<< "
                                                                                                        "std::endl;"
                                                                                                        "}" >>
                                                                                                        FWCore /
                                                                                                            MessageLogger /
                                                                                                            interface /
                                                                                                            MessageLogger
                                                                                                                .h)
                                                                                                          $(shell echo
                                                                                                            "    "
                                                                                                            "template "
                                                                                                            "<class T> "
                                                                                                            "LogVerbati"
                                                                                                            "m& "
                                                                                                            "operator<<"
                                                                                                            "(T const& "
                                                                                                            "t) "
                                                                                                            "{std::"
                                                                                                            "cout << "
                                                                                                            "t;return "
                                                                                                            "*this;}" >>
                                                                                                            FWCore /
                                                                                                                MessageLogger /
                                                                                                                interface /
                                                                                                                MessageLogger
                                                                                                                    .h)
                                                                                                              $(shell echo
                                                                                                                "    "
                                                                                                                "LogVer"
                                                                                                                "batim&"
                                                                                                                " opera"
                                                                                                                "tor<<("
                                                                                                                "std::"
                                                                                                                "ostrea"
                                                                                                                "m& "
                                                                                                                "(*f)("
                                                                                                                "std::"
                                                                                                                "ostrea"
                                                                                                                "m&)) "
                                                                                                                "{std::"
                                                                                                                "cout "
                                                                                                                "<< "
                                                                                                                "f;"
                                                                                                                "return"
                                                                                                                " *"
                                                                                                                "this;"
                                                                                                                "}" >>
                                                                                                                FWCore /
                                                                                                                    MessageLogger /
                                                                                                                    interface /
                                                                                                                    MessageLogger
                                                                                                                        .h) $(shell echo
                                                                                                                              "    "
                                                                                                                              "LogVer"
                                                                                                                              "batim&"
                                                                                                                              " opera"
                                                                                                                              "tor<<("
                                                                                                                              "std::"
                                                                                                                              "ios_"
                                                                                                                              "base& "
                                                                                                                              "(*f)("
                                                                                                                              "std::"
                                                                                                                              "ios_"
                                                                                                                              "base&)"
                                                                                                                              ") "
                                                                                                                              "{std::"
                                                                                                                              "cout "
                                                                                                                              "<< "
                                                                                                                              "f;"
                                                                                                                              "return"
                                                                                                                              " *"
                                                                                                                              "this;"
                                                                                                                              "}" >>
                                                                                                                              FWCore /
                                                                                                                                  MessageLogger /
                                                                                                                                  interface /
                                                                                                                                  MessageLogger
                                                                                                                                      .h)
                                                                                                                  $(shell echo
                                                                                                                    "  "
                                                                                                                    "}"
                                                                                                                    ";" >>
                                                                                                                    FWCore /
                                                                                                                        MessageLogger /
                                                                                                                        interface /
                                                                                                                        MessageLogger
                                                                                                                            .h)
#WARNING statement
                                                                                                                      $(shell echo
                                                                                                                        "  class LogPrint{" >>
                                                                                                                        FWCore /
                                                                                                                            MessageLogger /
                                                                                                                            interface /
                                                                                                                            MessageLogger
                                                                                                                                .h) $(shell echo
                                                                                                                                      "  public:" >>
                                                                                                                                      FWCore /
                                                                                                                                          MessageLogger /
                                                                                                                                          interface /
                                                                                                                                          MessageLogger
                                                                                                                                              .h)
                                                                                                                          $(shell echo
                                                                                                                            "    LogPrint(std::string type) { if(0) {std::cout<<type<<std::endl;} }" >>
                                                                                                                            FWCore /
                                                                                                                                MessageLogger /
                                                                                                                                interface /
                                                                                                                                MessageLogger
                                                                                                                                    .h)
                                                                                                                              $(shell echo
                                                                                                                                "    ~LogPrint() { std::cout << std::endl;}" >>
                                                                                                                                FWCore /
                                                                                                                                    MessageLogger /
                                                                                                                                    interface /
                                                                                                                                    MessageLogger
                                                                                                                                        .h) $(shell echo
                                                                                                                                              "    template <class T> LogPrint& operator<<(T const& t) {std::cout << t;return *this;}" >>
                                                                                                                                              FWCore /
                                                                                                                                                  MessageLogger /
                                                                                                                                                  interface /
                                                                                                                                                  MessageLogger
                                                                                                                                                      .h)
                                                                                                                                  $(shell echo
                                                                                                                                    "    LogPrint& operator<<(std::ostream& (*f)(std::ostream&)) {std::cout << f;return *this;}" >>
                                                                                                                                    FWCore /
                                                                                                                                        MessageLogger /
                                                                                                                                        interface /
                                                                                                                                        MessageLogger
                                                                                                                                            .h)
                                                                                                                                      $(shell echo
                                                                                                                                        "    LogPrint& operator<<(std::ios_base& (*f)(std::ios_base&)) {std::cout << f;return *this;}" >>
                                                                                                                                        FWCore /
                                                                                                                                            MessageLogger /
                                                                                                                                            interface /
                                                                                                                                            MessageLogger
                                                                                                                                                .h) $(shell echo
                                                                                                                                                      "  };" >>
                                                                                                                                                      FWCore /
                                                                                                                                                          MessageLogger /
                                                                                                                                                          interface /
                                                                                                                                                          MessageLogger
                                                                                                                                                              .h)
#ERROR statement
                                                                                                                                          $(shell echo
                                                                                                                                            "  class LogProblem{" >>
                                                                                                                                            FWCore /
                                                                                                                                                MessageLogger /
                                                                                                                                                interface /
                                                                                                                                                MessageLogger
                                                                                                                                                    .h)
                                                                                                                                              $(shell echo
                                                                                                                                                "  public:" >>
                                                                                                                                                FWCore /
                                                                                                                                                    MessageLogger /
                                                                                                                                                    interface /
                                                                                                                                                    MessageLogger
                                                                                                                                                        .h) $(shell echo
                                                                                                                                                              "    LogProblem(std::string type) { if(0) {std::cout<<type<<std::endl;} }" >>
                                                                                                                                                              FWCore /
                                                                                                                                                                  MessageLogger /
                                                                                                                                                                  interface /
                                                                                                                                                                  MessageLogger
                                                                                                                                                                      .h)
                                                                                                                                                  $(shell echo
                                                                                                                                                    "    ~LogProblem() { std::cout << std::endl;}" >>
                                                                                                                                                    FWCore /
                                                                                                                                                        MessageLogger /
                                                                                                                                                        interface /
                                                                                                                                                        MessageLogger
                                                                                                                                                            .h)
                                                                                                                                                      $(shell echo
                                                                                                                                                        "    template <class T> LogProblem& operator<<(T const& t) {std::cout << t;return *this;}" >>
                                                                                                                                                        FWCore /
                                                                                                                                                            MessageLogger /
                                                                                                                                                            interface /
                                                                                                                                                            MessageLogger
                                                                                                                                                                .h) $(shell echo
                                                                                                                                                                      "    LogProblem& operator<<(std::ostream& (*f)(std::ostream&)) {std::cout << f;return *this;}" >>
                                                                                                                                                                      FWCore /
                                                                                                                                                                          MessageLogger /
                                                                                                                                                                          interface /
                                                                                                                                                                          MessageLogger
                                                                                                                                                                              .h)
                                                                                                                                                          $(shell echo
                                                                                                                                                            "    LogProblem& operator<<(std::ios_base& (*f)(std::ios_base&)) {std::cout << f;return *this;}" >>
                                                                                                                                                            FWCore /
                                                                                                                                                                MessageLogger /
                                                                                                                                                                interface /
                                                                                                                                                                MessageLogger
                                                                                                                                                                    .h)
                                                                                                                                                              $(shell echo
                                                                                                                                                                "  };" >>
                                                                                                                                                                FWCore /
                                                                                                                                                                    MessageLogger /
                                                                                                                                                                    interface /
                                                                                                                                                                    MessageLogger
                                                                                                                                                                        .h)
#
                                                                                                                                                                  $(shell echo
                                                                                                                                                                    "};" >>
                                                                                                                                                                    FWCore /
                                                                                                                                                                        MessageLogger /
                                                                                                                                                                        interface /
                                                                                                                                                                        MessageLogger
                                                                                                                                                                            .h)
                                                                                                                                                                      $(shell echo
                                                                                                                                                                        "#endif" >>
                                                                                                                                                                        FWCore /
                                                                                                                                                                            MessageLogger /
                                                                                                                                                                            interface /
                                                                                                                                                                            MessageLogger
                                                                                                                                                                                .h)
                                                                                                                                                                          endif

#find out about the installed ROOT
                                                                                                                                                                              ROOTCONFIG
    : = $(ROOTSYS) / bin / root - config ROOTCFLAGS : = $(shell $(ROOTCONFIG)-- cflags) ROOTLDFLAGS
    : = $(shell $(ROOTCONFIG)-- ldflags) $(shell $(ROOTCONFIG)-- libs-- nonew-- glibs)

#append to options
          CXXFLAGS += $(ROOTCFLAGS) LDFLAGS += $(ROOTLDFLAGS)
#Needed for ROOT - Tree
#LDFLAGS += -L./ FPGAEvent_cxx.so

              LIBS = -lHtml

                         ifneq($(CMSSW_BASE), ) ROOTFIT_INCLUDE
    : = $(shell cd $(CMSSW_BASE); scram tool info roofitcore | grep INCLUDE = | sed 's|INCLUDE=||') ROOTFIT_LIBDIR
    : = $(shell cd $(CMSSW_BASE); scram tool info roofitcore | grep LIBDIR = | sed 's|LIBDIR=||')

          CINTINCLUDES : = -I$(ROOTFIT_INCLUDE) CXXFLAGS += -I$(ROOTFIT_INCLUDE) LDFLAGS += -L$(ROOTFIT_LIBDIR) endif
