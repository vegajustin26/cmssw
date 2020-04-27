#ifndef L1Trigger_TrackFindingTracklet_interface_Logger_h
#define L1Trigger_TrackFindingTracklet_interface_Logger_h


namespace Trklet {

  class LogVerbatim{
    
  public:
    
    LogVerbatim(std::string type) {
      cout << type;
    }
    
    template <class T>
      LogVerbatim& operator<<(T const& t) {
      cout << t;
      return *this;
    }
    LogVerbatim& operator<<(std::ostream& (*f)(std::ostream&)) {
      cout << f;
      return *this;
    }
    LogVerbatim& operator<<(std::ios_base& (*f)(std::ios_base&)) {
      cout << f;
      return *this;
    }    
  }; //end class LogVerbatim
  
};

#endif
