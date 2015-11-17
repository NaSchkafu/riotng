#include "adaptintval.h"

void Adapt::setup()
  {
#ifdef FILIB_VERSION
      static bool called = false;
      if( not called ) 
      {
	  filib::fp_traits<double>::setup();
	  called = true;
      }
#endif
  }
