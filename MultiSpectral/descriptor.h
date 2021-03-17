#ifndef _DESCRIPTOR_H_
#define _DESCRIPTOR_H_



#include "Coregistration.h"

class Descriptor {
  
public:
  
  Descriptor(uint8_t* I,int32_t width,int32_t height,int32_t bpl,bool half_resolution);
  
  ~Descriptor();
  
  uint8_t* I_desc;
  
private:

  void createDescriptor(uint8_t* I_du,uint8_t* I_dv,int32_t width,int32_t height,int32_t bpl,bool half_resolution);

};


#endif