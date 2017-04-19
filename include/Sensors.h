#ifndef SENSORS_DEFD
#define SENSORS_DEFD

#include"Requisites.h"
#include"Sunmodel.h"
void SunSensor(type t)
{
  sunvector_ECI(t/(365*24*60*60)); 
}

#endif
