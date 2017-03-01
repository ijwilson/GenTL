#include "StepStone.H"
#include <iostream>

std::ostream &StepStone::print(std::ostream &o) const {
  if (d.size()>0) {
    for (size_t i=0;i<d.size()-1;i++) o << d[i] << "-";
    o << d[d.size()-1];
  }
  return o;
}

std::ostream &operator<<(std::ostream &o, const StepStone &a)
{
  return a.print(o);
}
