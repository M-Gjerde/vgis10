//
// Created by magnus on 1/28/23.
//

#ifndef VGIS10_FRAME_H
#define VGIS10_FRAME_H

#include <iostream>

namespace VO {
  struct Frame {

      Frame(){
        std::cout << "Construct\n";
      }

      ~Frame(){
        std::cout << "Destruct\n";
      }

      int width = 0;
      int height = 0;
      float pixels = 0;

  };
};

#endif //VGIS10_FRAME_H
