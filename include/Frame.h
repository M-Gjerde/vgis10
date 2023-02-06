//
// Created by magnus on 1/28/23.
//

#ifndef VGIS10_FRAME_H
#define VGIS10_FRAME_H

#include <iostream>
#include <Eigen/Core>

#include "Logger.h"


#define PYR_LEVELS 6
#define setting_minGradHistCut 0.5f
#define setting_minGradHistAdd 5.0f
#define setting_gradDownweightPerLevel 0.74f
#define setting_selectDirectionDistribution true
#define setting_gradient_block_32 32

namespace VO {
  struct Frame {
      Frame(){
          //Log::Logger::getInstance()->info("Created Frame Object");
      }
      ~Frame(){
          //Log::Logger::getInstance()->info("Destroyed Frame Object");
      }
      int width = 0;
      int height = 0;
      std::vector<unsigned char> pixels{};
      std::vector<float> dataf{};
      std::vector<Eigen::Vector3f> color{};
      std::vector<std::vector<float>> absSquaredGrad{};
      std::vector<std::vector<Eigen::Vector3f>> pyramid{};

      std::vector<float> histogram{};
      std::vector<float> histogramSmoothed{};
      std::vector<int> gradientHistogram{};

      int widthLevel[PYR_LEVELS]{}, heightLevel[PYR_LEVELS]{};
      float fxLevel[PYR_LEVELS]{}, fyLevel[PYR_LEVELS]{}, cxLevel[PYR_LEVELS]{}, cyLevel[PYR_LEVELS]{};
      float fxiLevel[PYR_LEVELS]{}, fyiLevel[PYR_LEVELS]{}, cxiLevel[PYR_LEVELS]{}, cyiLevel[PYR_LEVELS]{};
      Eigen::Matrix3f KLevel[PYR_LEVELS], KiLevel[PYR_LEVELS], K;

  };
};

#endif //VGIS10_FRAME_H
