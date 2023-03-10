//
// Created by magnus on 1/28/23.
//

#ifndef VGIS10_CORE_H
#define VGIS10_CORE_H

#include "FrameClass.h"
#include "CameraCalibration.h"
#include "Tracker.h"

namespace VO {
  class Core {
  public:
      Core();

      [[nodiscard]] bool spin();

      void cleanUp();

  private:
      bool m_IsRunning = false;

      std::unique_ptr<FrameClass> m_FrameClass;
      std::unique_ptr<CameraCalibration> m_CamCal;
      std::unique_ptr<Tracker> m_Tracker;

      cv::Mat posImage;
  };
};

#endif //VGIS10_CORE_H
