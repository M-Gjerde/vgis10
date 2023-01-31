//
// Created by magnus on 1/28/23.
//

#include <iostream>
#include "Core.h"

namespace VO {
    Core::Core() {
        m_FrameClass = std::make_unique<FrameClass>("/home/magnus/CLionProjects/data/tum_mono_dataset/sequence_02/images");
        m_IsRunning = true;
    }

    bool Core::spin() {

        {
            std::shared_ptr<Frame> frame = m_FrameClass->getNextFrame();
            std::cout << "Use count: " << frame.use_count() << "\n";
        }


        return false;
    }

    void Core::cleanUp() {

    }


}