#include <iostream>


#include "FrameClass.h"
#include "Core.h"

int main() {
    VO::Core app;
    while (app.spin()){

    }
    app.cleanUp();

    return 0;
}
