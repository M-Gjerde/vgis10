/**
 * @file: MultiSense-Viewer/src/Tools/Logger.cpp
 *
 * Copyright 2022
 * Carnegie Robotics, LLC
 * 4501 Hatfield Street, Pittsburgh, PA 15201
 * http://www.carnegierobotics.com
 *
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *     * Neither the name of the Carnegie Robotics, LLC nor the
 *       names of its contributors may be used to endorse or promote products
 *       derived from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL CARNEGIE ROBOTICS, LLC BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * Significant history (date, user, action):
 *   2022-09-12, mgjerde@carnegierobotics.com, Created file.
 **/

// C++ Header File(s)
#include <iostream>
#include <ctime>
#include <vector>

#ifdef WIN32
#define semPost(x) SetEvent(x)
#define semWait(x, y) WaitForSingleObject(x, y)
#else
#include<bits/stdc++.h>
#define semPost(x) sem_post(x)
#endif

#include "Logger.h"
#include "ThreadPool.h"

// Code Specific Header Files(s)
namespace Log {


    Logger *Logger::m_Instance = nullptr;
    ThreadPool *Logger::m_ThreadPool = nullptr;

// Log file m_Name. File m_Name should be change from here only
    const std::string logFileName = "logger.log";

    Logger::Logger() {
        m_File.open(logFileName.c_str(), std::ios::out | std::ios::app);
        m_LogLevel = LOG_LEVEL_TRACE;
        m_LogType = CONSOLE;

    }


    Logger::~Logger() {
        m_File.close();
        delete m_Instance;
        delete m_ThreadPool;
    }

    Logger *Logger::getInstance() noexcept {
        if (m_Instance == nullptr) {
            m_Instance = new Logger();
            m_ThreadPool = new ThreadPool(1);
        }
        return m_Instance;
    }

    void Logger::logIntoFile(void* ctx, std::string &data) {
        auto * app = static_cast<Logger*> (ctx);
        app->m_Mutex.lock();
        app->m_File << getCurrentTime() << "  " << data << std::endl;
        app->m_Mutex.unlock();
    }

    void Logger::logOnConsole(std::string &data) {
        std::cout << getCurrentTime() << "  " << data << std::endl;
    }

    std::string Logger::getCurrentTime() {
        std::string currTime;
        //Current date/time based on current time
        time_t now = time(0);
        // Convert current time to string
        currTime.assign(ctime(&now));

        // Last charactor of currentTime is "\n", so remove it
        std::string currentTime = currTime.substr(0, currTime.size() - 1);
        return currentTime;
    }

// Interface for Error Log
    void Logger::errorInternal(const char *text) noexcept {
        std::string data;
        data.append("[ERROR]: ");
        data.append(text);

        // ERROR must be capture
        if (m_LogType == FILE_LOG) {
            m_ThreadPool->Push(Logger::logIntoFile, this, data);
        } else if (m_LogType == CONSOLE) {
            m_ThreadPool->Push(Logger::logOnConsole, data);
        }
    }


// Interface for Info Log
    void Logger::infoInternal(const char *text) noexcept {
        std::string data;
        data.append("[INFO]: ");
        data.append(text);

        if ((m_LogType == FILE_LOG) && (m_LogLevel >= LOG_LEVEL_INFO)) {
            m_ThreadPool->Push(Logger::logIntoFile, this, data);
        } else if ((m_LogType == CONSOLE) && (m_LogLevel >= LOG_LEVEL_INFO)) {
            m_ThreadPool->Push(Logger::logOnConsole, data);
        }

    }

// Interface for Always Log
    void Logger::always(std::string text) noexcept {
        std::string data;
        data.append("[ALWAYS]: ");
        data.append(text);

        // No check for ALWAYS logs
        if (m_LogType == FILE_LOG) {
            m_ThreadPool->Push(Logger::logIntoFile, this, data);
        } else if (m_LogType == CONSOLE) {
            m_ThreadPool->Push(Logger::logOnConsole, data);
        }
    }

};