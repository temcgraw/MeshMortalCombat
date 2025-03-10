#ifndef FRAME_RATE_MONITOR_H
#define FRAME_RATE_MONITOR_H

#include <chrono>
#include <list>
#include <GLFW/glfw3.h>

class FrameRateMonitor  {
public:
    FrameRateMonitor() {
        reset();
    }


    float getFPS() const {
        return fps;
    }
    float getAverageFPS() const {
        return average_fps;
    }

    bool isRealTime() const {
        return is_realtime;
    }

    float getFrameDeltaTime() const {
        return deltaTime;
    }

    float getFrameDeltaTimeMS() const {
        return deltaTime * 1000.0f;
    }

    void reset() {
        deltaTime = 0.0;
        lastFrame = 0.0;
        LastTime = 0.0;
        is_realtime = true;
        average_fps = 0.0;
        fps = 0.0;
        sliding_deltaTime = 0.0;
        num_frames_in_sliding_window = 0;
        frameTime_list.clear();
    }

    void update() {
        // per-frame time logic
        // --------------------
        float currentFrame = static_cast<float>(glfwGetTime());
        deltaTime = currentFrame - lastFrame;
        lastFrame = currentFrame;
        fps = 1.0f / deltaTime;
        LastTime += deltaTime;



        if (num_frames_in_sliding_window >= num_frames_to_average) {
            frameTime_list.pop_front();
            frameTime_list.push_back(currentFrame);
        }
        else {
            frameTime_list.push_back(fps);
            num_frames_in_sliding_window++;
        }
        average_fps = 1.0f * (num_frames_in_sliding_window - 1) / (frameTime_list.back() - frameTime_list.front());
        // based on the fps, we define whether the simulation is in real time or not
        // if the average fps is less than 24, we consider it as not real time
        is_realtime = fps > 24.0f;
    }

private:
    // timing
    float deltaTime = 0.0f;	// time between current frame and last frame
    float lastFrame = 0.0f;
    float LastTime = 0.0f;
    bool is_realtime = true; // the current simulation is in real time or not
    float average_fps = 0.0f;
    float fps = 0.0f;
    float sliding_deltaTime = 0.0f;
    const int num_frames_to_average = 100;
    int num_frames_in_sliding_window = 0;
    std::list<float> frameTime_list;
};





















#endif