#ifndef GPU_TIMER_H
#define GPU_TIMER_H

#include <GLFW/glfw3.h>
#include <glad/glad.h>
#include <string>
#include <iostream>
#include <vector>
#include <list>

// offer public interface for GPUTimer and AsyncGPUTimer
class GPUTimer {
public:
    GPUTimer(const std::string& name, bool print)
        : mName(name), mPrint(print), mEnabled(true), mLastTimeMs(0.0) {}
    virtual ~GPUTimer() {}

    // put them before and after the GPU function code block you want to measure
    virtual void Start() = 0;
    virtual void Stop() = 0;

    double GetLastTimeMs() const { return mLastTimeMs; }

    void SetPrint(bool print) { mPrint = print; }
    bool GetPrint() const { return mPrint; }
    void SetEnabled(bool enabled) { mEnabled = enabled; }
    bool GetEnabled() const { return mEnabled; }
    std::string GetName() const { return mName; }

    float getAverageTimeMS() {
        float sum = 0.0;
        for (auto &time : frameTime_list) {
            sum += time;
        }
        return sum / frameTime_list.size();
    }

protected:
    std::string mName;   
    bool mPrint = false;      
    bool mEnabled = true;    
    double mLastTimeMs; 
    std::list<float> frameTime_list; // store the frame time in the sliding window, for average time calculation
    const int num_frames_to_average = 100;
    int num_frames_in_sliding_window = 0;
};

// Async GPU Timer with fixed delay frames
class AsyncGPUTimer : public GPUTimer {
public:
    // if delay frame = 0, it is a SyncGpuTimer which immediately get the result at Stop()
    AsyncGPUTimer(const std::string& name, bool print, int delayFrames = 2)
        : GPUTimer(name, print),
          mDelayFrames(delayFrames),
          mBufferSize(delayFrames+1),
          mCurrentIndex(0)
    {
        mStartQueries.resize(mBufferSize, 0);
        mEndQueries.resize(mBufferSize, 0);
        glGenQueries(mBufferSize, mStartQueries.data());
        glGenQueries(mBufferSize, mEndQueries.data());
        framesRemaining = mDelayFrames;
    }

    ~AsyncGPUTimer() override {
        glDeleteQueries(mBufferSize, mStartQueries.data());
        glDeleteQueries(mBufferSize, mEndQueries.data());
    }

    void Start() override {
        if(!mEnabled){
            return;
        }
        glQueryCounter(mStartQueries[mCurrentIndex], GL_TIMESTAMP);
    }

    void Stop() override {
        if(!mEnabled){
            return;
        }
        glQueryCounter(mEndQueries[mCurrentIndex], GL_TIMESTAMP);
        // check whether we have enough async results yet
        if (framesRemaining > 0) {
            framesRemaining--;
            mCurrentIndex = (mCurrentIndex + 1) % mBufferSize;
            return;
        }

        int readIndex = (mCurrentIndex + mBufferSize - mDelayFrames) % mBufferSize;
        GLint available = 0;
        const int safety = 200000;
        int i = 0;
        do {
            glGetQueryObjectiv(mEndQueries[readIndex], GL_QUERY_RESULT_AVAILABLE, &available);
            i++;
        } while (!available && i < safety);
        if (!available) {
            if (mPrint)
                std::cout << "{" << mName << "}" << ": GPU query result not available after " << safety << " iterations." << std::endl;
        }
    
        GLuint64 startTime = 0, endTime = 0;
        glGetQueryObjectui64v(mStartQueries[readIndex], GL_QUERY_RESULT, &startTime);
        glGetQueryObjectui64v(mEndQueries[readIndex], GL_QUERY_RESULT, &endTime);
        mLastTimeMs = static_cast<double>(endTime - startTime) / 1e6; // convert nanoseconds to milliseconds
    
        if (mPrint) {
            std::cout << "{" << mName << "}" << ": " << mLastTimeMs << " ms" << std::endl;
        }
    
        mCurrentIndex = (mCurrentIndex + 1) % mBufferSize;

        // store the average time in the sliding window
        if (num_frames_in_sliding_window >= num_frames_to_average) {
            frameTime_list.pop_front();
            frameTime_list.push_back(mLastTimeMs);
        }
        else {
            frameTime_list.push_back(mLastTimeMs);
            num_frames_in_sliding_window++;
        }
    }


    

private:
    int mDelayFrames;    // the fixed delay frame number
    int mBufferSize;     // the size of the buffer, equal to mDelayFrames+1
    int mCurrentIndex;   // current index in the buffer

    std::vector<GLuint> mStartQueries; 
    std::vector<GLuint> mEndQueries;  
    // since it has delay frames, we need to ensure the readIndex in Stop() is valid
    int framesRemaining; // =0 means the result is ready to read
};




#endif // GPU_TIMER_H
